#include "photogrammetry.h"

// 注：外方位元素以【线元素、角元素】顺序排列

Eigen::IOFormat
fullFormat(		//全精度显示
	Eigen::FullPrecision, 0, ",\t", ";\n", "", "", "[", "]"
);


//单像空间后方交会（只支持φ-ω-κ转角系统）
double Resection(				//返回单位权中误差【单位：mm】
	const Matrix3Xd& groundPts_in_meter			//地面控制点坐标【单位：m】
	, const Matrix2Xd& photoPts_in_millimeter	//控制点对应的像点坐标【单位：mm】
	, const Vector3d innerElement_in_millimeter	//内方元素[x0,y0,f] 【单位：mm】
	, VectorXd& outerElements					//输出外方元素[线元素，角元素] 【单位：m 和 rad】
	, const double m							//航摄设计比例尺分母
	, bool printResult							//在屏幕打印结果
	, MatrixXd* in_weightMatrix					//输入权阵，默认单位阵
	, Matrix3d* out_rotMatrix					//输出旋转矩阵
	, MatrixXd* out_A_Matrix					//输出系数矩阵
	, MatrixXd* out_V_Matrix					//输出改正数矩阵
	, double anglesLimit						//外方位角元素限差【单位：rad】
	, unsigned int maxIterTimes					//最大迭代次数
)
{
	/**测试数据
	内方元素：
		f = 153.24	mm;	x0 = 0.0 mm,y0 = 0.0 mm
	地面坐标：X, Y, Z(m)
		36589.41		  25273.32  		  2195.17
		37631.08		  31324.51		  	  728.69
		39100.97		  24934.98  		  2386.50
		40426.54		  30319.81  		  757.31
	像点坐标：x, y(mm)
		-86.15		 -68.99
		-53.40	  	  82.21
		-14.78       -76.63
		10.46		  64.434
	*/

	//**对所输入的数据进行处理
	int num = groundPts_in_meter.cols();
	if (num != photoPts_in_millimeter.cols()) {
		cout << "ERROR\t控制点坐标数与对应像点数不等！(func: " << __func__ << ")" << endl;
		return -1.0;
	}
	if (num < 4) {
		cout << "ERROR\t单片后方交会至少需要4个控制点！(func: " << __func__ << ")" << endl;
		return -1.0;
	}
	double f = innerElement_in_millimeter(2) / 1000.0;
	const Matrix3Xd &ground = groundPts_in_meter;
	Matrix2Xd P;	P.resize(2, num);
	for (size_t i = 0; i < num; i++) {			//单位转化：mm -> m；并计算以像主点为原点的坐标
		P(0,i) = (photoPts_in_millimeter(0,i) - innerElement_in_millimeter(0)) / 1000.0;
		P(1,i) = (photoPts_in_millimeter(1,i) - innerElement_in_millimeter(1)) / 1000.0;
	}
	
	bool useWieghtMat = true;
	if (!in_weightMatrix ||
		(in_weightMatrix && (in_weightMatrix->rows() != 2 * num || in_weightMatrix->cols() != 2 * num))
		)
		useWieghtMat = false;

	//定义进行计算的变量
	BUNDLE model;
	model.f = f, model.x0 = model.y0 = 0.0;
	MatrixXd X;		X.resize(6, 1);			// 外方位元素的改正数，是需要求解的矩阵
	MatrixXd A;		A.resize(2 * num, 6);	// 误差方程中的偏导数
	MatrixXd L;		L.resize(2 * num, 1);	// 误差方程式中的常数项
	
	//设置外方元素的初值
	model.Xs = ground.row(0).mean();
	model.Ys = ground.row(1).mean();
	model.Zs = ground.row(2).mean();
	//认为航摄比例尺分母都是大于等于500的，不满足则Z方向上直接取均值
	model.Zs += f * (m >= 500.0 ? m : 0);
	model.phi = 0.0, model.omega = 0.0, model.kappa = 0.1;
	if (printResult) {
		printf_s("初值：------\n\t[Xs, Ys, Zs] = [%10lf, %10lf, %10lf] (m)\n", model.Xs, model.Ys, model.Zs);
		printf_s("\t[φ, ω, κ] = [%10lf, %10lf, %10lf] (rad)\n", model.phi, model.omega, model.kappa);
	}
	//迭代计算
	unsigned int t = 0;
	MatrixXd tmpA, tmpL;
	do {
		t++;
		//	计算旋转矩阵R
		model.caculateR();
		for (int i = 0; i < num; ++i) {
			//逐点计算像点坐标近似值 (approx,approy) 和误差方程式中的系数和常数项并组成法方程
			//	计算误差方程式的系数和常数项并组成法方程式
			model.photoPt(0) = P(0, i), model.photoPt(1) = P(1, i);
			model.groundPt(0) = ground(0, i), model.groundPt(1) = ground(1, i)
				, model.groundPt(2) = ground(2, i);
			//以下三行的顺序不可变
			model.caculateXYZ_();
			model.caculteApproXY();
			model.caculate_A_and_L(tmpA, tmpL, 0);
			A.block(2 * i, 0, 2, 6) = tmpA;
			L.block(2 * i, 0, 2, 1) = tmpL;
		}

		// 按照法方程式解的表达式求外方位元素改正值 Xs，Ys，Zs，φ，ω，κ
		if (useWieghtMat)
			X = (A.transpose()*(*in_weightMatrix) * A).inverse() * A.transpose()*(*in_weightMatrix)*L;
		else
			X = (A.transpose() * A).inverse() * A.transpose()*L;

		// 修正外方位元素的值
		model.Xs += X(0, 0), model.Ys += X(1, 0), model.Zs += X(2, 0);
		model.phi += X(3, 0), model.omega += X(4, 0), model.kappa += X(5, 0);
		if (printResult) {
			cout << "\n第" << t << "次迭代\n------------------------------------------" << endl;
			cout << "derta:\t" << X.transpose() << "\n" << "外方角元素:\t"
				<< model.Xs << "\t" << model.Ys << "\t" << model.Zs << "\t"
				<< model.phi << "\t" << model.omega << "\t" << model.kappa << endl;
		}
		// 检查所得改正数是否小于限值与迭代次数检验
	} while ((abs(X(3, 0)) > anglesLimit  ||  abs(X(4, 0)) > anglesLimit || abs(X(5, 0)) > anglesLimit) && t < maxIterTimes);


	//**输出
	if (t >= maxIterTimes) {
		cout << "WARNING\t达到最大迭代次数！( " << __func__ << ")" << endl;
		return -1.0;
	}
	outerElements.resize(6);
	outerElements(0) = model.Xs, outerElements(1) =  model.Ys, outerElements(2) =    model.Zs;
	outerElements(3) = model.phi, outerElements(4) = model.omega, outerElements(5) = model.kappa;
	MatrixXd V = A * X - L;
	if (out_A_Matrix)	*out_A_Matrix = A;
	if (out_V_Matrix)	*out_V_Matrix = V;

	if (in_weightMatrix)	
		V = V.transpose() * (*in_weightMatrix)* V;
	else
		V = V.transpose() * V;

	if (out_rotMatrix)	 model.caculateR(), *out_rotMatrix = model.R;
	
	double m0 = sqrt(V(0) / (2 * num - 6));
	m0 *= 1e6;	//米 -> 微米
	if (printResult) {
		cout << "\n\n旋转矩阵：\n" << model.R.format(fullFormat) << endl;
		printf_s("[Xs, Ys, Zs] = [%10lf, %10lf, %10lf] (m)\n", outerElements(0), outerElements(1), outerElements(2));
		printf_s("[φ, ω, κ] = [%10lf, %10lf, %10lf] (rad)\n", outerElements(3), outerElements(4), outerElements(5));
		cout << "单位权中误差:\t" << m0 << " 微米" << endl;
	}
	return m0;
}

//双像前方交会【双像点投影系数法】
Vector3d pointPrjIntersection(				//返回地面点坐标
	const Vector3d innerElement				//内方元素[x0,y0,f] 【单位：mm】
	, const VectorXd& outerElements1		//片1外方位元素【单位：m 和rad】
	, const Vector2d& pt1					//片1上的对应像点坐标【单位：mm】
	, const VectorXd& outerElements2		//片2外方位元素【单位：m 和rad】
	, const Vector2d& pt2					//片2上的对应像点坐标【单位：mm】
) 
{
	Matrix3d R1 = Rotation::rotate(outerElements1(3), outerElements1(4), outerElements1(5))
		, R2 = Rotation::rotate(outerElements2(3), outerElements2(4), outerElements2(5));
	Matrix<double, 3, 1> PT1, PT2, B;
	//由方位角元素计算相空间辅助坐标
	PT1 << (pt1.x() - innerElement.x()) / 1000.0, (pt1.y() - innerElement.y()) / 1000.0, -innerElement(2) / 1000.0;
	PT2 << (pt2.x() - innerElement.x()) / 1000.0, (pt2.y() - innerElement.y()) / 1000.0, -innerElement(2) / 1000.0;
	//cout << "点投影像空间坐标：\n" << "\t" << PT1.transpose() << "\n\t" << PT2.transpose() << endl;
	PT1 = R1 * PT1, PT2 = R2 * PT2;
	//cout << "旋转矩阵：\n"  << R1 << "\n" << R2 << endl;
	//cout << "点投影像空间辅助坐标：\n" << "\t" << PT1.transpose() << "\n\t" << PT2.transpose() << endl;
	//由外方线元素计算基线分量
	B << (outerElements2(0) - outerElements1(0)), (outerElements2(1) - outerElements1(1)), (outerElements2(2) - outerElements1(2));
	//计算点投影系数
	double N1 = (B(0)*PT2(2) - B(2)*PT2(0)) / (PT1(0)*PT2(2) - PT2(0)*PT1(2))
		, N2 = (B(0)*PT1(2) - B(2)*PT1(0)) / (PT1(0)*PT2(2) - PT2(0)*PT1(2));
	//计算地面点坐标
	Vector3d ground;
	ground.x() = outerElements1(0) + N1 * PT1(0);
	ground.y() = ((outerElements1(1) + N1 * PT1(1) + outerElements2(1) + N2 * PT2(1))) / 2.0;
	
	ground.z() =  outerElements1(2) + N1 * PT1(2);
		//z坐标：//这里改了一下，计算出其摄影测量坐标系下的坐标【可能存在问题：没理由】
		//(outerElements1(2)+ outerElements2(2))/2.0  - (outerElements1(2) + N1 * PT1(2));
	return ground;
}

//多片前方交会-迭代方式
double Intersection(					//返回单位权中误差
	Vector3d &geo_pt					//地面点坐标【单位：m】
	, const Vector3d innerElement		//内方元素[x0,y0,f] 【单位：mm】
	, const Matrix6Xd& outerElements6X	//外方元素		【单位：m、rad】
	, const Matrix2Xd& points2X			//对应像平面坐标【单位：mm】
	, bool printResult					//在屏幕打印结果
	, MatrixXd* in_weightMatrix			//输入权阵，默认单位阵
	, MatrixXd* out_A_Matrix			//输出系数矩阵
	, MatrixXd* out_V_Matrix			//输出改正数矩阵
	, double limit_plane				//平面精度
	, double limit_elevation			//高程精度
	, unsigned int maxIterTimes			//最大迭代次数
)
{

	if (outerElements6X.cols() != points2X.cols()) {
		cout << "ERROR\t外方元素个数与对影像点个数不等！(func: " << __func__ << ")" << endl;
		return -1;
	}
	if (outerElements6X.cols() < 2) {
		cout << "ERROR\t至少需要2对同名像点与对应外方元素！(func: " << __func__ << ")" << endl;
		return -1;
	}
	
	size_t num = outerElements6X.cols();
	//地面点坐标初始值
	VectorXd outE1, outE2;
	MatrixXd tmp = outerElements6X.block(0, 0, 6, 1);
	outE1 = Eigen::Map<Eigen::RowVectorXd>(tmp.data(), tmp.size());
	tmp = outerElements6X.block(0, 1, 6, 1);
	outE2 = Eigen::Map<Eigen::RowVectorXd>(tmp.data(), tmp.size());
	tmp = points2X.block(0, 0, 2, 1);
	Vector2d p1(points2X(0, 0), points2X(1, 0)), p2(points2X(0, 1), points2X(1, 1));
	geo_pt = pointPrjIntersection(innerElement, outE1, p1, outE2, p2);
//	geo_pt(2) = 0.0;
	if(printResult)
		cout << "【点投影法】设置的地面点坐标初始值: " << geo_pt.transpose() << endl;
	//数据准备
	double f = innerElement(2) / 1000.0;
	Matrix2Xd P;	P.resize(2, num);
	for (size_t i = 0; i < num; i++) {			//单位转化：mm -> m；并计算以像主点为原点的坐标
		P(0, i) = (points2X(0, i) - innerElement.x()) / 1000.0;
		P(1, i) = (points2X(1, i) - innerElement.y()) / 1000.0;
	}

	bool useWieghtMat = true;
	if (!in_weightMatrix ||
		(in_weightMatrix && (in_weightMatrix->rows() != 2 * num || in_weightMatrix->cols() != 2 * num))
		)
		useWieghtMat = false;

	//定义进行计算的变量
	BUNDLE model;
	model.f = f, model.x0 = model.y0 = 0.0;
	MatrixXd X;		X.resize(3, 1);			// 地面点坐标改正数，是需要求解的矩阵
	MatrixXd A;		A.resize(2 * num, 3);	// 误差方程中的偏导数
	MatrixXd L;		L.resize(2 * num, 1);	// 误差方程式中的常数项
	//初值
	model.groundPt = geo_pt;
	
	//迭代计算
	unsigned int t = 0;
	MatrixXd tmpA, tmpL;
	double plane_Acu = 100.0;
	do {
		t++;
		for (int i = 0; i < num; ++i) {
			model.photoPt(0) = P(0, i), model.photoPt(1) = P(1, i);

			model.Xs = outerElements6X(0, i), model.Ys = outerElements6X(1, i), model.Zs = outerElements6X(2, i);
			model.phi = outerElements6X(3, i), model.omega = outerElements6X(4, i), model.kappa = outerElements6X(5, i);
			model.caculateR();

			//以下三行的顺序不可变
			model.caculateXYZ_();
			model.caculteApproXY();
			model.caculate_A_and_L(tmpA, tmpL, 3);
			A.block(2 * i, 0, 2, 3) = tmpA;
			L.block(2 * i, 0, 2, 1) = tmpL;
		}

		// 按照法方程式解的表达式求地面点坐标改正数
		if (useWieghtMat)
			X = (A.transpose()*(*in_weightMatrix) * A).inverse() * A.transpose()*(*in_weightMatrix)*L;
		else
			X = (A.transpose() * A).inverse() * A.transpose()*L;

		// 修正地面点坐标
		model.groundPt(0) += X(0, 0), model.groundPt(1) += X(1, 0), model.groundPt(2) += X(2, 0);
		if (printResult) {
			cout << "\n第" << t << "次迭代\n------------------------------------------" << endl;
			cout << "derta:\t\t" << X.transpose().format(fullFormat) << endl
				<< "ground pt:\t" << model.groundPt.transpose().format(fullFormat) << endl;
		}
		// 检查所得改正数是否小于限值与迭代次数检验
		plane_Acu = sqrt(X(0, 0)*X(0, 0) + X(1, 0)*X(1, 0));
	} while ((plane_Acu > limit_plane || abs(X(2, 0)) > limit_elevation) && t < maxIterTimes);

	//**输出
	if (t >= maxIterTimes) {
		printf_s("\nWARNING\t达到最大迭代次数！(fun:%s)\n", __func__);
		cout << "derta:\t\t" << X.transpose().format(fullFormat) << endl
			<< "ground pt:\t" << model.groundPt.transpose().format(fullFormat) << endl;
		return -1;
	}
	geo_pt = model.groundPt;
	MatrixXd V = A * X - L;
	if (out_A_Matrix)	*out_A_Matrix = A;
	if (out_V_Matrix)	*out_V_Matrix = V;
	if (in_weightMatrix)
		V = V.transpose() * (*in_weightMatrix)* V;
	else
		V = V.transpose() * V;

	double m0 = sqrt(V(0) / (2 * num - 3));
	if (printResult) {
		cout << "\nResult: ----------------------------------------" << endl
			<< "[X, Y, Z] = [" << geo_pt.transpose().format(fullFormat) << "]" << endl
			<< "m0:\t" << m0 * 1000.0 << " mm" << endl;
	}
	return m0 * 1000.0;//转为毫米
}

//前方交会-线性方式
double IntersectionLiner(					//返回单位权中误差
	Vector3d &geo_pt					//地面点坐标【单位：m】
	, const Vector3d innerElement		//内方元素[x0,y0,f] 【单位：mm】
	, const Matrix6Xd& outerElements6X	//外方元素		【单位：m、rad】
	, const Matrix2Xd& points2X			//对应像平面坐标【单位：mm】
	, bool printResult					//在屏幕打印结果
	, MatrixXd* out_A_Matrix			//输出系数矩阵
	, MatrixXd* out_V_Matrix			//输出改正数矩阵
)
{
	if (outerElements6X.cols() != points2X.cols()) {
		cout << "ERROR\t外方元素个数与对影像点个数不等！(func: " << __func__ << ")" << endl;
		return -1;
	}
	if (outerElements6X.cols() < 2) {
		cout << "ERROR\t至少需要2对同名像点与对应外方元素！(func: " << __func__ << ")" << endl;
		return -1;
	}

	//数据准备
	size_t num = outerElements6X.cols();
	double f = innerElement(2) / 1000.0;
	Matrix2Xd P;	P.resize(2, num);
	for (size_t i = 0; i < num; i++) {			//单位转化：mm -> m；并计算以像主点为原点的坐标
		P(0, i) = (points2X(0, i) - innerElement.x()) / 1000.0;
		P(1, i) = (points2X(1, i) - innerElement.y()) / 1000.0;
	}

	//定义进行计算的变量
	MatrixXd A;		A.resize(2 * num, 3);	// 误差方程中的偏导数
	MatrixXd L = MatrixXd::Zero(2 * num, 1);// 误差方程式中的常数项
	Matrix3d R;								//旋转矩阵

	//计算系数阵及常数项阵
	for (int k = 0; k < num; k++) {
		R = Rotation::rotate(outerElements6X(3, k), outerElements6X(4, k), outerElements6X(5, k));
		for (int j = 0; j < 2; j++){
			for (int i = 0; i < 3; i++)	{
				A(2 * k + j, i) = f * R(i, j) + P(j, k)*R(i, 2);
				L(2 * k + j, 0) += A(2 * k + j, i)*outerElements6X(i, k);
			}
		}
	}
	MatrixXd X = (A.transpose() * A).inverse() * A.transpose()*L;
	
	geo_pt = matrix2vectorD(X);
	MatrixXd V = A * X - L;
	if (out_A_Matrix)	*out_A_Matrix = A;
	if (out_V_Matrix)	*out_V_Matrix = V;
	V = V.transpose() * V;

	double m0 = sqrt(V(0) / (2 * num - 3));
	if (printResult) {
		cout << "\nResult: ----------------------------------------" << endl
			<< "[X, Y, Z] = [" << geo_pt.transpose().format(fullFormat) << "]" << endl
			<< "m0:\t" << m0 * 1000.0 << " mm" << endl;
	}
	return m0 * 1000.0;//转为毫米
}

//对内定向参数的运用
void internalTrans(
	const VectorXd& prmVec				//内定向参数
	, const Matrix2Xd& pts_in_pixel		//量测的像素坐标【单位：pixel】
	, const double pxl_size				//像素大小【单位：毫末】
	, int type							//内定向方法【0 : 仿射; 1:正形；2:双线性】
	, Matrix2Xd& trgPts_inn_millimeter	//返回实际坐标【单位：mm】
)
{
	MatrixXd xs = pts_in_pixel.row(0) * pxl_size
		, ys = pts_in_pixel.row(1) * pxl_size;
	int num = pts_in_pixel.cols();
	trgPts_inn_millimeter.resizeLike(pts_in_pixel);
	if (type == 0) {	//仿射 x = a0 + a1*x + a2*y
		trgPts_inn_millimeter.block(0, 0, 1, num)
			= (prmVec(0)*MatrixXd::Ones(1, num)) + prmVec(1)*xs + prmVec(2)*ys;
		trgPts_inn_millimeter.block(1, 0, 1, num)
			= (prmVec(3)*MatrixXd::Ones(1, num)) + prmVec(4)*xs + prmVec(5)*ys;
	}
	else if (type == 1) {
		// 正形 [a0, a1, a2, b0]
		//	X = a0 + a1*x - a2*y
		//  Y = b0 + a2*x + a1*y
		trgPts_inn_millimeter.block(0, 0, 1, num)
			= (prmVec(0)*MatrixXd::Ones(1, num)) + prmVec(1)*xs - prmVec(2)*ys;
		trgPts_inn_millimeter.block(1, 0, 1, num)
			= (prmVec(3)*MatrixXd::Ones(1, num)) + prmVec(2)*xs + prmVec(1)*ys;
	}
	else {
		//一般在8框标时用
		trgPts_inn_millimeter.block(0, 0, 1, num)
			= (prmVec(0)*MatrixXd::Ones(1, num)) + prmVec(1)*xs + prmVec(2)*ys
			+ prmVec(3)*multi(xs, ys);
		trgPts_inn_millimeter.block(1, 0, 1, num)
			= (prmVec(4)*MatrixXd::Ones(1, num)) + prmVec(5)*xs + prmVec(6)*ys
			+ prmVec(7)*multi(xs, ys);
	}
}

//内定向
Vector2d internalOritation(					//返回x、y方向上的残差rms【单位：微米】
	const Matrix2Xd &trgPts_in_millimeter	//理论坐标【单位：mm】
	, const Matrix2Xd &pts_in_pixel			//实际量测的坐标【单位：pixel】
	, VectorXd& coefVec						//记录内定向的参数
	, const double pxl_size					//像素大小【单位：mm】
	, Matrix2Xd* RMSs						//每个点上的残差【单位：微米】
	, double* sigma0						//单位权中误差【单位：微米】
	, int type								//内定向方法【0:仿射;1:正形；2:双线性】
	, bool printResult						//打印中间结果
)
{
	if (trgPts_in_millimeter.cols() != pts_in_pixel.cols()) {
		printf_s("ERROR\t内定向理论坐标数与实际坐标数不等！(fun:%s)\n", __func__);
		exit(-1);
	}
	if (trgPts_in_millimeter.cols() < 4) {
		printf_s("ERROR\t内定向坐标对数太少，至少需要4对坐标！(fun:%s)\n", __func__);
		exit(-1);
	}

	int num = trgPts_in_millimeter.cols();
	const Matrix2Xd &tgpts = trgPts_in_millimeter;
	MatrixXd A, L, pts = pts_in_pixel * pxl_size;
	L.resize(2 * num, 1);
	if (type == 0) {	
		//仿射 x = a0 + a1*x + a2*y
		A = MatrixXd::Zero(2 * num, 6);
		for (int i = 0; i < num; i++) {
			A(2*i, 0) = 1.0;
			A(2*i, 1) = pts(0, i);
			A(2*i, 2) = pts(1, i);
			L(2*i, 0) = tgpts(0, i);
			A(2*i + 1, 3) = A(2*i, 0);
			A(2*i + 1, 4) = A(2*i, 1);
			A(2*i + 1, 5) = A(2*i, 2);
			L(2*i + 1, 0) = tgpts(1, i);
		}
		if (printResult)
			cout << "\n--- 内定向方法：仿射变换[a0, a1 ,a2, b0, b1, b2]" << endl;
	}
	else if (type == 1) {
		// 正形 [a0, a1, a2, b0]
		//	X = a0 + a1*x - a2*y
		//  Y = b0 + a2*x + a1*y
		A = MatrixXd::Zero(2 * num, 4);
		for (size_t i = 0; i < num; i++) {
			A(2*i, 0) = 1.0;
			A(2*i, 1) = pts(0, i);
			A(2*i, 2) = - pts(1, i);
			L(2*i, 0) = tgpts(0, i);
			A(2*i + 1, 1) = pts(1, i);
			A(2*i + 1, 2) = pts(0, i);
			A(2*i + 1, 3) = 1.0;
			L(2*i + 1, 0) = tgpts(1, i);
		}
		if (printResult)
			cout << "\n--- 内定向方法：正形变换[a0, a1 ,a2, b0]" << endl;
	}
	else{	
		//双线性 X = a0 + a1*x +a2*y +a3*x*y
		//	该方式一般在量测了8个框标后使用
		A = MatrixXd::Zero(2 * num, 8);
		for (size_t i = 0; i < num; i ++) {
			A(2*i, 0) = 1.0;
			A(2*i, 1) = pts(0, i);
			A(2*i, 2) = pts(1, i);
			A(2*i, 3) = pts(0, i)*pts(1, i);
			L(2*i, 0) = tgpts(0, i);
			A(2*i + 1, 4) = A(2*i, 0);
			A(2*i + 1, 5) = A(2*i, 1);
			A(2*i + 1, 6) = A(2*i, 2);
			A(2*i + 1, 7) = A(2*i, 3);
			L(2*i + 1, 0) = tgpts(1, i);
		}
		if (printResult)
			cout << "\n--- 内定向方法：双线性变换[a0, a1 ,a2, a3, b0, b1, b2, b3]" << endl;
	}

	MatrixXd X = (A.transpose()*A).inverse()*A.transpose()*L
		, V = A * X - L;
	//输出
	coefVec = matrix2vectorD(X);
	
	//单位权中误差【单位：mm】
	double m0 = 1000.0 * sqrt((V.transpose() * V)(0) / (X.rows()));
	if (sigma0)*sigma0 = m0;
	//计算RMS
	Matrix2Xd rms;
	internalTrans(coefVec, pts_in_pixel, pxl_size, type, rms);
	rms -= trgPts_in_millimeter;

	rms *= 1000.0; //毫米-->微米
	if (RMSs)*RMSs = rms;
	A = multi(rms, rms);
	Vector2d rms_xy(A.row(0).sum(), A.row(1).sum());
	rms_xy(0) = sqrt(rms_xy(0) / num), rms_xy(1) = sqrt(rms_xy(1) / num);

	if (printResult) {
		cout << "coefVec:\t" << coefVec.transpose().format(fullFormat) << endl
			<< "改正数:\t" << V.transpose().format(fullFormat) << endl
			<< "单位权中误差：" << m0 << " 微米\t【像素大小:" << pxl_size * 1000.0 << " 微米】" << endl
			<< "\n残差: \t【微米】\n" << rms.transpose().format(fullFormat) << endl
			<< "x_rms: " << rms_xy(0) << "\ty_rms: " << rms_xy(1)
			<< "【单位:微米】" << endl;
	}
	return rms_xy;
}

//对三位相似变换参数进行运用
void similarityTrans3D(
	const VectorXd& prmVec		//七参数向量（X0,Y0,Z0,Lamda,Phi,Omega,Kappa）【单位：m、rad】
	, const Matrix3Xd& relPts	//相对定向后的坐标
	, Matrix3Xd& geoPts			//地面控制点坐标
)
{
	Matrix3d R = Rotation::rotate(prmVec(4), prmVec(5), prmVec(6));
	Matrix<double, 3, 1> tmp;
	tmp << prmVec(0), prmVec(1), prmVec(2);
	geoPts = prmVec(3)*R*relPts;
	for (int i = 0; i < relPts.cols(); i++) {
		geoPts.block(0, i, 3, 1) += tmp;
	}
}

//重心化坐标
void centerPoints(
	MatrixXd& points					//输入、输出坐标【单位一致】
	, VectorXd& center					//输出重心
	, bool colCoordinate				//每一列是一个坐标
	, const VectorXd* centerSetted		//以此为重心进行重心化
)
{
	int dim = 0, num = 0;
	if (colCoordinate) {
		dim = points.rows();
		num = points.cols();
		center.resize(dim);
		if (centerSetted && centerSetted->size() == center.size()) {
			for (int i = 0; i < dim; i++) {
				center(i) = (*centerSetted)(i);
				for (int j = 0; j < num; j++)
					points(i, j) -= center(i);
			}
		}
		else {
			for (int i = 0; i < dim; i++) {
				center(i) = points.row(i).mean();
				for (int j = 0; j < num; j++)
					points(i, j) -= center(i);
			}
		}
	}
	else {
		num = points.rows();
		dim = points.cols();
		center.resize(dim);
		if (centerSetted && centerSetted->size() == center.size()) {
			for (int i = 0; i < dim; i++) {
				center(i) = (*centerSetted)(i);
				for (int j = 0; j < num; j++)
					points(j, i) -= center(i);
			}
		}
		else {
			for (int i = 0; i < dim; i++) {
				center(i) = points.col(i).mean();
				for (int j = 0; j < num; j++)
					points(j, i) -= center(i);
			}
		}
	}
}


//绝对定向
double absoluteOrientation(				//返回绝对定向中误差
	const Matrix3Xd& phtPoints			//绝对定向后的坐标
	, const Matrix3Xd& ctrlPoints		//对应地面点坐标
	, VectorXd& coefVec					//绝对定向元素:七参数向量（X0,Y0,Z0,Lamda,Phi,Omega,Kappa）【单位：m、rad】
	, bool printResult					//在屏幕打印结果
	, MatrixXd* in_weightMatrix			//输入权阵，默认单位阵
	, MatrixXd* out_A_Matrix			//输出系数矩阵
	, MatrixXd* out_V_Matrix			//输出改正数矩阵
	, double limit_angle				//角度精度【单位：rad】
	, unsigned int maxIterTimes			//最大迭代次数
)
{
	if (phtPoints.cols() != ctrlPoints.cols()) {
		cout << "ERROR\t两组点个数不相等！(func: " << __func__ << ")" << endl;
		return -1.0;
	}
	if (ctrlPoints.cols() < 3) {
		cout << "ERROR\t绝对定向至少需要3个点！(func: " << __func__ << ")" << endl;
		return -1.0;
	}

	//权阵
	int num = ctrlPoints.cols();
	bool useWieghtMat = true;
	if (!in_weightMatrix ||
		(in_weightMatrix && (in_weightMatrix->rows() != 3 * num || in_weightMatrix->cols() != 3 * num))
	)
		useWieghtMat = false;

	//重心化
	MatrixXd phtPts = phtPoints, ctlPts = ctrlPoints;
	VectorXd phtCenter, ctlCenter;
	centerPoints(phtPts, phtCenter);
	centerPoints(ctlPts, ctlCenter);

	if (printResult) {
		cout << "\n重心化后：" << endl
			<< "像坐标:\t\tcenter ->" << phtCenter.transpose().format(fullFormat) << endl
			<< phtPts.format(fullFormat) << endl
			<< "对应地面点坐标:\tcenter ->" << ctlCenter.transpose().format(fullFormat) << endl
			<< ctlPts.format(fullFormat) << endl;
	}
	//初始化
	//		X0,Y0,Z0,lamda,Phi,Omega,Kappa
	coefVec.resize(7, 1);
	coefVec << 0., 0., 0., 1., 0., 0., 0.;

	//参数
	MatrixXd A, L, X;

	int t = 0;
	double tmp = 0.0;
	do {
		t++;
		//获取系数和常数项
		similarityTransformationCaculate(phtPts, ctlPts, coefVec, A, L);

		// 求解
		if (useWieghtMat)
			X = (A.transpose()*(*in_weightMatrix) * A).inverse() * A.transpose()*(*in_weightMatrix)*L;
		else
			X = (A.transpose() * A).inverse() * A.transpose()*L;

		// 修正七参数
		coefVec += matrix2vectorD(X);
		if (printResult) {
			cout << "\n第" << t << "次迭代\n------------------------------------------" << endl;
			cout << "derta:  \t" << X.transpose().format(fullFormat) << endl
				<<  "coefVec:\t" << coefVec.transpose().format(fullFormat) << endl;
		}
	} while (
		(abs(X(4)) > limit_angle || abs(X(5)) > limit_angle || abs(X(6)) > limit_angle)
		&& (t < maxIterTimes)
	);

	//反重心化
	Matrix3Xd p;
	similarityTrans3D(coefVec, vector2matrixD(phtCenter), p);
	p = vector2matrixD(ctlCenter) - p;
	coefVec(0) += p(0), coefVec(1) += p(1), coefVec(2) += p(2);

	//**输出
	if (t >= maxIterTimes) {
		printf_s("\nWARNING\t达到最大迭代次数！(fun:%s)\n", __func__);
		cout << "derta:  \t" << X.transpose().format(fullFormat) << endl
			<< "coefVec:\t" << coefVec.transpose().format(fullFormat) << endl;
		return -1;
	}
	MatrixXd V = A * X - L;
	if (out_A_Matrix)	*out_A_Matrix = A;
	if (out_V_Matrix)	*out_V_Matrix = V;
	if (in_weightMatrix)
		V = V.transpose() * (*in_weightMatrix)* V;
	else
		V = V.transpose() * V;

	double m0 = sqrt(V(0) / (3 * num - 7));
	if (printResult) {
		cout << "\nResult: ----------------------------------------" << endl;
		cout << "[X0, Y0, Z0, lamda, phi, omega, kappa] =\n"
			 << coefVec.transpose().format(fullFormat) << endl;
		std::cout << "m0:\t" << m0 << " (Unit is same as ctrlPoints')" << std::endl;
	}
	return m0;
}

//连续法相对定向(双像)
double relativeOrientation(				//返回像点量测中误差【单位：mm】
	const double f						//主距(左右像相同)	【单位：mm】
	, const Matrix2Xd& ptsL				//左像像平面坐标	【单位：mm】
	, const Matrix2Xd& ptsR				//右像像平面坐标	【单位：mm】
	, VectorXd& coefVec					//相对定向元素(phi,omega,kappa,u,v)
	, bool printResult					//在屏幕打印结果
	, MatrixXd* in_weightMatrix			//输入权阵，默认单位阵
	, MatrixXd* out_A_Matrix			//输出系数矩阵
	, MatrixXd* out_V_Matrix			//输出改正数矩阵
	, double limit_angle				//角度精度			【单位：rad】
	, unsigned int maxIterTimes			//最大迭代次数
)
{
	if (ptsL.cols() != ptsR.cols()) {
		cout << "ERROR\t两组像点个数不相等！ (func: " << __func__ << ")" << endl;
		return -1.0;
	}
	if (ptsL.cols() < 5) {
		cout << "ERROR\t相对定向至少需要5组点对！(func: " << __func__ << ")" << endl;
		return -1.0;
	}

	//权阵
	int num = ptsL.cols();
	bool useWieghtMat = true;
	if (!in_weightMatrix || (in_weightMatrix && (in_weightMatrix->rows() != num || in_weightMatrix->cols() != num))) {
		useWieghtMat = false;
	}

	//计算量
	MatrixXd A, L, X, coefMat;
	A.resize(num, 5), L.resize(num, 1), coefMat.resize(5, 1);
	coefMat.fill(0.0);	//初始化
	
	//准备
	double Bx = 0.0, By = 0.0, Bz = 0.0;
	double N1 = 0.0, N2 = 0.0;
	MatrixXd pts1, pts2, rpts2;
	pts1.resize(3, num), pts2.resize(3, num);
	pts1.block(0, 0, 2, num) = ptsL, pts1.row(2).fill(-f);
	pts2.block(0, 0, 2, num) = ptsR, pts2.row(2).fill(-f);

	int t = 0;
	do {
		t++;
		//右像旋转后坐标
		rpts2 = Rotation::rotate(coefMat(0), coefMat(1), coefMat(2)) * pts2;
		for (int i = 0; i < num; i++) {
			//基线计算  ****** X基线设为x2-x1时
			Bx = pts2(0, i) - pts1(0, i);
			By = Bx * coefMat(3), Bz = Bx * coefMat(4);

			//计算N1，N2
			N1 = (Bx*rpts2(2, i) - Bz * rpts2(0, i)) / (pts1(0, i)*rpts2(2, 0) - rpts2(0, 0)*pts1(2, i));
			N2 = (Bx*pts1(2, i) - Bz * pts1(0, i)) / (pts1(0, i)*rpts2(2, 0) - rpts2(0, 0)*pts1(2, i));
			
			//计算误差系数(phi,omega,kappa,u,v)
			A(i,0) = -rpts2(0, i)*rpts2(1, i)*N2 / rpts2(2, i);
			A(i,1) = -(rpts2(2, i) + rpts2(1, i)*rpts2(1, i) / rpts2(2, i))*N2;
			A(i,2) = rpts2(0, i)*N2;
			A(i,3) = Bx;
			A(i,4) = -rpts2(1, i)*Bx / rpts2(2, i);
			//计算上下视差 Q,也即误差方程常数项L
			L(i, 0) = N1 * pts1(1, i) - N2 * rpts2(1, i) - By;
		}
		// 求解
		if (useWieghtMat)
			X = (A.transpose()*(*in_weightMatrix) * A).inverse() * A.transpose()*(*in_weightMatrix)*L;
		else
			X = (A.transpose() * A).inverse() * A.transpose()*L;
		//修正
		coefMat += X;
		if (printResult) {
			cout << "\n第" << t << "次迭代\n------------------------------------------" << endl;
			cout << "derta:  \t" << X.transpose().format(fullFormat) << endl
				<< "coefVec:\t" << coefMat.transpose().format(fullFormat) << endl;
		}
	} while (((abs(X(0)) > limit_angle) || (abs(X(1)) > limit_angle) || (abs(X(2)) > limit_angle)) && (t < maxIterTimes));

	coefVec = matrix2vectorD(coefMat);
	//**输出
	if (t >= maxIterTimes) {
		cout << "\nWARNING\t达到最大迭代次数！(func: " << __func__ << ")" << endl
			<< "derta:  \t" << X.transpose().format(fullFormat) << endl
			<< "coefVec:\t" << coefVec.transpose().format(fullFormat) << endl;
		return -1;
	}
	MatrixXd V = A * X - L;
	if (out_A_Matrix)	*out_A_Matrix = A;
	if (out_V_Matrix)	*out_V_Matrix = V;
	if (in_weightMatrix)
		V = V.transpose() * (*in_weightMatrix)* V;
	else
		V = V.transpose() * V;
	double m0 = sqrt(V(0) / (num - 5));
	if (printResult) {
		cout << "\nResult: ----------------------------------------" << endl;
		cout << "[phi, omega, kappa, u, v] =\n"
			<< coefVec.transpose().format(fullFormat) << endl;
		cout << "m0:\t" << m0 * 1000 << endl;
	}
	return m0*1000.0;
}