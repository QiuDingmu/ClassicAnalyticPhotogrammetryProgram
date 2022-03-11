#include "photogrammetry.h"

// ע���ⷽλԪ���ԡ���Ԫ�ء���Ԫ�ء�˳������

Eigen::IOFormat
fullFormat(		//ȫ������ʾ
	Eigen::FullPrecision, 0, ",\t", ";\n", "", "", "[", "]"
);


//����ռ�󷽽��ᣨֻ֧�֦�-��-��ת��ϵͳ��
double Resection(				//���ص�λȨ������λ��mm��
	const Matrix3Xd& groundPts_in_meter			//������Ƶ����꡾��λ��m��
	, const Matrix2Xd& photoPts_in_millimeter	//���Ƶ��Ӧ��������꡾��λ��mm��
	, const Vector3d innerElement_in_millimeter	//�ڷ�Ԫ��[x0,y0,f] ����λ��mm��
	, VectorXd& outerElements					//����ⷽԪ��[��Ԫ�أ���Ԫ��] ����λ��m �� rad��
	, const double m							//������Ʊ����߷�ĸ
	, bool printResult							//����Ļ��ӡ���
	, MatrixXd* in_weightMatrix					//����Ȩ��Ĭ�ϵ�λ��
	, Matrix3d* out_rotMatrix					//�����ת����
	, MatrixXd* out_A_Matrix					//���ϵ������
	, MatrixXd* out_V_Matrix					//�������������
	, double anglesLimit						//�ⷽλ��Ԫ���޲��λ��rad��
	, unsigned int maxIterTimes					//����������
)
{
	/**��������
	�ڷ�Ԫ�أ�
		f = 153.24	mm;	x0 = 0.0 mm,y0 = 0.0 mm
	�������꣺X, Y, Z(m)
		36589.41		  25273.32  		  2195.17
		37631.08		  31324.51		  	  728.69
		39100.97		  24934.98  		  2386.50
		40426.54		  30319.81  		  757.31
	������꣺x, y(mm)
		-86.15		 -68.99
		-53.40	  	  82.21
		-14.78       -76.63
		10.46		  64.434
	*/

	//**������������ݽ��д���
	int num = groundPts_in_meter.cols();
	if (num != photoPts_in_millimeter.cols()) {
		cout << "ERROR\t���Ƶ����������Ӧ��������ȣ�(func: " << __func__ << ")" << endl;
		return -1.0;
	}
	if (num < 4) {
		cout << "ERROR\t��Ƭ�󷽽���������Ҫ4�����Ƶ㣡(func: " << __func__ << ")" << endl;
		return -1.0;
	}
	double f = innerElement_in_millimeter(2) / 1000.0;
	const Matrix3Xd &ground = groundPts_in_meter;
	Matrix2Xd P;	P.resize(2, num);
	for (size_t i = 0; i < num; i++) {			//��λת����mm -> m����������������Ϊԭ�������
		P(0,i) = (photoPts_in_millimeter(0,i) - innerElement_in_millimeter(0)) / 1000.0;
		P(1,i) = (photoPts_in_millimeter(1,i) - innerElement_in_millimeter(1)) / 1000.0;
	}
	
	bool useWieghtMat = true;
	if (!in_weightMatrix ||
		(in_weightMatrix && (in_weightMatrix->rows() != 2 * num || in_weightMatrix->cols() != 2 * num))
		)
		useWieghtMat = false;

	//������м���ı���
	BUNDLE model;
	model.f = f, model.x0 = model.y0 = 0.0;
	MatrixXd X;		X.resize(6, 1);			// �ⷽλԪ�صĸ�����������Ҫ���ľ���
	MatrixXd A;		A.resize(2 * num, 6);	// �����е�ƫ����
	MatrixXd L;		L.resize(2 * num, 1);	// ����ʽ�еĳ�����
	
	//�����ⷽԪ�صĳ�ֵ
	model.Xs = ground.row(0).mean();
	model.Ys = ground.row(1).mean();
	model.Zs = ground.row(2).mean();
	//��Ϊ��������߷�ĸ���Ǵ��ڵ���500�ģ���������Z������ֱ��ȡ��ֵ
	model.Zs += f * (m >= 500.0 ? m : 0);
	model.phi = 0.0, model.omega = 0.0, model.kappa = 0.1;
	if (printResult) {
		printf_s("��ֵ��------\n\t[Xs, Ys, Zs] = [%10lf, %10lf, %10lf] (m)\n", model.Xs, model.Ys, model.Zs);
		printf_s("\t[��, ��, ��] = [%10lf, %10lf, %10lf] (rad)\n", model.phi, model.omega, model.kappa);
	}
	//��������
	unsigned int t = 0;
	MatrixXd tmpA, tmpL;
	do {
		t++;
		//	������ת����R
		model.caculateR();
		for (int i = 0; i < num; ++i) {
			//����������������ֵ (approx,approy) ������ʽ�е�ϵ���ͳ������ɷ�����
			//	��������ʽ��ϵ���ͳ������ɷ�����ʽ
			model.photoPt(0) = P(0, i), model.photoPt(1) = P(1, i);
			model.groundPt(0) = ground(0, i), model.groundPt(1) = ground(1, i)
				, model.groundPt(2) = ground(2, i);
			//�������е�˳�򲻿ɱ�
			model.caculateXYZ_();
			model.caculteApproXY();
			model.caculate_A_and_L(tmpA, tmpL, 0);
			A.block(2 * i, 0, 2, 6) = tmpA;
			L.block(2 * i, 0, 2, 1) = tmpL;
		}

		// ���շ�����ʽ��ı��ʽ���ⷽλԪ�ظ���ֵ Xs��Ys��Zs���գ��أ���
		if (useWieghtMat)
			X = (A.transpose()*(*in_weightMatrix) * A).inverse() * A.transpose()*(*in_weightMatrix)*L;
		else
			X = (A.transpose() * A).inverse() * A.transpose()*L;

		// �����ⷽλԪ�ص�ֵ
		model.Xs += X(0, 0), model.Ys += X(1, 0), model.Zs += X(2, 0);
		model.phi += X(3, 0), model.omega += X(4, 0), model.kappa += X(5, 0);
		if (printResult) {
			cout << "\n��" << t << "�ε���\n------------------------------------------" << endl;
			cout << "derta:\t" << X.transpose() << "\n" << "�ⷽ��Ԫ��:\t"
				<< model.Xs << "\t" << model.Ys << "\t" << model.Zs << "\t"
				<< model.phi << "\t" << model.omega << "\t" << model.kappa << endl;
		}
		// ������ø������Ƿ�С����ֵ�������������
	} while ((abs(X(3, 0)) > anglesLimit  ||  abs(X(4, 0)) > anglesLimit || abs(X(5, 0)) > anglesLimit) && t < maxIterTimes);


	//**���
	if (t >= maxIterTimes) {
		cout << "WARNING\t�ﵽ������������( " << __func__ << ")" << endl;
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
	m0 *= 1e6;	//�� -> ΢��
	if (printResult) {
		cout << "\n\n��ת����\n" << model.R.format(fullFormat) << endl;
		printf_s("[Xs, Ys, Zs] = [%10lf, %10lf, %10lf] (m)\n", outerElements(0), outerElements(1), outerElements(2));
		printf_s("[��, ��, ��] = [%10lf, %10lf, %10lf] (rad)\n", outerElements(3), outerElements(4), outerElements(5));
		cout << "��λȨ�����:\t" << m0 << " ΢��" << endl;
	}
	return m0;
}

//˫��ǰ�����᡾˫���ͶӰϵ������
Vector3d pointPrjIntersection(				//���ص��������
	const Vector3d innerElement				//�ڷ�Ԫ��[x0,y0,f] ����λ��mm��
	, const VectorXd& outerElements1		//Ƭ1�ⷽλԪ�ء���λ��m ��rad��
	, const Vector2d& pt1					//Ƭ1�ϵĶ�Ӧ������꡾��λ��mm��
	, const VectorXd& outerElements2		//Ƭ2�ⷽλԪ�ء���λ��m ��rad��
	, const Vector2d& pt2					//Ƭ2�ϵĶ�Ӧ������꡾��λ��mm��
) 
{
	Matrix3d R1 = Rotation::rotate(outerElements1(3), outerElements1(4), outerElements1(5))
		, R2 = Rotation::rotate(outerElements2(3), outerElements2(4), outerElements2(5));
	Matrix<double, 3, 1> PT1, PT2, B;
	//�ɷ�λ��Ԫ�ؼ�����ռ丨������
	PT1 << (pt1.x() - innerElement.x()) / 1000.0, (pt1.y() - innerElement.y()) / 1000.0, -innerElement(2) / 1000.0;
	PT2 << (pt2.x() - innerElement.x()) / 1000.0, (pt2.y() - innerElement.y()) / 1000.0, -innerElement(2) / 1000.0;
	//cout << "��ͶӰ��ռ����꣺\n" << "\t" << PT1.transpose() << "\n\t" << PT2.transpose() << endl;
	PT1 = R1 * PT1, PT2 = R2 * PT2;
	//cout << "��ת����\n"  << R1 << "\n" << R2 << endl;
	//cout << "��ͶӰ��ռ丨�����꣺\n" << "\t" << PT1.transpose() << "\n\t" << PT2.transpose() << endl;
	//���ⷽ��Ԫ�ؼ�����߷���
	B << (outerElements2(0) - outerElements1(0)), (outerElements2(1) - outerElements1(1)), (outerElements2(2) - outerElements1(2));
	//�����ͶӰϵ��
	double N1 = (B(0)*PT2(2) - B(2)*PT2(0)) / (PT1(0)*PT2(2) - PT2(0)*PT1(2))
		, N2 = (B(0)*PT1(2) - B(2)*PT1(0)) / (PT1(0)*PT2(2) - PT2(0)*PT1(2));
	//������������
	Vector3d ground;
	ground.x() = outerElements1(0) + N1 * PT1(0);
	ground.y() = ((outerElements1(1) + N1 * PT1(1) + outerElements2(1) + N2 * PT2(1))) / 2.0;
	
	ground.z() =  outerElements1(2) + N1 * PT1(2);
		//z���꣺//�������һ�£����������Ӱ��������ϵ�µ����꡾���ܴ������⣺û���ɡ�
		//(outerElements1(2)+ outerElements2(2))/2.0  - (outerElements1(2) + N1 * PT1(2));
	return ground;
}

//��Ƭǰ������-������ʽ
double Intersection(					//���ص�λȨ�����
	Vector3d &geo_pt					//��������꡾��λ��m��
	, const Vector3d innerElement		//�ڷ�Ԫ��[x0,y0,f] ����λ��mm��
	, const Matrix6Xd& outerElements6X	//�ⷽԪ��		����λ��m��rad��
	, const Matrix2Xd& points2X			//��Ӧ��ƽ�����꡾��λ��mm��
	, bool printResult					//����Ļ��ӡ���
	, MatrixXd* in_weightMatrix			//����Ȩ��Ĭ�ϵ�λ��
	, MatrixXd* out_A_Matrix			//���ϵ������
	, MatrixXd* out_V_Matrix			//�������������
	, double limit_plane				//ƽ�澫��
	, double limit_elevation			//�߳̾���
	, unsigned int maxIterTimes			//����������
)
{

	if (outerElements6X.cols() != points2X.cols()) {
		cout << "ERROR\t�ⷽԪ�ظ������Ӱ���������ȣ�(func: " << __func__ << ")" << endl;
		return -1;
	}
	if (outerElements6X.cols() < 2) {
		cout << "ERROR\t������Ҫ2��ͬ��������Ӧ�ⷽԪ�أ�(func: " << __func__ << ")" << endl;
		return -1;
	}
	
	size_t num = outerElements6X.cols();
	//����������ʼֵ
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
		cout << "����ͶӰ�������õĵ���������ʼֵ: " << geo_pt.transpose() << endl;
	//����׼��
	double f = innerElement(2) / 1000.0;
	Matrix2Xd P;	P.resize(2, num);
	for (size_t i = 0; i < num; i++) {			//��λת����mm -> m����������������Ϊԭ�������
		P(0, i) = (points2X(0, i) - innerElement.x()) / 1000.0;
		P(1, i) = (points2X(1, i) - innerElement.y()) / 1000.0;
	}

	bool useWieghtMat = true;
	if (!in_weightMatrix ||
		(in_weightMatrix && (in_weightMatrix->rows() != 2 * num || in_weightMatrix->cols() != 2 * num))
		)
		useWieghtMat = false;

	//������м���ı���
	BUNDLE model;
	model.f = f, model.x0 = model.y0 = 0.0;
	MatrixXd X;		X.resize(3, 1);			// ��������������������Ҫ���ľ���
	MatrixXd A;		A.resize(2 * num, 3);	// �����е�ƫ����
	MatrixXd L;		L.resize(2 * num, 1);	// ����ʽ�еĳ�����
	//��ֵ
	model.groundPt = geo_pt;
	
	//��������
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

			//�������е�˳�򲻿ɱ�
			model.caculateXYZ_();
			model.caculteApproXY();
			model.caculate_A_and_L(tmpA, tmpL, 3);
			A.block(2 * i, 0, 2, 3) = tmpA;
			L.block(2 * i, 0, 2, 1) = tmpL;
		}

		// ���շ�����ʽ��ı��ʽ���������������
		if (useWieghtMat)
			X = (A.transpose()*(*in_weightMatrix) * A).inverse() * A.transpose()*(*in_weightMatrix)*L;
		else
			X = (A.transpose() * A).inverse() * A.transpose()*L;

		// �������������
		model.groundPt(0) += X(0, 0), model.groundPt(1) += X(1, 0), model.groundPt(2) += X(2, 0);
		if (printResult) {
			cout << "\n��" << t << "�ε���\n------------------------------------------" << endl;
			cout << "derta:\t\t" << X.transpose().format(fullFormat) << endl
				<< "ground pt:\t" << model.groundPt.transpose().format(fullFormat) << endl;
		}
		// ������ø������Ƿ�С����ֵ�������������
		plane_Acu = sqrt(X(0, 0)*X(0, 0) + X(1, 0)*X(1, 0));
	} while ((plane_Acu > limit_plane || abs(X(2, 0)) > limit_elevation) && t < maxIterTimes);

	//**���
	if (t >= maxIterTimes) {
		printf_s("\nWARNING\t�ﵽ������������(fun:%s)\n", __func__);
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
	return m0 * 1000.0;//תΪ����
}

//ǰ������-���Է�ʽ
double IntersectionLiner(					//���ص�λȨ�����
	Vector3d &geo_pt					//��������꡾��λ��m��
	, const Vector3d innerElement		//�ڷ�Ԫ��[x0,y0,f] ����λ��mm��
	, const Matrix6Xd& outerElements6X	//�ⷽԪ��		����λ��m��rad��
	, const Matrix2Xd& points2X			//��Ӧ��ƽ�����꡾��λ��mm��
	, bool printResult					//����Ļ��ӡ���
	, MatrixXd* out_A_Matrix			//���ϵ������
	, MatrixXd* out_V_Matrix			//�������������
)
{
	if (outerElements6X.cols() != points2X.cols()) {
		cout << "ERROR\t�ⷽԪ�ظ������Ӱ���������ȣ�(func: " << __func__ << ")" << endl;
		return -1;
	}
	if (outerElements6X.cols() < 2) {
		cout << "ERROR\t������Ҫ2��ͬ��������Ӧ�ⷽԪ�أ�(func: " << __func__ << ")" << endl;
		return -1;
	}

	//����׼��
	size_t num = outerElements6X.cols();
	double f = innerElement(2) / 1000.0;
	Matrix2Xd P;	P.resize(2, num);
	for (size_t i = 0; i < num; i++) {			//��λת����mm -> m����������������Ϊԭ�������
		P(0, i) = (points2X(0, i) - innerElement.x()) / 1000.0;
		P(1, i) = (points2X(1, i) - innerElement.y()) / 1000.0;
	}

	//������м���ı���
	MatrixXd A;		A.resize(2 * num, 3);	// �����е�ƫ����
	MatrixXd L = MatrixXd::Zero(2 * num, 1);// ����ʽ�еĳ�����
	Matrix3d R;								//��ת����

	//����ϵ���󼰳�������
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
	return m0 * 1000.0;//תΪ����
}

//���ڶ������������
void internalTrans(
	const VectorXd& prmVec				//�ڶ������
	, const Matrix2Xd& pts_in_pixel		//������������꡾��λ��pixel��
	, const double pxl_size				//���ش�С����λ����ĩ��
	, int type							//�ڶ��򷽷���0 : ����; 1:���Σ�2:˫���ԡ�
	, Matrix2Xd& trgPts_inn_millimeter	//����ʵ�����꡾��λ��mm��
)
{
	MatrixXd xs = pts_in_pixel.row(0) * pxl_size
		, ys = pts_in_pixel.row(1) * pxl_size;
	int num = pts_in_pixel.cols();
	trgPts_inn_millimeter.resizeLike(pts_in_pixel);
	if (type == 0) {	//���� x = a0 + a1*x + a2*y
		trgPts_inn_millimeter.block(0, 0, 1, num)
			= (prmVec(0)*MatrixXd::Ones(1, num)) + prmVec(1)*xs + prmVec(2)*ys;
		trgPts_inn_millimeter.block(1, 0, 1, num)
			= (prmVec(3)*MatrixXd::Ones(1, num)) + prmVec(4)*xs + prmVec(5)*ys;
	}
	else if (type == 1) {
		// ���� [a0, a1, a2, b0]
		//	X = a0 + a1*x - a2*y
		//  Y = b0 + a2*x + a1*y
		trgPts_inn_millimeter.block(0, 0, 1, num)
			= (prmVec(0)*MatrixXd::Ones(1, num)) + prmVec(1)*xs - prmVec(2)*ys;
		trgPts_inn_millimeter.block(1, 0, 1, num)
			= (prmVec(3)*MatrixXd::Ones(1, num)) + prmVec(2)*xs + prmVec(1)*ys;
	}
	else {
		//һ����8���ʱ��
		trgPts_inn_millimeter.block(0, 0, 1, num)
			= (prmVec(0)*MatrixXd::Ones(1, num)) + prmVec(1)*xs + prmVec(2)*ys
			+ prmVec(3)*multi(xs, ys);
		trgPts_inn_millimeter.block(1, 0, 1, num)
			= (prmVec(4)*MatrixXd::Ones(1, num)) + prmVec(5)*xs + prmVec(6)*ys
			+ prmVec(7)*multi(xs, ys);
	}
}

//�ڶ���
Vector2d internalOritation(					//����x��y�����ϵĲв�rms����λ��΢�ס�
	const Matrix2Xd &trgPts_in_millimeter	//�������꡾��λ��mm��
	, const Matrix2Xd &pts_in_pixel			//ʵ����������꡾��λ��pixel��
	, VectorXd& coefVec						//��¼�ڶ���Ĳ���
	, const double pxl_size					//���ش�С����λ��mm��
	, Matrix2Xd* RMSs						//ÿ�����ϵĲв��λ��΢�ס�
	, double* sigma0						//��λȨ������λ��΢�ס�
	, int type								//�ڶ��򷽷���0:����;1:���Σ�2:˫���ԡ�
	, bool printResult						//��ӡ�м���
)
{
	if (trgPts_in_millimeter.cols() != pts_in_pixel.cols()) {
		printf_s("ERROR\t�ڶ���������������ʵ�����������ȣ�(fun:%s)\n", __func__);
		exit(-1);
	}
	if (trgPts_in_millimeter.cols() < 4) {
		printf_s("ERROR\t�ڶ����������̫�٣�������Ҫ4�����꣡(fun:%s)\n", __func__);
		exit(-1);
	}

	int num = trgPts_in_millimeter.cols();
	const Matrix2Xd &tgpts = trgPts_in_millimeter;
	MatrixXd A, L, pts = pts_in_pixel * pxl_size;
	L.resize(2 * num, 1);
	if (type == 0) {	
		//���� x = a0 + a1*x + a2*y
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
			cout << "\n--- �ڶ��򷽷�������任[a0, a1 ,a2, b0, b1, b2]" << endl;
	}
	else if (type == 1) {
		// ���� [a0, a1, a2, b0]
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
			cout << "\n--- �ڶ��򷽷������α任[a0, a1 ,a2, b0]" << endl;
	}
	else{	
		//˫���� X = a0 + a1*x +a2*y +a3*x*y
		//	�÷�ʽһ����������8������ʹ��
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
			cout << "\n--- �ڶ��򷽷���˫���Ա任[a0, a1 ,a2, a3, b0, b1, b2, b3]" << endl;
	}

	MatrixXd X = (A.transpose()*A).inverse()*A.transpose()*L
		, V = A * X - L;
	//���
	coefVec = matrix2vectorD(X);
	
	//��λȨ������λ��mm��
	double m0 = 1000.0 * sqrt((V.transpose() * V)(0) / (X.rows()));
	if (sigma0)*sigma0 = m0;
	//����RMS
	Matrix2Xd rms;
	internalTrans(coefVec, pts_in_pixel, pxl_size, type, rms);
	rms -= trgPts_in_millimeter;

	rms *= 1000.0; //����-->΢��
	if (RMSs)*RMSs = rms;
	A = multi(rms, rms);
	Vector2d rms_xy(A.row(0).sum(), A.row(1).sum());
	rms_xy(0) = sqrt(rms_xy(0) / num), rms_xy(1) = sqrt(rms_xy(1) / num);

	if (printResult) {
		cout << "coefVec:\t" << coefVec.transpose().format(fullFormat) << endl
			<< "������:\t" << V.transpose().format(fullFormat) << endl
			<< "��λȨ����" << m0 << " ΢��\t�����ش�С:" << pxl_size * 1000.0 << " ΢�ס�" << endl
			<< "\n�в�: \t��΢�ס�\n" << rms.transpose().format(fullFormat) << endl
			<< "x_rms: " << rms_xy(0) << "\ty_rms: " << rms_xy(1)
			<< "����λ:΢�ס�" << endl;
	}
	return rms_xy;
}

//����λ���Ʊ任������������
void similarityTrans3D(
	const VectorXd& prmVec		//�߲���������X0,Y0,Z0,Lamda,Phi,Omega,Kappa������λ��m��rad��
	, const Matrix3Xd& relPts	//��Զ���������
	, Matrix3Xd& geoPts			//������Ƶ�����
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

//���Ļ�����
void centerPoints(
	MatrixXd& points					//���롢������꡾��λһ�¡�
	, VectorXd& center					//�������
	, bool colCoordinate				//ÿһ����һ������
	, const VectorXd* centerSetted		//�Դ�Ϊ���Ľ������Ļ�
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


//���Զ���
double absoluteOrientation(				//���ؾ��Զ��������
	const Matrix3Xd& phtPoints			//���Զ���������
	, const Matrix3Xd& ctrlPoints		//��Ӧ���������
	, VectorXd& coefVec					//���Զ���Ԫ��:�߲���������X0,Y0,Z0,Lamda,Phi,Omega,Kappa������λ��m��rad��
	, bool printResult					//����Ļ��ӡ���
	, MatrixXd* in_weightMatrix			//����Ȩ��Ĭ�ϵ�λ��
	, MatrixXd* out_A_Matrix			//���ϵ������
	, MatrixXd* out_V_Matrix			//�������������
	, double limit_angle				//�ǶȾ��ȡ���λ��rad��
	, unsigned int maxIterTimes			//����������
)
{
	if (phtPoints.cols() != ctrlPoints.cols()) {
		cout << "ERROR\t������������ȣ�(func: " << __func__ << ")" << endl;
		return -1.0;
	}
	if (ctrlPoints.cols() < 3) {
		cout << "ERROR\t���Զ���������Ҫ3���㣡(func: " << __func__ << ")" << endl;
		return -1.0;
	}

	//Ȩ��
	int num = ctrlPoints.cols();
	bool useWieghtMat = true;
	if (!in_weightMatrix ||
		(in_weightMatrix && (in_weightMatrix->rows() != 3 * num || in_weightMatrix->cols() != 3 * num))
	)
		useWieghtMat = false;

	//���Ļ�
	MatrixXd phtPts = phtPoints, ctlPts = ctrlPoints;
	VectorXd phtCenter, ctlCenter;
	centerPoints(phtPts, phtCenter);
	centerPoints(ctlPts, ctlCenter);

	if (printResult) {
		cout << "\n���Ļ���" << endl
			<< "������:\t\tcenter ->" << phtCenter.transpose().format(fullFormat) << endl
			<< phtPts.format(fullFormat) << endl
			<< "��Ӧ���������:\tcenter ->" << ctlCenter.transpose().format(fullFormat) << endl
			<< ctlPts.format(fullFormat) << endl;
	}
	//��ʼ��
	//		X0,Y0,Z0,lamda,Phi,Omega,Kappa
	coefVec.resize(7, 1);
	coefVec << 0., 0., 0., 1., 0., 0., 0.;

	//����
	MatrixXd A, L, X;

	int t = 0;
	double tmp = 0.0;
	do {
		t++;
		//��ȡϵ���ͳ�����
		similarityTransformationCaculate(phtPts, ctlPts, coefVec, A, L);

		// ���
		if (useWieghtMat)
			X = (A.transpose()*(*in_weightMatrix) * A).inverse() * A.transpose()*(*in_weightMatrix)*L;
		else
			X = (A.transpose() * A).inverse() * A.transpose()*L;

		// �����߲���
		coefVec += matrix2vectorD(X);
		if (printResult) {
			cout << "\n��" << t << "�ε���\n------------------------------------------" << endl;
			cout << "derta:  \t" << X.transpose().format(fullFormat) << endl
				<<  "coefVec:\t" << coefVec.transpose().format(fullFormat) << endl;
		}
	} while (
		(abs(X(4)) > limit_angle || abs(X(5)) > limit_angle || abs(X(6)) > limit_angle)
		&& (t < maxIterTimes)
	);

	//�����Ļ�
	Matrix3Xd p;
	similarityTrans3D(coefVec, vector2matrixD(phtCenter), p);
	p = vector2matrixD(ctlCenter) - p;
	coefVec(0) += p(0), coefVec(1) += p(1), coefVec(2) += p(2);

	//**���
	if (t >= maxIterTimes) {
		printf_s("\nWARNING\t�ﵽ������������(fun:%s)\n", __func__);
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

//��������Զ���(˫��)
double relativeOrientation(				//�����������������λ��mm��
	const double f						//����(��������ͬ)	����λ��mm��
	, const Matrix2Xd& ptsL				//������ƽ������	����λ��mm��
	, const Matrix2Xd& ptsR				//������ƽ������	����λ��mm��
	, VectorXd& coefVec					//��Զ���Ԫ��(phi,omega,kappa,u,v)
	, bool printResult					//����Ļ��ӡ���
	, MatrixXd* in_weightMatrix			//����Ȩ��Ĭ�ϵ�λ��
	, MatrixXd* out_A_Matrix			//���ϵ������
	, MatrixXd* out_V_Matrix			//�������������
	, double limit_angle				//�ǶȾ���			����λ��rad��
	, unsigned int maxIterTimes			//����������
)
{
	if (ptsL.cols() != ptsR.cols()) {
		cout << "ERROR\t��������������ȣ� (func: " << __func__ << ")" << endl;
		return -1.0;
	}
	if (ptsL.cols() < 5) {
		cout << "ERROR\t��Զ���������Ҫ5���ԣ�(func: " << __func__ << ")" << endl;
		return -1.0;
	}

	//Ȩ��
	int num = ptsL.cols();
	bool useWieghtMat = true;
	if (!in_weightMatrix || (in_weightMatrix && (in_weightMatrix->rows() != num || in_weightMatrix->cols() != num))) {
		useWieghtMat = false;
	}

	//������
	MatrixXd A, L, X, coefMat;
	A.resize(num, 5), L.resize(num, 1), coefMat.resize(5, 1);
	coefMat.fill(0.0);	//��ʼ��
	
	//׼��
	double Bx = 0.0, By = 0.0, Bz = 0.0;
	double N1 = 0.0, N2 = 0.0;
	MatrixXd pts1, pts2, rpts2;
	pts1.resize(3, num), pts2.resize(3, num);
	pts1.block(0, 0, 2, num) = ptsL, pts1.row(2).fill(-f);
	pts2.block(0, 0, 2, num) = ptsR, pts2.row(2).fill(-f);

	int t = 0;
	do {
		t++;
		//������ת������
		rpts2 = Rotation::rotate(coefMat(0), coefMat(1), coefMat(2)) * pts2;
		for (int i = 0; i < num; i++) {
			//���߼���  ****** X������Ϊx2-x1ʱ
			Bx = pts2(0, i) - pts1(0, i);
			By = Bx * coefMat(3), Bz = Bx * coefMat(4);

			//����N1��N2
			N1 = (Bx*rpts2(2, i) - Bz * rpts2(0, i)) / (pts1(0, i)*rpts2(2, 0) - rpts2(0, 0)*pts1(2, i));
			N2 = (Bx*pts1(2, i) - Bz * pts1(0, i)) / (pts1(0, i)*rpts2(2, 0) - rpts2(0, 0)*pts1(2, i));
			
			//�������ϵ��(phi,omega,kappa,u,v)
			A(i,0) = -rpts2(0, i)*rpts2(1, i)*N2 / rpts2(2, i);
			A(i,1) = -(rpts2(2, i) + rpts2(1, i)*rpts2(1, i) / rpts2(2, i))*N2;
			A(i,2) = rpts2(0, i)*N2;
			A(i,3) = Bx;
			A(i,4) = -rpts2(1, i)*Bx / rpts2(2, i);
			//���������Ӳ� Q,Ҳ�����̳�����L
			L(i, 0) = N1 * pts1(1, i) - N2 * rpts2(1, i) - By;
		}
		// ���
		if (useWieghtMat)
			X = (A.transpose()*(*in_weightMatrix) * A).inverse() * A.transpose()*(*in_weightMatrix)*L;
		else
			X = (A.transpose() * A).inverse() * A.transpose()*L;
		//����
		coefMat += X;
		if (printResult) {
			cout << "\n��" << t << "�ε���\n------------------------------------------" << endl;
			cout << "derta:  \t" << X.transpose().format(fullFormat) << endl
				<< "coefVec:\t" << coefMat.transpose().format(fullFormat) << endl;
		}
	} while (((abs(X(0)) > limit_angle) || (abs(X(1)) > limit_angle) || (abs(X(2)) > limit_angle)) && (t < maxIterTimes));

	coefVec = matrix2vectorD(coefMat);
	//**���
	if (t >= maxIterTimes) {
		cout << "\nWARNING\t�ﵽ������������(func: " << __func__ << ")" << endl
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