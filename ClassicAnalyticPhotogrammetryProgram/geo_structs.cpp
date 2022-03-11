#include"geo_structs.h"

/*********************************************************************************************/
//							运算辅助函数
/*********************************************************************************************/
MatrixXd vector2matrixD(const VectorXd& vec)
{
	MatrixXd mat;
	mat.resize(vec.size(), 1);
	for (Eigen::Index i = 0; i < vec.size(); i++)
		mat(i) = vec(i);
	return mat;
}
VectorXd matrix2vectorD(const MatrixXd& mat)
{
	VectorXd vec;
	vec.resize(mat.size());
	for (Eigen::Index i = 0; i < mat.size(); i++)
		vec(i) = mat(i);
	return vec;
}
MatrixXd multi(const MatrixXd& mat1, const MatrixXd& mat2)
{
	assert(mat1.rows() == mat2.rows() && mat1.cols() == mat2.cols());
	MatrixXd m = mat1;
	for (Eigen::Index i = 0; i < m.rows(); i++){
		for (Eigen::Index j = 0; j < m.cols(); j++){
			m(i, j) *= mat2(i, j);
		}
	}
	return m;
}

/*******************************************************************************************/
//					共线方程模型【光束模型】（三点一线）
/**************************************************************************************/
void CBundleModel::caculateR() {						//旋转矩阵计算
	this->R = Rotation::rotate(phi, omega, kappa);
}

void CBundleModel::caculateXYZ_() {					//共线方程的分子、分母计算
	this->X_ = (R(0, 0) * (groundPt.x() - Xs) + R(1, 0) * (groundPt.y() - Ys) + R(2, 0) * (groundPt.z() - Zs));
	this->Y_ = (R(0, 1) * (groundPt.x() - Xs) + R(1, 1) * (groundPt.y() - Ys) + R(2, 1) * (groundPt.z() - Zs));
	this->Z_ = (R(0, 2) * (groundPt.x() - Xs) + R(1, 2) * (groundPt.y() - Ys) + R(2, 2) * (groundPt.z() - Zs));
}

void CBundleModel::caculteApproXY() {				//近似坐标计算
	approx = -f * X_ / Z_;
	approy = -f * Y_ / Z_;
}

void CBundleModel::caculate_A_and_L(MatrixXd &A, MatrixXd &L, int type)	//系数矩阵计算2x3大小
{
	//系数矩阵计算2x6 或 2*3大小
	//偏导顺序: 【线元素、角元素】
	//	2*6系数阵			--> type == 0
	//	Xs,Ys,Zs			--> type == 1
	//	phi,omega,kappa		--> type == 2
	//	X,Y,Z				--> type == 3
	//	[,f,,x0,y0]			--> type == 4

	////下面三行顺序不可改变
	//if (!this->RIsOK)this->caculateR();
	//if (!this->XYZ_IsOK)this->caculateXYZ_();
	//if (!this->approxyIsOK)this->caculteApproXY();

	assert(type >= 0 && type < 4);

	// [Xs,Ys,Zs],[X,Y,Z]
	if (type == 0 || type == 3) {
		double tmp = 1.0;
		if (type == 3) { 
			tmp = -1.0;
			A.resize(2, 3);
		}
		else A.resize(2, 6);
		A(0, 0) = tmp*((R(0, 0) * f + R(0, 2) * (photoPt.x() - x0)) / Z_);
		A(1, 0) = tmp*((R(0, 1) * f + R(0, 2) * (photoPt.y() - y0)) / Z_);
		A(0, 1) = tmp*((R(1, 0) * f + R(1, 2) * (photoPt.x() - x0)) / Z_);
		A(1, 1) = tmp*((R(1, 1) * f + R(1, 2) * (photoPt.y() - y0)) / Z_);
		A(0, 2) = tmp*((R(2, 0) * f + R(2, 2) * (photoPt.x() - x0)) / Z_);
		A(1, 2) = tmp*((R(2, 1) * f + R(2, 2) * (photoPt.y() - y0)) / Z_);
	}
	if (type != 3 && type != 1) {
		int tmp1 = 0;
		if (type == 0) tmp1 = 3; 
		else A.resize(2, 3);
		A(0, 0 + tmp1) = (photoPt.y() - y0) * sin(omega) - ((photoPt.x() - x0) * ((photoPt.x() - x0) * cos(kappa)
			- (photoPt.y() - y0) * sin(kappa)) / f + f * cos(kappa)) * cos(omega);
		A(1, 0 + tmp1) = -(photoPt.x() - x0) * sin(omega) - ((photoPt.y() - y0) * ((photoPt.x() - x0) * cos(kappa)
			- (photoPt.y() - y0) * sin(kappa)) / f - f * sin(kappa)) * cos(omega);
		A(0, 1 + tmp1) = -f * sin(kappa) - (photoPt.x() - x0) * ((photoPt.x() - x0) * sin(kappa) + (photoPt.y() - y0) * cos(kappa)) / f;
		A(1, 1 + tmp1) = -f * cos(kappa) - ((photoPt.y() - y0) * ((photoPt.x() - x0) * sin(kappa) + (photoPt.y() - y0) * cos(kappa))) / f;
		A(0, 2 + tmp1) = (photoPt.y() - y0);
		A(1, 2 + tmp1) = -(photoPt.x() - x0);
	}
	//else if(type==4){
	//	printf_s("EROOR:\tUnknow coefficient type!(func: %s)\n", __func__);
	//	exit(-1);
	//	//A(0, 6) = (photoPt.x - x0) / f;
	//	//A(1, 6) = (photoPt.y - y0) / f;
	//	//A(0, 7) = 1;
	//	//A(1, 7) = 0;
	//	//A(0, 8) = 0;
	//	//A(1, 8) = 1;
	//}
	// l11,l21
	L.resize(2, 1);
	L(0, 0) = photoPt.x() - x0 - approx;
	L(1, 0) = photoPt.y() - y0 - approy;
}


/*******************************************************************************************/
//					三维相似变换
/**************************************************************************************/
void similarityTransformationCaculate(
	const Matrix3Xd& relPts						//相对定向后的坐标【单位：m】
	, const Matrix3Xd& groudPts					//地面控制点坐标【单位：m】
	, const VectorXd& prmVec					//七参数向量（X0,Y0,Z0,Lamda,Phi,Omega,Kappa）【单位：m、rad】
	, MatrixXd& A								//相似变换方程系数矩阵
	, MatrixXd& L								//相似变换方程常数项矩阵
)
{	
	assert(relPts.cols() >= 5 && (relPts.cols() == groudPts.cols()) && (prmVec.size() == 7));

	//数据准备*****************************************
	//相对定向坐标旋转后
	Matrix3Xd rotedRelPts = Rotation::rotate(prmVec(4), prmVec(5), prmVec(6)) * relPts;
	//线偏移
	MatrixXd tmp= vector2matrixD(Vector3d(prmVec(0), prmVec(1), prmVec(2)));
	//近似地面点坐标
	MatrixXd approxGPts = prmVec(3) * rotedRelPts;
	for (int i = 0; i < approxGPts.cols(); i++){
		approxGPts.block(0, i, 3, 1) += tmp;
	}

	//求解*************************************
	A.resize(3 * relPts.cols(), 7);
	L.resize(3 * relPts.cols(), 1);
	A.fill(0.0);
	for (size_t i = 0; i < relPts.cols(); i++) {
		//偏导顺序
		//		X,	Y,	Z,	lamda,	phi,	omega,	kappa
		// X,Y,Z
		A(i * 3, 0) = A(i * 3 + 1, 1) = A(i * 3 + 2, 2) = 1.0;
		//lamda
		A.block(i * 3, 3, 3, 1) = rotedRelPts.block(0, i, 3, 1);
		//phi
		A(i * 3+0, 4) = -prmVec(3) * rotedRelPts(2,i);
		A(i * 3+2, 4) = prmVec(3) * rotedRelPts(0,i);
		//omega
		A(i * 3 + 0, 5) = -prmVec(3) * sin(prmVec(4))* rotedRelPts(1, i);
		A(i * 3 + 1, 5) = prmVec(3) * rotedRelPts(0, i)*sin(prmVec(4)) - prmVec(3) * rotedRelPts(2, i)*cos(prmVec(4));
		A(i * 3 + 2, 5) = prmVec(3) * rotedRelPts(1, i)*cos(prmVec(4));
		//kappa
		A(i * 3 + 0, 6) = -prmVec(3) *rotedRelPts(1, i)*cos(prmVec(4))* cos(prmVec(5)) - prmVec(3) * rotedRelPts(2, i)*sin(prmVec(5));
		A(i * 3 + 1, 6) = prmVec(3) * rotedRelPts(0, i)*cos(prmVec(4))*cos(prmVec(5))* +prmVec(3) * rotedRelPts(2, i)*sin(prmVec(4))*cos(prmVec(5));
		A(i * 3 + 2, 6) = prmVec(3) * rotedRelPts(0, i)*sin(prmVec(5)) - prmVec(3) * rotedRelPts(1, i)*cos(prmVec(5))* sin(prmVec(4));
		//求常数项矩阵
		L.block(i * 3, 0, 3, 1) = groudPts.block(0, i, 3, 1) - approxGPts.block(0, i, 3, 1);
	}
}



/*******************************************************************************/
//				转角矩阵
/*******************************************************************************/
inline 
Matrix3d Rotation::rotate(double phi, double omega, double kappa) {
	Matrix3d R;
	R(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
	R(0, 1) = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
	R(0, 2) = -sin(phi) * cos(omega);
	R(1, 0) = cos(omega) * sin(kappa);
	R(1, 1) = cos(omega) * cos(kappa);
	R(1, 2) = -sin(omega);
	R(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
	R(2, 1) = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
	R(2, 2) = cos(phi) * cos(omega);
	return R;
}
inline
Matrix3d Rotation::rotX(double angle) {
	Matrix3d tmp = Matrix3d::Zero();
	tmp(0, 0) = 1.0;
	tmp(1, 1) = cos(angle), tmp(1, 2) = -sin(angle);
	tmp(2, 1) = sin(angle), tmp(2, 2) = cos(angle);
	return tmp;
}
inline
Matrix3d Rotation::rotY(double angle) {
	Matrix3d tmp = Matrix3d::Zero();
	tmp(1, 1) = 1.0;
	tmp(0, 0) = cos(angle), tmp(0, 2) = -sin(angle);
	tmp(2, 0) = sin(angle), tmp(2, 2) = cos(angle);
	return tmp;
}
inline
Matrix3d Rotation::rotZ(double angle) {
	Matrix3d tmp = Matrix3d::Zero();
	tmp(2, 2) = 1.0;
	tmp(0, 0) = cos(angle), tmp(0, 1) = -sin(angle);
	tmp(1, 0) = sin(angle), tmp(1, 1) = cos(angle);
	return tmp;
}
inline
Matrix3d Rotation::rot_dX(double angle) {
	Matrix3d tmp = Matrix3d::Zero();
	tmp(1, 1) = -sin(angle), tmp(1, 2) = -sin(angle);
	tmp(2, 1) = cos(angle), tmp(2, 2) = -sin(angle);
	return tmp;
}
inline
Matrix3d Rotation::rot_dY(double angle) {
	Matrix3d tmp = Matrix3d::Zero();
	tmp(0, 0) = -sin(angle), tmp(0, 2) = -cos(angle);
	tmp(2, 0) = cos(angle), tmp(2, 2) = -sin(angle);
	return tmp;
}
inline
Matrix3d Rotation::rot_dZ(double angle) {
	Matrix3d tmp = Matrix3d::Zero();
	tmp(0, 0) = -sin(angle), tmp(0, 1) = -cos(angle);
	tmp(1, 0) = cos(angle), tmp(1, 1) = -sin(angle);
	return tmp;
}
inline
Matrix3d Rotation::rotMat(double angle1, double angle2, double angle3, int sysType) {
	switch (sysType) {
	case 0:
		return (Rotation::rotY(angle1)*Rotation::rotX(angle2)*Rotation::rotZ(angle3));
	case 1:
		return (Rotation::rotX(angle1)*Rotation::rotY(angle2)*Rotation::rotZ(angle3));
	case 2:
		return (Rotation::rotZ(angle1)*Rotation::rotY(angle2)*Rotation::rotZ(angle3));
	default:
		printf_s("ERROR\tfunc: %s -> Unknow angle system!\n",  __func__);
		return Matrix3d();
	}
}
