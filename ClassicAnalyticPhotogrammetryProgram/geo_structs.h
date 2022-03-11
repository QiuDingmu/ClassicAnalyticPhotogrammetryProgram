#pragma once
#ifndef GEO_STRUCTS_H
#define GEO_STRUCTS_H
#include"typeDefined.h"
#include<math.h>
#include<assert.h>



/***	运算辅助函数			********************************************************/
MatrixXd vector2matrixD(const VectorXd& vec);
VectorXd matrix2vectorD(const MatrixXd& mat);
MatrixXd multi(const MatrixXd& mat1, const MatrixXd& mat2);


/***	旋转角矩阵		*****************************************************************/
class Rotation {
public:
	Rotation(){}
	~Rotation(){}
	static inline Matrix3d rotate(double phi, double omega, double kappa);
	static inline Matrix3d rotX(double angle);
	static inline Matrix3d rotY(double angle);
	static inline Matrix3d rotZ(double angle);
	static inline Matrix3d rot_dX(double angle);
	static inline Matrix3d rot_dY(double angle);
	static inline Matrix3d rot_dZ(double angle);

	/*
		@type:转角系统
			0	PhiOmegaKappa		以Y为主轴的转角系统
			1	Omega_Phi_Kappa_	以X为主轴的转角系统
			2	AAlphaKappa			以Z为主轴的转角系统
	*/
	static inline Matrix3d rotMat(double angle1, double angle2, double angle3, int sysType = 0);

};



//***	共线方程模型【光束模型】（三点一线）	**************************************************
//	注：其中所涉及到的坐标单位均相同【米】
//		角度单位【rad】
typedef class CBundleModel {
public:
	CBundleModel():f(.0),x0(.0),y0(.0)
		, phi(.0), omega(.0), kappa(.0), Xs(.0), Ys(.0), Zs(.0)
		, X_(.0), Y_(.0), Z_(.0)
		,approx(.0),approy(.0)
	{}
	CBundleModel(const CBundleModel& model) {
		this->f = model.f, this->x0 = model.x0, this->y0 = model.y0;
		this->Xs = model.Xs, this->Ys = model.Ys, this->Zs = model.Zs;
		this->groundPt = model.groundPt, this->photoPt = model.photoPt;
	}
	void operator=(const CBundleModel& model) {
		this->f = model.f, this->x0 = model.x0, this->y0 = model.y0;
		this->Xs = model.Xs, this->Ys = model.Ys, this->Zs = model.Zs;
		this->groundPt = model.groundPt, this->photoPt = model.photoPt;
	}

	void caculateR();						//计算旋转矩阵
	void caculateXYZ_();					//共线方程的分子、分母计算
	void caculteApproXY();					//近似坐标计算
	void caculate_A_and_L(					//系数矩阵计算2x3大小
		MatrixXd &A, MatrixXd &L,int type
	);
public:
	double f, x0, y0;						//焦距，像主点
	double phi, omega, kappa, Xs, Ys, Zs;	//外方元素
	double X_, Y_, Z_;						//共线方程的分子、分母
	double approx, approy;					//计算所得的近似坐标
	Vector2d photoPt;						//像点
	Vector3d groundPt;						//地面点
	Matrix3d R;								//旋转矩阵
}BUNDLE;



//三维相似变换模型
//	注：其中所涉及到的坐标单位均相同【米】
//		角度单位【rad】
void similarityTransformationCaculate(
	const Matrix3Xd& relPts						//相对定向后的坐标【单位：m】
	, const Matrix3Xd& groudPts					//地面控制点坐标【单位：m】
	, const VectorXd& prmVec					//七参数向量（X0,Y0,Z0,Lamda,Phi,Omega,Kappa）【单位：m、rad】
	, MatrixXd& A								//相似变换方程系数矩阵
	, MatrixXd& L								//相似变换方程常数项矩阵
);


#endif // !GEO_STRUCTS_H
