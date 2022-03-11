#pragma once
#ifndef GEO_STRUCTS_H
#define GEO_STRUCTS_H
#include"typeDefined.h"
#include<math.h>
#include<assert.h>



/***	���㸨������			********************************************************/
MatrixXd vector2matrixD(const VectorXd& vec);
VectorXd matrix2vectorD(const MatrixXd& mat);
MatrixXd multi(const MatrixXd& mat1, const MatrixXd& mat2);


/***	��ת�Ǿ���		*****************************************************************/
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
		@type:ת��ϵͳ
			0	PhiOmegaKappa		��YΪ�����ת��ϵͳ
			1	Omega_Phi_Kappa_	��XΪ�����ת��ϵͳ
			2	AAlphaKappa			��ZΪ�����ת��ϵͳ
	*/
	static inline Matrix3d rotMat(double angle1, double angle2, double angle3, int sysType = 0);

};



//***	���߷���ģ�͡�����ģ�͡�������һ�ߣ�	**************************************************
//	ע���������漰�������굥λ����ͬ���ס�
//		�Ƕȵ�λ��rad��
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

	void caculateR();						//������ת����
	void caculateXYZ_();					//���߷��̵ķ��ӡ���ĸ����
	void caculteApproXY();					//�����������
	void caculate_A_and_L(					//ϵ���������2x3��С
		MatrixXd &A, MatrixXd &L,int type
	);
public:
	double f, x0, y0;						//���࣬������
	double phi, omega, kappa, Xs, Ys, Zs;	//�ⷽԪ��
	double X_, Y_, Z_;						//���߷��̵ķ��ӡ���ĸ
	double approx, approy;					//�������õĽ�������
	Vector2d photoPt;						//���
	Vector3d groundPt;						//�����
	Matrix3d R;								//��ת����
}BUNDLE;



//��ά���Ʊ任ģ��
//	ע���������漰�������굥λ����ͬ���ס�
//		�Ƕȵ�λ��rad��
void similarityTransformationCaculate(
	const Matrix3Xd& relPts						//��Զ��������꡾��λ��m��
	, const Matrix3Xd& groudPts					//������Ƶ����꡾��λ��m��
	, const VectorXd& prmVec					//�߲���������X0,Y0,Z0,Lamda,Phi,Omega,Kappa������λ��m��rad��
	, MatrixXd& A								//���Ʊ任����ϵ������
	, MatrixXd& L								//���Ʊ任���̳��������
);


#endif // !GEO_STRUCTS_H
