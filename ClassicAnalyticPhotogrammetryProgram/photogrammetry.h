#pragma once
#ifndef PHOTOGRAMMETRY_H
#define PHOTOGRAMMETRY_H
#include"geo_structs.h"

#define LIMIT_ANGLE (3E-5)
#define MAX_ITER_TIMES	(20)


/*����ռ�󷽽��ᣨֻ֧�֦�-��-��ת��ϵͳ��
	���ص�λȨ�����(Ϊ��ֵʱ����ʧ��)		���ص�λȨ������λ��mm��
	@groundPts_in_meter						������Ƶ����꡾��λ��m��
	@photoPts_in_millimeter					���Ƶ��Ӧ��������꡾��λ��mm��
	@innerElement_in_millimeter				�ڷ�Ԫ��[x0,y0,f] ����λ��mm��
	@outerElements							����ⷽԪ��[��Ԫ�أ���Ԫ��] ����λ��m �� rad��
	@m										������Ʊ����߷�ĸ
	@printResult = false					����Ļ��ӡ���
	@*in_weightMatrix = NULL				����Ȩ��Ĭ�ϵ�λ��
	@*out_rotMatrix = NULL					�����ת����
	@anglesLimit = 3E-5						���ϵ������
	@maxIterTimes = 20U						�������������
*/											
double Resection(				
	const Matrix3Xd& groundPts_in_meter
	, const Matrix2Xd& photoPts_in_millimeter
	, const Vector3d innerElement_in_millimeter	
	, VectorXd& outerElements
	, const double m = -1.0
	, bool printResult = false
	, MatrixXd* in_weightMatrix = NULL
	, Matrix3d* out_rotMatrix = NULL
	, MatrixXd* out_A_Matrix = NULL
	, MatrixXd* out_V_Matrix = NULL
	, double anglesLimit = LIMIT_ANGLE
	, unsigned int maxIterTimes = MAX_ITER_TIMES
);

/*˫��ǰ�����᡾˫���ͶӰϵ������
	������ά���Ƶ����꡾��λ��m��
	@innerElement			�ڷ�Ԫ��[x0,y0,f] ����λ��mm��
	@outerElements1			Ƭ1�ⷽλԪ�ء���λ��m ��rad��
	@pt1					Ƭ1�ϵĶ�Ӧ������꡾��λ��mm��
	@outerElements2			Ƭ2�ⷽλԪ�ء���λ��m ��rad��
	@pt2					Ƭ2�ϵĶ�Ӧ������꡾��λ��mm��
*/
Vector3d pointPrjIntersection(
	const Vector3d innerElement
	, const VectorXd& outerElements1, const Vector2d& pt1
	, const VectorXd& outerElements2, const Vector2d& pt2
);

/*��Ƭǰ������-������ʽ
	�����������������λ:mm��
	@geo_pt						�����������꡾��λ��m��
	@innerElement				�ڷ�Ԫ��[x0,y0,f] ����λ��mm��
	@outerElements6X			�ⷽԪ��		����λ��m��rad��
	@points2X					��Ӧ��ƽ�����꡾��λ��mm��
	@printResult = false		����Ļ��ӡ���
	@in_weightMatrix = NULL		����Ȩ��Ĭ�ϵ�λ��
	@out_A_Matrix = NULL		���ϵ������
	@out_V_Matrix = NULL		�������������
	@limit_plane = 0.5			ƽ�澫�ȡ���λ��m��
	@limit_elevation = 1		�߳̾��ȡ���λ��m��
	@maxIterTimes = 20U			����������
*/
double Intersection(
	Vector3d &geo_pt
	, const Vector3d innerElement
	, const Matrix6Xd& outerElements6X
	, const Matrix2Xd& points2X
	, bool printResult = false
	, MatrixXd* in_weightMatrix = NULL
	, MatrixXd* out_A_Matrix = NULL
	, MatrixXd* out_V_Matrix = NULL
	, double limit_plane = 0.5
	, double limit_elevation = 1
	, unsigned int maxIterTimes = MAX_ITER_TIMES
);
//ǰ������-���Է�ʽ
double IntersectionLiner(				//���ص�λȨ�����
	Vector3d &geo_pt					//��������꡾��λ��m��
	, const Vector3d innerElement		//�ڷ�Ԫ��[x0,y0,f] ����λ��mm��
	, const Matrix6Xd& outerElements6X	//�ⷽԪ��		����λ��m��rad��
	, const Matrix2Xd& points2X			//��Ӧ��ƽ�����꡾��λ��mm��
	, bool printResult					//����Ļ��ӡ���
	, MatrixXd* out_A_Matrix = NULL		//���ϵ������
	, MatrixXd* out_V_Matrix = NULL		//�������������
);

/*���ڶ������������
	@prmVec					�ڶ������
	@pts_in_pixel			������������꡾��λ��pixel��
	@pxl_size				���ش�С����λ����ĩ��
	@type					�ڶ��򷽷���0:����;1:���Σ�2:˫���ԡ�
	@trgPts_inn_millimeter	����ʵ�����꡾��λ��mm��
*/
void internalTrans(
	const VectorXd& prmVec
	, const Matrix2Xd& pts_in_pixel
	, const double pxl_size
	, int type
	, Matrix2Xd& trgPts_inn_millimeter
);


 /*�ڶ���
	����x��y�����ϵĲв�rms����λ��΢�ס�
	@trgPts_in_millimeter	//�������꡾��λ��mm��
	@pts_in_pixel			//ʵ����������꡾��λ��pixel��
	@coefVec				//��¼�ڶ���Ĳ���
	@pxl_size				//���ش�С����λ��mm��
	@RMSs = NULL			//ÿ�����ϵĲв��λ��΢�ס�
	@*sigma0 = NULL			//��λȨ������λ��΢�ס�
	@type = 0				//�ڶ��򷽷���0:����;1:���Σ�2:˫���ԡ�
	@printResult = false	//��ӡ�м���
 */
 Vector2d internalOritation(			
	 const Matrix2Xd &trgPts_in_milimetr
	 , const Matrix2Xd &pts_in_pixel
	 , VectorXd& coefVec
	 , const double pxl_size
	 , Matrix2Xd* RMSs = NULL
	 , double* sigma0 = NULL
	 , int type = 0
	 , bool printResult = false
 );

/*����λ���Ʊ任������������
	@prmVec		//�߲���������X0,Y0,Z0,Lamda,Phi,Omega,Kappa������λ��m��rad��
	@relPts		//��Զ���������
	@geoPts		//������Ƶ�����
*/
void similarityTrans3D(
	const VectorXd& prmVec
	, const Matrix3Xd& relPts
	, Matrix3Xd& geoPts
);

/*���Ļ�����
	@points							���롢������꡾��λһ�¡�
	@center							�������
	@colCoordinate = true			true:ÿһ����һ������
	@centerSetted = NULL			�Դ�Ϊ���Ľ������Ļ�
*/
void centerPoints(
	MatrixXd& points
	, VectorXd& center
	, bool colCoordinate = true
	, const VectorXd* centerSetted = NULL
);



/*���Զ���						���ؾ��Զ��������
	@phtPoints					���Զ���������
	@ctrlPoints					��Ӧ���������
	@coefVec					���Զ���Ԫ��:�߲���������X0,Y0,Z0,Lamda,Phi,Omega,Kappa������λ��m��rad��
	@printResult = false		����Ļ��ӡ���
	@*in_weightMatrix = NULL	����Ȩ��Ĭ�ϵ�λ��
	@*out_A_Matrix = NULL		���ϵ������
	@*out_V_Matrix = NULL		�������������
	@limit_angle = 3E-5			�ǶȾ��ȡ���λ��rad��
	@int maxIterTimes = 20U		����������
*/
double absoluteOrientation(
	const Matrix3Xd& phtPoints
	, const Matrix3Xd& ctrlPoints
	, VectorXd& coefVec
	, bool printResult = false
	, MatrixXd* in_weightMatrix = NULL
	, MatrixXd* out_A_Matrix = NULL
	, MatrixXd* out_V_Matrix = NULL
	, double limit_angle = LIMIT_ANGLE
	, unsigned int maxIterTimes = MAX_ITER_TIMES
);


/*��������Զ���(˫��)			�����������������λ��΢�ס�
	@f							����(��������ͬ)	����λ��mm��
	@ptsL						������ƽ������	����λ��mm��
	@ptsR						������ƽ������	����λ��mm��
	@coefVec					��Զ���Ԫ��(phi,omega,kappa,mu,v)
	@printResult = false		����Ļ��ӡ���
	@in_weightMatrix = NULL		����Ȩ��Ĭ�ϵ�λ��
	@out_A_Matrix = NULL		���ϵ������
	@out_V_Matrix = NULL		�������������
	@limit_angle = 3E-5			�ǶȾ���			����λ��rad��
	@maxIterTimes = 20U			����������
*/
double relativeOrientation(
	const double f
	, const Matrix2Xd& ptsL
	, const Matrix2Xd& ptsR
	, VectorXd& coefVec
	, bool printResult = false
	, MatrixXd* in_weightMatrix = NULL
	, MatrixXd* out_A_Matrix = NULL
	, MatrixXd* out_V_Matrix = NULL
	, double limit_angle = LIMIT_ANGLE
	, unsigned int maxIterTimes = MAX_ITER_TIMES
);


#endif // !PHOTOGRAMMETRY_H
