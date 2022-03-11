#pragma once
#ifndef PHOTOGRAMMETRY_H
#define PHOTOGRAMMETRY_H
#include"geo_structs.h"

#define LIMIT_ANGLE (3E-5)
#define MAX_ITER_TIMES	(20)


/*单像空间后方交会（只支持φ-ω-κ转角系统）
	返回单位权中误差(为负值时程序失败)		返回单位权中误差【单位：mm】
	@groundPts_in_meter						地面控制点坐标【单位：m】
	@photoPts_in_millimeter					控制点对应的像点坐标【单位：mm】
	@innerElement_in_millimeter				内方元素[x0,y0,f] 【单位：mm】
	@outerElements							输出外方元素[线元素，角元素] 【单位：m 和 rad】
	@m										航摄设计比例尺分母
	@printResult = false					在屏幕打印结果
	@*in_weightMatrix = NULL				输入权阵，默认单位阵
	@*out_rotMatrix = NULL					输出旋转矩阵
	@anglesLimit = 3E-5						输出系数矩阵
	@maxIterTimes = 20U						输出改正数矩阵
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

/*双像前方交会【双像点投影系数法】
	返回三维控制点坐标【单位：m】
	@innerElement			内方元素[x0,y0,f] 【单位：mm】
	@outerElements1			片1外方位元素【单位：m 和rad】
	@pt1					片1上的对应像点坐标【单位：mm】
	@outerElements2			片2外方位元素【单位：m 和rad】
	@pt2					片2上的对应像点坐标【单位：mm】
*/
Vector3d pointPrjIntersection(
	const Vector3d innerElement
	, const VectorXd& outerElements1, const Vector2d& pt1
	, const VectorXd& outerElements2, const Vector2d& pt2
);

/*多片前方交会-迭代方式
	返回像点量测中误差【单位:mm】
	@geo_pt						所求地面点坐标【单位：m】
	@innerElement				内方元素[x0,y0,f] 【单位：mm】
	@outerElements6X			外方元素		【单位：m、rad】
	@points2X					对应像平面坐标【单位：mm】
	@printResult = false		在屏幕打印结果
	@in_weightMatrix = NULL		输入权阵，默认单位阵
	@out_A_Matrix = NULL		输出系数矩阵
	@out_V_Matrix = NULL		输出改正数矩阵
	@limit_plane = 0.5			平面精度【单位：m】
	@limit_elevation = 1		高程精度【单位：m】
	@maxIterTimes = 20U			最大迭代次数
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
//前方交会-线性方式
double IntersectionLiner(				//返回单位权中误差
	Vector3d &geo_pt					//地面点坐标【单位：m】
	, const Vector3d innerElement		//内方元素[x0,y0,f] 【单位：mm】
	, const Matrix6Xd& outerElements6X	//外方元素		【单位：m、rad】
	, const Matrix2Xd& points2X			//对应像平面坐标【单位：mm】
	, bool printResult					//在屏幕打印结果
	, MatrixXd* out_A_Matrix = NULL		//输出系数矩阵
	, MatrixXd* out_V_Matrix = NULL		//输出改正数矩阵
);

/*对内定向参数的运用
	@prmVec					内定向参数
	@pts_in_pixel			量测的像素坐标【单位：pixel】
	@pxl_size				像素大小【单位：毫末】
	@type					内定向方法【0:仿射;1:正形；2:双线性】
	@trgPts_inn_millimeter	返回实际坐标【单位：mm】
*/
void internalTrans(
	const VectorXd& prmVec
	, const Matrix2Xd& pts_in_pixel
	, const double pxl_size
	, int type
	, Matrix2Xd& trgPts_inn_millimeter
);


 /*内定向
	返回x、y方向上的残差rms【单位：微米】
	@trgPts_in_millimeter	//理论坐标【单位：mm】
	@pts_in_pixel			//实际量测的坐标【单位：pixel】
	@coefVec				//记录内定向的参数
	@pxl_size				//像素大小【单位：mm】
	@RMSs = NULL			//每个点上的残差【单位：微米】
	@*sigma0 = NULL			//单位权中误差【单位：微米】
	@type = 0				//内定向方法【0:仿射;1:正形；2:双线性】
	@printResult = false	//打印中间结果
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

/*对三位相似变换参数进行运用
	@prmVec		//七参数向量（X0,Y0,Z0,Lamda,Phi,Omega,Kappa）【单位：m、rad】
	@relPts		//相对定向后的坐标
	@geoPts		//地面控制点坐标
*/
void similarityTrans3D(
	const VectorXd& prmVec
	, const Matrix3Xd& relPts
	, Matrix3Xd& geoPts
);

/*重心化坐标
	@points							输入、输出坐标【单位一致】
	@center							输出重心
	@colCoordinate = true			true:每一列是一个坐标
	@centerSetted = NULL			以此为重心进行重心化
*/
void centerPoints(
	MatrixXd& points
	, VectorXd& center
	, bool colCoordinate = true
	, const VectorXd* centerSetted = NULL
);



/*绝对定向						返回绝对定向中误差
	@phtPoints					绝对定向后的坐标
	@ctrlPoints					对应地面点坐标
	@coefVec					绝对定向元素:七参数向量（X0,Y0,Z0,Lamda,Phi,Omega,Kappa）【单位：m、rad】
	@printResult = false		在屏幕打印结果
	@*in_weightMatrix = NULL	输入权阵，默认单位阵
	@*out_A_Matrix = NULL		输出系数矩阵
	@*out_V_Matrix = NULL		输出改正数矩阵
	@limit_angle = 3E-5			角度精度【单位：rad】
	@int maxIterTimes = 20U		最大迭代次数
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


/*连续法相对定向(双像)			返回像点量测中误差【单位：微米】
	@f							主距(左右像相同)	【单位：mm】
	@ptsL						左像相平面坐标	【单位：mm】
	@ptsR						右像相平面坐标	【单位：mm】
	@coefVec					相对定向元素(phi,omega,kappa,mu,v)
	@printResult = false		在屏幕打印结果
	@in_weightMatrix = NULL		输入权阵，默认单位阵
	@out_A_Matrix = NULL		输出系数矩阵
	@out_V_Matrix = NULL		输出改正数矩阵
	@limit_angle = 3E-5			角度精度			【单位：rad】
	@maxIterTimes = 20U			最大迭代次数
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
