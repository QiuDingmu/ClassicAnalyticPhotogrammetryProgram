#include"photogrammetry.h"
#include"ReadFiles.h"
#include"Timer.h"
using namespace std;


Eigen::IOFormat
fullFmt(		//全精度显示
	Eigen::FullPrecision, 0, ",\t", ";\n", "", "", "[", "]"
);


int func1(bool printResult = false);	//后方交会
int func2(bool printResult = false);	//前方交会
int func3(bool printResult = false);	//绝对定向
int func4(bool printResult = false);	//相对定向
int func5(bool printResult = false);	//  内定向

int main() {
	int n = 0, m = 0, k = 0;
	cout << "------\t\t选择显示的方法；是否显示中间过程\t\t------" << endl
		<< "选择显示的方法：" << endl
		<< "\t1  -> 后方交会" << endl
		<< "\t2  -> 前方交会" << endl
		<< "\t3  -> 绝对定向" << endl
		<< "\t4  -> 相对定向" << endl
		<< "\t5  -> 内定向" << endl
		<< "是否显示中间过程：" << endl
		<< "\t0 -> 否;\t 1 -> 是" << endl
		<< "输入:";
	cin >> n >> m;

	switch (n)	{
	case 1:k = func1(m ? true : false); break;//后方交会
	case 2:k = func2(m ? true : false); break;//前方交会
	case 3:k = func3(m ? true : false); break;//绝对定向
	case 4:k = func4(m ? true : false); break;//相对定向
	case 5:k = func5(m ? true : false); break;//内定向
	default:
		cout << "------\t\t未选择方法，退出!\t\t---------" << endl;
		break;
	}
	if (k < 0) {
		cout << "-----\t\t方法运行失败\t\t----------" << endl;
	}

	system("pause");
	return 0;
}

int func1(bool printResult) {
	Matrix3Xd G;
	Matrix2Xd P;
	VectorXd outerElements;
	Vector3d innerElement_in_milimeter(0.0, 0.0, 153.24);
	string fn = "后方交会.txt";
	if (!readResection(fn, innerElement_in_milimeter, G, P)) {
		return -1;
	}

	Timer tm;
	MatrixXd V;
	tm.Start();
	double m =
		Resection(
			G, P, innerElement_in_milimeter
			, outerElements, -1.0, printResult
			, NULL, NULL, NULL, &V
			//, 0.001
		);
	tm.Stop();
	if (m > 0) {
		cout << "外方元素:" << outerElements.transpose().format(fullFmt) << endl
		<< "单位权中误差：" << m << endl;
		cout << "V:" << V.transpose()*1000.0 << "微米" << endl;
	}
	else
		cout << "解算失败！" << endl;
	printf("\n\n用时: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}

int func2(bool printResult){	//前方交会
	string filename = "前方交会.txt", ptName = "unknow_pt";
	vector<string> photoNames;
	Vector3d innerElement;							//内方元素[x0,y0,f] 【单位：mm】
	Matrix6Xd outerElements6X;						//外方元素【单位：m、rad】
	Matrix2Xd points2X;								//对应像平面坐标【单位：mm】

	if (!readInesection(filename, ptName, photoNames, innerElement, outerElements6X, points2X)) {
		return -1;
	}
	
	double m = 0;
	cout << "\n前方交会计算方式：\n\t1:共线方程-迭代方式\n\t0:共线方程-线性方式" << endl;
	cin >> m;
	
	int cc = 0;
	cout << "\n是否进行重心化: 1(yes);0(no)" << endl;
	cin >> cc;
	
	bool centerG = false;
	if (cc)centerG = true;
	VectorXd pt0;
	if (centerG) {
		MatrixXd linerElem = outerElements6X.block(0, 0, 3, photoNames.size());
		centerPoints(linerElem, pt0);	//重心化
		outerElements6X.block(0, 0, 3, photoNames.size()) = linerElem;
	}
	if (printResult) {
		cout << "\n--------------\t前方交会 -> 输入数据显示\t-----------------" << endl;
		cout << "地面点点名:\t" << ptName << endl;
		cout << "内方元素：\t" << innerElement.transpose().format(fullFmt) << endl;
		for (size_t i = 0; i < photoNames.size(); i++) {
			cout << "像点 " << i + 1 << endl
				<< "\t所在片名：" << photoNames[i] << endl
				<< "\t外方元素：" << outerElements6X.block(0, i, 6, 1).transpose().format(fullFmt) << endl
				<< "\t坐标值：" << points2X.block(0, i, 2, 1).transpose().format(fullFmt) << endl;
		}
		if (centerG) {
			cout << "重心坐标： " << pt0.transpose().format(fullFmt) << endl;
		}
	}

	Vector3d geopt;
	Timer tm;
	tm.Start();
	if (m) {
		cout << "\n\n###   前方交会-迭代方式      #####################################" << endl;
		m = Intersection(
			geopt, innerElement, outerElements6X, points2X, printResult
		);
	}
	else {
		cout << "\n\n###   前方交会-线性方式      #####################################" << endl;
		m = IntersectionLiner(
			geopt, innerElement, outerElements6X, points2X, printResult
		);
	}
	tm.Stop();
	if (m > 0) {
		cout << "\n前方交会结果显示：--------------------------------" << endl
			<< "单位权中误差:\t" << m << " mm" << endl;
		if (centerG) {
			cout << "地面点重心化坐标:" << geopt.transpose().format(fullFmt) << endl
				<< "\t地面点坐标：" << (geopt + pt0).transpose().format(fullFmt) << endl;
		}
		else {
			cout << "地面点坐标:" << geopt.transpose().format(fullFmt) << endl;
		}
	}
	else
		cout << "\n前方交会解算失败！" << endl;
	printf("\n\n用时: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}

int func3(bool printResult) {
	string filename = "绝对定向.txt";
	vector<string> ptNames;
	Matrix3Xd phtPoints;
	Matrix3Xd ctlPoints;
	if (!readAbOriention(filename, ptNames, phtPoints, ctlPoints)) {
		return -1;
	}
	
	if (printResult) {
		cout << "\n绝对定向 -> 输入数据显示 ---------------------------------------------------" << endl;
		for (int i = 0; i < ptNames.size(); i++) {
			cout << ptNames[i] << ":" << endl
				<< phtPoints(0, i) << "\t" << phtPoints(1, i) << "\t" << phtPoints(2, i) << endl
				<< ctlPoints(0, i) << "\t" << ctlPoints(1, i) << "\t" << ctlPoints(2, i) << endl;
		}
	}
	VectorXd coefVec;
	Timer tm;
	tm.Start();
	double m0 = absoluteOrientation(
		phtPoints, ctlPoints, coefVec, printResult
	);
	tm.Stop();
	if (m0 > 0) {
		cout << "\n绝对定向结算结果显示： -----------------------------------------------" << endl;
		cout << "单位权中误差：" << m0 << "\t(单位与ctrlPoints相同)" << endl
			<< "七参数[X0,Y0,Z0,Lamda,Phi,Omega,Kappa]:\n\t"
			<< coefVec.transpose().format(fullFmt) << endl;
		Matrix3Xd re;
		similarityTrans3D(coefVec, phtPoints, re);
		cout << "残差：\n" << (ctlPoints - re).transpose().format(fullFmt) << endl;
	}
	else
		cout << "\n绝对定向解算错误!" << endl;
	printf("\n\n用时: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}

int func4(bool printResult) {
	const string filename = "相对定向.txt";
	vector<string> ptNames;
	string phtName1, phtName2;
	Matrix2Xd ptsL, ptsR;
	if (!readRelOriention(filename,ptNames,phtName1,phtName2,ptsL,ptsR)) {
		return -1;
	}
	if (printResult) {
		cout << "\n相对定向 -> 输入数据显示 ---------------------------------------------------" << endl;
		for (int i = 0; i < ptNames.size(); i++) {
			cout << "点号:" << ptNames[i] << "  --------------------" << endl
				<< "\t" << phtName1 << ":\t" << ptsL(0, i) << "\t" << ptsL(1, i) << endl
				<< "\t" << phtName2 << ":\t" << ptsR(0, i) << "\t" << ptsR(1, i) << endl << endl;
		}
	}

	VectorXd coefVec;
	MatrixXd  V;
	double f = 153.84;
	Timer tm;
	tm.Start();
	double  m0 = relativeOrientation(
		f, ptsL, ptsR, coefVec, printResult
		, NULL, NULL, &V
	);
	tm.Stop();
	cout << "\n相对定向解算结果显示 -------------------------------------" << endl
		<< "相对定向元素：[phi,omega,kappa,u,v]\n"
		<< coefVec.transpose().format(fullFmt) << endl
		<< "\n单位权中误差：" << m0 << " 微米 " << endl;
	V *= 1000.0;
	cout << "\n\n改正数/残差 ------------------" << endl;
	cout << "V: [微米]\n" << V.transpose().format(fullFmt) << endl;

	printf("\n\n用时: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}

int func5(bool printResult) {
	string filename = "内定向.txt";
	Matrix2Xd trgpts			//理论坐标
		, pts;					//实际坐标
	double pxl_size = 0;		//像素大小
	if (!readInnerOriention(filename,trgpts,pts,pxl_size)) {
		return -1;
	}

	if (printResult) {
		cout << "内定向像素大小:" << pxl_size << " mm" << endl
			<< "理论坐标: x\ty (mm)\t 量测坐标: x\ty (pxl)" << endl;
		for (int i = 0; i < pts.cols(); i++) {
			cout << trgpts(0, i) << "\t" << trgpts(1, i) << "\t"
				<< pts(0, i) << "\t" << pts(1, i) << endl;
		}
	}

	int n = 0;
	cout << "\n选择内定向的方法：[共有" << trgpts.cols() << "组点]" << endl
		<< "\t0:仿射;\t1:正形；\t2:双线性\n" << "选择:";
	cin >> n;

	VectorXd coefVec;
	Matrix2Xd RMSs;
	Timer tm;
	tm.Start();
	Vector2d
		m = internalOritation(
			trgpts, pts, coefVec, pxl_size, &RMSs, NULL
			,n ,printResult
		);
	tm.Stop();
	cout << "\nx、y方向上残差:\t[" << m.transpose().format(fullFmt) << "] 微米" << endl;
	printf("\n\n用时: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}