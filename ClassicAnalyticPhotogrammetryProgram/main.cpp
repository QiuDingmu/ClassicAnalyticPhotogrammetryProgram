#include"photogrammetry.h"
#include"ReadFiles.h"
#include"Timer.h"
using namespace std;


Eigen::IOFormat
fullFmt(		//ȫ������ʾ
	Eigen::FullPrecision, 0, ",\t", ";\n", "", "", "[", "]"
);


int func1(bool printResult = false);	//�󷽽���
int func2(bool printResult = false);	//ǰ������
int func3(bool printResult = false);	//���Զ���
int func4(bool printResult = false);	//��Զ���
int func5(bool printResult = false);	//  �ڶ���

int main() {
	int n = 0, m = 0, k = 0;
	cout << "------\t\tѡ����ʾ�ķ������Ƿ���ʾ�м����\t\t------" << endl
		<< "ѡ����ʾ�ķ�����" << endl
		<< "\t1  -> �󷽽���" << endl
		<< "\t2  -> ǰ������" << endl
		<< "\t3  -> ���Զ���" << endl
		<< "\t4  -> ��Զ���" << endl
		<< "\t5  -> �ڶ���" << endl
		<< "�Ƿ���ʾ�м���̣�" << endl
		<< "\t0 -> ��;\t 1 -> ��" << endl
		<< "����:";
	cin >> n >> m;

	switch (n)	{
	case 1:k = func1(m ? true : false); break;//�󷽽���
	case 2:k = func2(m ? true : false); break;//ǰ������
	case 3:k = func3(m ? true : false); break;//���Զ���
	case 4:k = func4(m ? true : false); break;//��Զ���
	case 5:k = func5(m ? true : false); break;//�ڶ���
	default:
		cout << "------\t\tδѡ�񷽷����˳�!\t\t---------" << endl;
		break;
	}
	if (k < 0) {
		cout << "-----\t\t��������ʧ��\t\t----------" << endl;
	}

	system("pause");
	return 0;
}

int func1(bool printResult) {
	Matrix3Xd G;
	Matrix2Xd P;
	VectorXd outerElements;
	Vector3d innerElement_in_milimeter(0.0, 0.0, 153.24);
	string fn = "�󷽽���.txt";
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
		cout << "�ⷽԪ��:" << outerElements.transpose().format(fullFmt) << endl
		<< "��λȨ����" << m << endl;
		cout << "V:" << V.transpose()*1000.0 << "΢��" << endl;
	}
	else
		cout << "����ʧ�ܣ�" << endl;
	printf("\n\n��ʱ: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}

int func2(bool printResult){	//ǰ������
	string filename = "ǰ������.txt", ptName = "unknow_pt";
	vector<string> photoNames;
	Vector3d innerElement;							//�ڷ�Ԫ��[x0,y0,f] ����λ��mm��
	Matrix6Xd outerElements6X;						//�ⷽԪ�ء���λ��m��rad��
	Matrix2Xd points2X;								//��Ӧ��ƽ�����꡾��λ��mm��

	if (!readInesection(filename, ptName, photoNames, innerElement, outerElements6X, points2X)) {
		return -1;
	}
	
	double m = 0;
	cout << "\nǰ��������㷽ʽ��\n\t1:���߷���-������ʽ\n\t0:���߷���-���Է�ʽ" << endl;
	cin >> m;
	
	int cc = 0;
	cout << "\n�Ƿ�������Ļ�: 1(yes);0(no)" << endl;
	cin >> cc;
	
	bool centerG = false;
	if (cc)centerG = true;
	VectorXd pt0;
	if (centerG) {
		MatrixXd linerElem = outerElements6X.block(0, 0, 3, photoNames.size());
		centerPoints(linerElem, pt0);	//���Ļ�
		outerElements6X.block(0, 0, 3, photoNames.size()) = linerElem;
	}
	if (printResult) {
		cout << "\n--------------\tǰ������ -> ����������ʾ\t-----------------" << endl;
		cout << "��������:\t" << ptName << endl;
		cout << "�ڷ�Ԫ�أ�\t" << innerElement.transpose().format(fullFmt) << endl;
		for (size_t i = 0; i < photoNames.size(); i++) {
			cout << "��� " << i + 1 << endl
				<< "\t����Ƭ����" << photoNames[i] << endl
				<< "\t�ⷽԪ�أ�" << outerElements6X.block(0, i, 6, 1).transpose().format(fullFmt) << endl
				<< "\t����ֵ��" << points2X.block(0, i, 2, 1).transpose().format(fullFmt) << endl;
		}
		if (centerG) {
			cout << "�������꣺ " << pt0.transpose().format(fullFmt) << endl;
		}
	}

	Vector3d geopt;
	Timer tm;
	tm.Start();
	if (m) {
		cout << "\n\n###   ǰ������-������ʽ      #####################################" << endl;
		m = Intersection(
			geopt, innerElement, outerElements6X, points2X, printResult
		);
	}
	else {
		cout << "\n\n###   ǰ������-���Է�ʽ      #####################################" << endl;
		m = IntersectionLiner(
			geopt, innerElement, outerElements6X, points2X, printResult
		);
	}
	tm.Stop();
	if (m > 0) {
		cout << "\nǰ����������ʾ��--------------------------------" << endl
			<< "��λȨ�����:\t" << m << " mm" << endl;
		if (centerG) {
			cout << "��������Ļ�����:" << geopt.transpose().format(fullFmt) << endl
				<< "\t��������꣺" << (geopt + pt0).transpose().format(fullFmt) << endl;
		}
		else {
			cout << "���������:" << geopt.transpose().format(fullFmt) << endl;
		}
	}
	else
		cout << "\nǰ���������ʧ�ܣ�" << endl;
	printf("\n\n��ʱ: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}

int func3(bool printResult) {
	string filename = "���Զ���.txt";
	vector<string> ptNames;
	Matrix3Xd phtPoints;
	Matrix3Xd ctlPoints;
	if (!readAbOriention(filename, ptNames, phtPoints, ctlPoints)) {
		return -1;
	}
	
	if (printResult) {
		cout << "\n���Զ��� -> ����������ʾ ---------------------------------------------------" << endl;
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
		cout << "\n���Զ����������ʾ�� -----------------------------------------------" << endl;
		cout << "��λȨ����" << m0 << "\t(��λ��ctrlPoints��ͬ)" << endl
			<< "�߲���[X0,Y0,Z0,Lamda,Phi,Omega,Kappa]:\n\t"
			<< coefVec.transpose().format(fullFmt) << endl;
		Matrix3Xd re;
		similarityTrans3D(coefVec, phtPoints, re);
		cout << "�в\n" << (ctlPoints - re).transpose().format(fullFmt) << endl;
	}
	else
		cout << "\n���Զ���������!" << endl;
	printf("\n\n��ʱ: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}

int func4(bool printResult) {
	const string filename = "��Զ���.txt";
	vector<string> ptNames;
	string phtName1, phtName2;
	Matrix2Xd ptsL, ptsR;
	if (!readRelOriention(filename,ptNames,phtName1,phtName2,ptsL,ptsR)) {
		return -1;
	}
	if (printResult) {
		cout << "\n��Զ��� -> ����������ʾ ---------------------------------------------------" << endl;
		for (int i = 0; i < ptNames.size(); i++) {
			cout << "���:" << ptNames[i] << "  --------------------" << endl
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
	cout << "\n��Զ����������ʾ -------------------------------------" << endl
		<< "��Զ���Ԫ�أ�[phi,omega,kappa,u,v]\n"
		<< coefVec.transpose().format(fullFmt) << endl
		<< "\n��λȨ����" << m0 << " ΢�� " << endl;
	V *= 1000.0;
	cout << "\n\n������/�в� ------------------" << endl;
	cout << "V: [΢��]\n" << V.transpose().format(fullFmt) << endl;

	printf("\n\n��ʱ: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}

int func5(bool printResult) {
	string filename = "�ڶ���.txt";
	Matrix2Xd trgpts			//��������
		, pts;					//ʵ������
	double pxl_size = 0;		//���ش�С
	if (!readInnerOriention(filename,trgpts,pts,pxl_size)) {
		return -1;
	}

	if (printResult) {
		cout << "�ڶ������ش�С:" << pxl_size << " mm" << endl
			<< "��������: x\ty (mm)\t ��������: x\ty (pxl)" << endl;
		for (int i = 0; i < pts.cols(); i++) {
			cout << trgpts(0, i) << "\t" << trgpts(1, i) << "\t"
				<< pts(0, i) << "\t" << pts(1, i) << endl;
		}
	}

	int n = 0;
	cout << "\nѡ���ڶ���ķ�����[����" << trgpts.cols() << "���]" << endl
		<< "\t0:����;\t1:���Σ�\t2:˫����\n" << "ѡ��:";
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
	cout << "\nx��y�����ϲв�:\t[" << m.transpose().format(fullFmt) << "] ΢��" << endl;
	printf("\n\n��ʱ: <%4.2lf> ms  --------------------\n", tm.ElapsedTime());
	return 0;
}