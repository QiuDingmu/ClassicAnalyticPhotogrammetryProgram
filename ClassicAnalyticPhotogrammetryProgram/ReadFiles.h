#pragma once
#include<vector>
#include<string>
#include<fstream>
#include"typeDefined.h"

using std::string;
using std::vector;

#define TT (0.017453292519943295)//度转弧度乘以该数

bool readResection(
	const string& filename
	, Vector3d& innerElement_in_milimeter
	, Matrix3Xd& groundPts_in_meter
	, Matrix2Xd& photoPts_in_milimeter
) {
	std::ifstream fp(filename);
	if (fp.fail()) {
		cout << "前方交会文件“" << filename << "”读取失败!\n" << endl;
		return false;
	}
	// 读取像主点，主距
	fp >> innerElement_in_milimeter(0) 
		>> innerElement_in_milimeter(1) 
		>> innerElement_in_milimeter(2);
	int num = 0;
	fp >> num;
	groundPts_in_meter.resize(3, num);
	photoPts_in_milimeter.resize(2, num);
	for (size_t i = 0; i < num; i++) {
		fp >> groundPts_in_meter(0, i) >> groundPts_in_meter(1, i) >> groundPts_in_meter(2, i)
			>> photoPts_in_milimeter(0, i) >> photoPts_in_milimeter(1, i);
	}
	fp.close();
	return true;
}


bool readInesection(
	const string& filename
	, string &ptName
	, vector<string>& photoNames
	, Vector3d& innerElement							//内方元素[x0,y0,f] 【单位：mm】
	, Matrix<double, 6, Eigen::Dynamic>& outerElements6X//外方元素【单位：m、rad】
	, Matrix2Xd& points2X								//对应像平面坐标【单位：mm】
) {
	std::ifstream fp(filename);
	if (fp.fail()) {
		cout << "前方交会文件“" << filename << "”读取失败!\n" << endl;
		return false;
	}
	// 读取像主点，主距
	fp >> innerElement(0) >> innerElement(1) >> innerElement(2);
	//读取点名,像片个数
	int num = 0, tmp = 0, s = 1;
	fp >> ptName >> num;
	photoNames.resize(num, "unknow");
	outerElements6X.resize(6, num);
	points2X.resize(2, num);
	string str;


	for (int i = 0; i < num; i++) {
		fp >> photoNames[i];
		for (int j = 3; j <= 5; j++) {
			fp >> str;
			if (str[0] == '-')
				s = -1, str.erase(str.begin());
			else
				s = 1;
			tmp = str.find('.');//小数点
			outerElements6X(j, i) = (double)std::stoi(str.substr(0, tmp));
			outerElements6X(j, i) += ((double)std::stoi(str.substr(tmp + 1, 2)) / 60.0);
			outerElements6X(j, i) += ((double)std::stoi(str.substr(tmp + 3, 2)) / 3600.0);
			outerElements6X(j, i) *= TT*s;
		}

		fp	>> outerElements6X(1, i) >> outerElements6X(0, i) >> outerElements6X(2, i)
			>> points2X(0, i) >> points2X(1, i);
	}

	fp.close();
	return true;
}

bool readAbOriention(
	const string& filename
	, vector<string>& ptNames
	, Matrix3Xd& phtPoints
	, Matrix3Xd& ctlPoints
) 
{
	std::ifstream fp(filename);
	if (fp.fail()) {
		cout << "绝对定向文件“" << filename << "”读取失败!\n" << endl;
		return false;
	}
	int num = 0;
	fp >> num;
	ptNames.resize(num);
	phtPoints.resize(3, num);
	ctlPoints.resize(3, num);
	for (int i = 0; i < num; i++){
		fp >> ptNames[i]
			>> ctlPoints(0, i) >> ctlPoints(1, i) >> ctlPoints(2, i)
			>> phtPoints(0, i) >> phtPoints(1, i) >> phtPoints(2, i)
			;
	}
	fp.close();
	return true;
}

bool readRelOriention(
	const string& filename
	, vector<string>& ptNames
	, string& phtName1
	, string& phtName2
	, Matrix2Xd& ptsL
	, Matrix2Xd& ptsR
)
{
	std::ifstream fp(filename);
	if (fp.fail()) {
		cout << "相对定向文件“" << filename << "”读取失败!\n" << endl;
		return false;
	}
	int num = 0;
	fp >> num;
	ptNames.resize(num);
	ptsL.resize(Eigen::NoChange, num), ptsR.resize(Eigen::NoChange, num);
	string tmp;
	//先读一个点
	fp >> ptNames[0] >> phtName1 >> ptsL(0, 0) >> ptsL(1, 0)
		>> tmp >> phtName2 >> ptsR(0, 0) >> ptsR(1, 0);
	for (int i = 1; i < num; i++) {
		fp >> ptNames[i] >> tmp >> ptsL(0, i) >> ptsL(1, i)
			>> tmp >> tmp >> ptsR(0, i) >> ptsR(1, i);
	}

	fp.close();
	return true;
}

bool readInnerOriention(
	const string& filename
	, Matrix2Xd& trgpts			//理论坐标
	, Matrix2Xd& pts			//实际坐标
	, double& pxl_size			//像素大小	
) 
{
	std::ifstream fp(filename);
	if (fp.fail()) {
		cout << "内定向文件“" << filename << "”读取失败!\n" << endl;
		return false;
	}
	//读取像素大小
	fp >> pxl_size;
	int num = 0;
	fp >> num;
	trgpts.resize(Eigen::NoChange, num);
	pts.resize(Eigen::NoChange, num);
	for (int i = 0; i < num; i++){
		fp >> trgpts(0, i) >> trgpts(1, i)
			>> pts(0, i) >> pts(1, i);
	}
	fp.close();
	return true;
}