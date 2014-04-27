// SE3Tools.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <Eigen/Dense>
#include <iostream>
#include "SE3.h"
#include "Windows.h"

using Eigen::MatrixXd;
using Eigen::Vector3d;
using SO3::expm;
using SE3::expm;
using SO3::hat3;
using SE3::hat6;
using SE3::vee6;
using SE3::Adjoint;
using SE3::Vector6d;
using std::cout;
using std::endl;

LARGE_INTEGER t1,t2,freq;

int _tmain(int argc, _TCHAR* argv[])
{
	QueryPerformanceFrequency(&freq);

	Vector6d xi;
	xi << 1.0,0.0,0.5,0.0,0.0,1.0;
	QueryPerformanceCounter(&t1);


	Eigen::Matrix4d X = Eigen::Matrix4d::Identity();
	
	for (int i = 0; i < 10; ++i) {
		X = expm(hat6(xi))*X;
	}
	QueryPerformanceCounter(&t2);

	cout << X << endl;
	cout << endl << "Elapsed time: " << (double(t2.QuadPart)-double(t1.QuadPart))/(double(freq.QuadPart)) << endl;
	cout << SO3::log( (expm(hat6(xi))).block<3,3>(0,0) ) << endl;
	cout << endl;
	cout << SE3::log( (expm(hat6(xi))) ) << endl;
	return 0;
}
