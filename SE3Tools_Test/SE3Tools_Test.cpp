// SE3Tools_Test.cpp : Defines the entry point for the console application.
//

#include <Eigen/Dense>
#include <iostream>
#include "../SE3Tools/SE3.h"
#include "Windows.h"
#include <tchar.h>
#include <iomanip>

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
const static double PI = 3.14159265358979323846264;

int _tmain(int argc, _TCHAR* argv[])
{
	QueryPerformanceFrequency(&freq);

	Vector6d xi;
	xi << 1.0,0.0,0.5,0.0,PI-1e-10,0.0;
	QueryPerformanceCounter(&t1);


	Eigen::Matrix4d X = Eigen::Matrix4d::Identity();
	
	for (int i = 0; i < 10000; ++i) {
		X = expm(hat6(xi))*X;
	}
	QueryPerformanceCounter(&t2);

	cout << X << endl;
	cout << endl << "Elapsed time: " << (double(t2.QuadPart)-double(t1.QuadPart))/(double(freq.QuadPart)) << endl;
	cout << SE3::expm( hat6(xi) ) << endl;
	cout << SO3::log( (expm(hat6(xi))).block<3,3>(0,0) ) << endl;
	cout << SO3::RotationAxis( (expm(hat6(xi))).block<3,3>(0,0) ) << endl;
	cout << endl;
	cout << "hat6(xi) = " << endl << hat6(xi) << endl;
	cout << "Log(expm(hat6(xi))) = " << endl << SE3::log( SE3::expm(hat6(xi)) ) << endl;

	cout << "Error: " << endl << SE3::expm(SE3::log( SE3::expm(hat6(xi)) )) - SE3::expm(hat6(xi)) << endl;

	cout << "Testing small angles: " << endl;
	Eigen::Vector3d wSmall;
	wSmall(0) = 1e-1;
	wSmall(1) = 1e-1;
	wSmall(2) = 1e-1;

	cout << "Actual Angle norm(w): " << std::setprecision(10) << wSmall.norm() << endl;
	cout << "RotationAngle(expm(hat3(w))): " << std::setprecision(10) << SO3::RotationAngle( SO3::expm( SO3::hat3( wSmall ) ) ) << endl;
	return 0;
}

