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

	cout << "Testing angles near zero: " << endl;
	double err = 0;
	double errLog = 0;
	double dMax = 0;
	for (double d = 0.5; d >= 0; d -= 1e-6) {
		Eigen::Vector3d wSmall;
		wSmall(0) = d/sqrt(3.0);
		wSmall(1) = d/sqrt(3.0);
		wSmall(2) = d/sqrt(3.0);
		double e = fabs(wSmall.norm() - SO3::RotationAngle( SO3::expm( SO3::hat3( wSmall ) ) ) );
		if (e/wSmall.norm() > err) {
			err = e/wSmall.norm();
			dMax = d;
		}

		e = ( SO3::hat3(wSmall) - SO3::log( SO3::expm( SO3::hat3( wSmall ) ) ) ).norm();
		if (e > errLog)
			errLog = e/wSmall.norm();
	}

	cout << "Largest relative Angle Error: " << std::setprecision(10) << err << endl;
	cout << "d = " << dMax << endl;
	cout << "Largest relative Log Error: " << std::setprecision(10) << errLog << endl;

	cout << "Testing angles near pi: " << endl;
	err = 0;
	errLog = 0;
	Eigen::Matrix3d wHatActual;
	Eigen::Matrix3d wHatLog;
	double dMaxErr = 0;
	for (double d = 0.5; d >= 0; d -= 1e-6) {
		Eigen::Vector3d wPi;
		wPi(0) = -M_PI/sqrt(3.0) + d/sqrt(3.0);
		wPi(1) = -M_PI/sqrt(3.0) + d/sqrt(3.0);
		wPi(2) = M_PI/sqrt(3.0) - d/sqrt(3.0);
		double e = fabs(wPi.norm() - SO3::RotationAngle( SO3::expm( SO3::hat3( wPi ) ) ) );
		if (e/wPi.norm() > err)
			err = e/wPi.norm();

		e = ( SO3::hat3(wPi) - SO3::log( SO3::expm( SO3::hat3( wPi ) ) ) ).norm();
		if (e > errLog) {
			errLog = e/wPi.norm();

			wHatActual = SO3::hat3(wPi);
			wHatLog = SO3::log( SO3::expm( SO3::hat3( wPi ) ) );
			dMaxErr = d;
		}
	}

	cout << "Largest relative angle error: " << err << endl;
	cout << "Largest relative Log Error: " << std::setprecision(10) << errLog << endl;
	cout.precision(14);
	cout << "Max Log Error wHat: " << endl << wHatActual << endl;
	cout << "Max Log Error wHatLog: " << endl << wHatLog << endl;
	cout << "d = " << dMaxErr << endl;
	return 0;
}

