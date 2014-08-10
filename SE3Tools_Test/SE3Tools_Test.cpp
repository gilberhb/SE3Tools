// SE3Tools_Test.cpp : Defines the entry point for the console application.
//

#include <Eigen/Dense>
#include <iostream>
#include "../SE3Tools/SE3.h"
#include "Windows.h"
#undef min
#include <tchar.h>
#include <iomanip>
#include <limits>

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

const static double PI = 3.14159265358979323846264;

int _tmain(int argc, _TCHAR* argv[])
{
	try {

		Vector6d xi;
		xi << 1.0,0.0,0.5,0.0,1e-50,1e-50;


		Eigen::Matrix4d X = Eigen::Matrix4d::Identity();

		for (int i = 0; i < 10000; ++i) {
			X = expm(hat6(xi))*X;
		}

		cout << X << endl;
		cout << SE3::expm( hat6(xi) ) << endl;
		cout << SO3::log( (expm(hat6(xi))).block<3,3>(0,0) ) << endl;
		cout << SO3::RotationAxis( (expm(hat6(xi))).block<3,3>(0,0) ) << endl;
		cout << endl;
		cout << "hat6(xi) = " << endl << hat6(xi) << endl;
		cout << "expm(hat6(xi)) = " << endl << SE3::expm(SE3::hat6(xi)) << endl;
		cout << "Log(expm(hat6(xi))) = " << endl << SE3::log( SE3::expm(hat6(xi)) ) << endl;

		cout << "Error: " << endl << SE3::expm(SE3::log( SE3::expm(hat6(xi)) )) - SE3::expm(hat6(xi)) << endl;

		cout << "Testing angles near zero: " << endl;
		double err = 0;
		double errLog = 0;
		double dMax = 0;
		for (double d = 0.5; d >= 1e-127; d /= 1.001) {
			Eigen::Vector3d wSmall;
			wSmall(0) = d/sqrt(3.0);
			wSmall(1) = d/sqrt(3.0);
			wSmall(2) = d/sqrt(3.0);
			double e = fabs(wSmall.stableNorm() - SO3::RotationAngle( SO3::expm( SO3::hat3( wSmall ) ) ) );
			if (e/wSmall.stableNorm() > err) {
				err = e/wSmall.stableNorm();
				dMax = d;
			}

			e = ( SO3::hat3(wSmall) - SO3::log( SO3::expm( SO3::hat3( wSmall ) ) ) ).norm();
			if (e > errLog)
				errLog = e/wSmall.stableNorm();
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
		for (double d = 0.5; d >= 1e-50; d /= 1.001) {
			Eigen::Vector3d wPi;
			wPi(0) = -M_PI/sqrt(3.0) + d/sqrt(3.0);
			wPi(1) = -M_PI/sqrt(3.0) + d/sqrt(3.0);
			wPi(2) = M_PI/sqrt(3.0) - d/sqrt(3.0);
			double e = fabs(wPi.stableNorm() - SO3::RotationAngle( SO3::expm( SO3::hat3( wPi ) ) ) );
			if (e/wPi.stableNorm() > err)
				err = e/wPi.stableNorm();

			e = ( SO3::hat3(wPi) - SO3::log( SO3::expm( SO3::hat3( wPi ) ) ) ).norm();
			if (e > errLog) {
				errLog = e/wPi.stableNorm();

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

		{
			cout << "Testing Axis of Rotation for small rotations" << endl;
			Eigen::Vector3d w;
			w = Eigen::Vector3d::Zero();
			w(0) = 1e-100;
			w(1) = 1e-100;
			w(2) = 1e-100;
			cout << "Correct axis" << endl << w.normalized() << endl;
			cout << "Norm of correct axis: " << w.normalized().stableNorm() << endl;
			cout << "Axis of RotationAxis(expm(hat3(w))): " << endl << SO3::RotationAxis( SO3::expm( SO3::hat3( w ) ) ) << endl;
			cout << "Norm of RotationAxis(expm(hat3(w))): " << SO3::RotationAxis( SO3::expm( SO3::hat3( w ) ) ).stableNorm() << endl;
			try {
				cout << "Rotation Axis of Identity: " << endl << SO3::RotationAxis( Eigen::Matrix3d::Identity() ) << endl;
			} catch (std::runtime_error& e)	{
				cout << e.what() << endl;
			}

		}

		{
			Eigen::Vector3d w = Eigen::Vector3d::UnitX() + Eigen::Vector3d::UnitY();
			w = w.normalized().eval();
			cout.precision(20);
			cout << "Testing AxisAngle" << endl;
			cout << "AxisAngle(w,1e-320) " << endl << SO3::AxisAngle(w, 1e-320) << endl;
			cout << "expm(1e-320*hat3(w)) " << endl << SO3::expm( 1e-320*hat3(w) ) << endl;
			cout << "AxisAngle(w,0.1) - expm(0.1*hat3(w))" << endl << 
				(SO3::AxisAngle( w, 0.1 ) - SO3::expm( 0.1 * SO3::hat3(w) )) << endl;
			cout << "AxisAngle(w,0.1)" << endl <<
				SO3::AxisAngle(w,0.1) << endl;

			cout << "expm(0.1*hat3(w))" << endl << SO3::expm(0.1*hat3(w)) << endl;
		}

		{ 
			try {
				cout << "AxisAngle with the zero vector: " << endl;
				Eigen::Vector3d w = Eigen::Vector3d::Zero();
				cout << SO3::AxisAngle(w, 10) << endl;
			} catch (std::runtime_error& e) {
				cout << e.what() << endl;
			}
		}
	} catch (std::runtime_error& e) {
		cout << e.what() << endl;
	}

	return 0;
}

