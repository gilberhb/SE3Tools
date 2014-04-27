#include <Eigen/Dense>
#include "SE3.h"

using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector3d;

namespace SE3 {

EXPORT_SYM
Eigen::Matrix4d hat6(const Vector6d& xi)
{
	Matrix4d XiHat = Eigen::Matrix4d::Zero();
	XiHat.block(0,0,3,3) = SO3::hat3(xi.block(3,0,3,1));
	XiHat.block(0,3,3,1) = xi.block(0,0,3,1);
	return XiHat;
}


EXPORT_SYM
Vector6d vee6(const Eigen::Matrix4d& XiHat)
{
	Vector6d Xi;
	Xi.block<3,1>(0,0) = XiHat.block<3,1>(0,3);
	Xi.block<3,1>(3,0) = SO3::vee3( XiHat.block<3,3>(0,0) );
	return Xi;
}

static
double expmSE3_A(double theta)
{
	if (fabs(theta) > 1e-8) //return sin(theta)/theta
		return sin(theta)/theta; 
	else //return the truncated taylor series expansion
		return 1.0 - theta*theta/6.0 + theta*theta*theta*theta/120.0;
}

static 
double expmSE3_B(double theta)
{
	if (fabs(theta) > 1e-8) //return B factor
		return (1.0-cos(theta))/theta/theta;
	else //return the truncated taylor series for the B factor
		return 1.0/2.0 - theta*theta/24.0 + theta*theta*theta*theta/720.0;
}

static 
double expmSE3_C(double theta)
{
	if (fabs(theta) > 1e-8) //return the C factor
		return (1.0-expmSE3_A(theta))/theta/theta;
	else //return the truncated taylor series for the C factor
		return 1.0/6.0 - theta*theta/120.0 + theta*theta*theta*theta/5040.0;
}


Eigen::Matrix3d expmSE3_V(const Eigen::Matrix3d& What)
{
	double theta = 1.0/sqrt(2.0)*What.norm();
	return Matrix3d::Identity() + expmSE3_B(theta)*What + expmSE3_C(theta)*What;
}

EXPORT_SYM
Eigen::Matrix4d expm(const Eigen::Matrix4d& XiHat)
{
	Eigen::Matrix4d g = Eigen::Matrix4d::Identity();
	g.block<3,3>(0,0) = SO3::expm(XiHat.block(0,0,3,3));
	Eigen::Matrix3d V = expmSE3_V(XiHat.block(0,0,3,3));
	g.block<3,1>(0,3) = V*XiHat.block(0,3,3,1);
	return g;
}

EXPORT_SYM
Eigen::Matrix<double,6,6> Adjoint(const Eigen::Matrix4d& g)
{
	Eigen::Matrix<double,6,6>	Ad_g = Eigen::Matrix<double,6,6>::Identity();
	const Eigen::Block< const Eigen::Matrix<double,4,4>, 3, 3> R = g.block<3,3>(0,0);
	const Eigen::Block< const Eigen::Matrix<double,4,4>, 3, 1> p = g.block<3,1>(0,3);
	Ad_g.block<3,3>(0,0) = R;
	Ad_g.block<3,3>(3,3) = R;
	Ad_g.block<3,3>(0,3) = SO3::hat3( p ) * R;
	return Ad_g;
}

EXPORT_SYM
Eigen::Matrix4d log(const Eigen::Matrix4d& g)
{
	Eigen::Matrix4d Xi = Eigen::Matrix4d::Zero();
	Eigen::Matrix3d What = SO3::log( g.block<3,3>(0,0) );
	Xi.block<3,3>(0,0) = What;
	Xi.block<3,1>(0,3) = 
		expmSE3_V(What).partialPivLu().solve( g.block<3,1>(0,3) );
	return Xi;
}

}