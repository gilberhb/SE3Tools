#include "SO3.h"

using Eigen::Vector3d;
using Eigen::Matrix3d;

namespace SO3 {

EXPORT_SYM	
Eigen::Matrix3d hat3(const Eigen::Vector3d& w)
{
	Matrix3d What;
	What << 0.0, -w(2), w(1), 
		    w(2), 0.0, -w(0),
			-w(1), w(0), 0.0;
	return What;
}

EXPORT_SYM	
Eigen::Vector3d vee3(const Eigen::Matrix3d& What)
{
	Vector3d w (3);
	w(0) = (What(2,1) - What(1,2))/2.0;
	w(1) = (What(0,2) - What(2,0))/2.0;
	w(2) = (What(1,0) - What(0,1))/2.0;
	return w;
}

double ProtectedAcos(double arg)
{
	if (arg >= 1.0)
		return 0.0;
	else if (arg <= -1.0)
		return M_PI;
	else
		return acos( arg );
}

EXPORT_SYM	
double RotationAngle(const Eigen::Matrix3d& R)
{
	return ProtectedAcos( (R.trace() - 1.0)/2.0 );
}

double logSO3_A(double theta)
{
	if ( fabs(theta) < 1e-8 )
		return 1.0/2.0 + theta*theta/12.0 + 7.0/720.0*theta*theta*theta*theta;
	else
		return theta/2.0/sin(theta);
}

double sgn(double a)
{
	return (a > 0) ? 1 : -1;
}



//This handles the case when the rotation
//has an angle near pi
// TODO: NOT CURRENTLY IMPLEMENTED
Eigen::Matrix3d logSpecial(const Eigen::Matrix3d& R)
{
	double theta = RotationAngle(R);
	Eigen::Matrix3d	B = 1.0/2.0*( R + R.transpose() );
}

EXPORT_SYM	
Eigen::Matrix3d log(const Eigen::Matrix3d& R)
{
	double theta = RotationAngle(R);
	if ( fabs(theta - M_PI) > 1e-6 )
		return logSO3_A(theta)*( R - R.transpose() );
	else
		return logSpecial(R);
}

static
double expmSO3_A(double theta)
{
	//cutoff at 1e-3, because taylor series truncation will be 6th order 1e-3, ~1e-18
	if (fabs(theta) > 1e-3) //return sin(theta)/theta
		return sin(theta)/theta; 
	else //return the truncated taylor series expansion
		return 1.0 - theta*theta/6.0 + theta*theta*theta*theta/120.0;
}

static 
double expmSO3_B(double theta)
{
	if (fabs(theta) > 1e-3) //return B factor
		return (1.0-cos(theta))/theta/theta;
	else //return the truncated taylor series for the B factor
		return 1.0/2.0 - theta*theta/24.0 + theta*theta*theta*theta/720.0;
}

EXPORT_SYM	
Eigen::Matrix3d expm(const Eigen::Matrix3d& What)
{
	double theta = 1.0/sqrt(2.0)*What.norm();
	return Matrix3d::Identity() + expmSO3_A(theta)*What + expmSO3_B(theta)*What*What;
}

EXPORT_SYM  
Eigen::Matrix3d AxisAngle(const Eigen::Vector3d& axis, double angle)
{
	assert( axis.dot(axis) != 0 );
	return SO3::expm( angle/sqrt(axis.dot(axis)) * SO3::hat3(axis) );
}

EXPORT_SYM  
Eigen::Vector3d RotationAxis(const Eigen::Matrix3d& R)
{
	if (RotationAngle(R) == 0.0) {
		return Eigen::Vector3d::Zero();
	} else {
		return 1.0/RotationAngle(R) * SO3::vee3( SO3::log( R ) );
	}
}

}