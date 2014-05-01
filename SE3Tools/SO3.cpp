#include "SO3.h"
#include <iostream>

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
	if ( fabs(3 - R.trace()) > 1e-6 ) 
		return ProtectedAcos( (R.trace() - 1.0)/2.0 );
	else
		return 1.0/sqrt(2.0)*sqrt( R(0,2)*R(0,2) + R(0,1)*R(0,1) + R(1,2)*R(1,2) + R(1,0)*R(1,0) + R(2,0)*R(2,0) + R(2,1)*R(2,1) ); 
}

double logSO3_A(double theta)
{
	if ( fabs(theta) < 1e-8 )
		return 1.0/2.0 + theta*theta/12.0 + 7.0/720.0*theta*theta*theta*theta;
	else
		return theta/2.0/sin(theta);
}

//signum function, returns 1 for positive values
//and -1 for negative values, with 0 included in the positive half
double sgn(double a)
{
	return (a >= 0) ? 1 : -1;
}

Eigen::Matrix3d DecomposeWhat2(const Eigen::Matrix3d& R)
{
	return 1.0/expmSO3_B(RotationAngle(R)) * ( 1.0/2.0*(R + R.transpose()) - Eigen::Matrix3d::Identity() );
}

//Restricts answers to positive values,
//returns 0 on negative inputs
double ClampPos(double x)
{
	return (x >= 0) ? x : 0;
}

double logSpecial_AxisXMag(const Eigen::Matrix3d& what2)
{
	double gamma = -what2(2,2);
	double omega = -what2(0,0);
	double theta = -what2(1,1);
	return sqrt( ClampPos( (gamma - omega + theta)/2.0 ) );
}

double logSpecial_AxisYMag(const Eigen::Matrix3d& what2)
{
	double gamma = -what2(2,2);
	double omega = -what2(0,0);
	double theta = -what2(1,1);
	return sqrt( ClampPos( (gamma + omega - theta)/2.0 ) );
}

double logSpecial_AxisZMag(const Eigen::Matrix3d& what2)
{
	double gamma = -what2(2,2);
	double omega = -what2(0,0);
	double theta = -what2(1,1);
	return sqrt( ClampPos( (omega - gamma + theta)/2.0 ) );
}

Eigen::Vector3d logSpecial_AxisMagnitudes(const Eigen::Matrix3d& R)
{
	//We will get the values of |wx|, |wy|, |wz| from
	//what2
	Eigen::Matrix3d what2 = DecomposeWhat2(R);
	Eigen::Vector3d wmags;
	wmags(0) = logSpecial_AxisXMag(what2);
	wmags(1) = logSpecial_AxisYMag(what2);
	wmags(2) = logSpecial_AxisZMag(what2);

	return wmags;
}

//This handles the case when the rotation
//has an angle near pi
// TODO: NOT CURRENTLY IMPLEMENTED
Eigen::Matrix3d logSpecial(const Eigen::Matrix3d& R)
{
	//get the rotation vector magnitudes
	Eigen::Vector3d wmags = logSpecial_AxisMagnitudes(R);

	//don't want to divide by A(theta) becuase it is very small
	//used only to determine signs
	Eigen::Matrix3d Awhat = 1.0/2.0*( R - R.transpose() ); //equal to A(theta)*what
	Eigen::Vector3d Aw = SO3::vee3(Awhat);

	//now determine the signs individually from terms in Awhat
	Eigen::Vector3d w;
	w(0) = wmags(0) * sgn( Aw(0) );
	w(1) = wmags(1) * sgn( Aw(1) );
	w(2) = wmags(2) * sgn( Aw(2) );
	
	return SO3::hat3(w);
}

EXPORT_SYM	
Eigen::Matrix3d log(const Eigen::Matrix3d& R)
{
	//If the angle is not near PI, use the "simple" decomposition
	//that is easy and fast to compute. If we're near PI, then
	//we have to use a more complicated method, as logSO3_A(PI)
	//is zero, which gives the wrong answer. 
	double theta = RotationAngle(R);
	//use a tolerance that is generous, because 
	//logSpecial will work anywhere near theta = pi, but
	//has trouble around zero due to the computation of square roots
	//of tiny numbers
	if ( fabs(theta - M_PI) > 1e-3 ) 
		return logSO3_A(theta)*( R - R.transpose() );
	else
		return logSpecial(R); 
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