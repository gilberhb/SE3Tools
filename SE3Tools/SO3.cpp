#include "SO3.h"
#include <iostream>

using Eigen::Vector3d;
using Eigen::Matrix3d;

namespace SO3 {

//normalize a 3-vector
Eigen::Vector3d Normalize3d(const Eigen::Vector3d& in)
{
	//stableNorm uses a hypot function to compute length
	// to avoid potential underflows or overflows
	return in/in.stableNorm(); 
}

//convert 3-vector to skew-symmetric matrix
EXPORT_SYM	
Eigen::Matrix3d hat3(const Eigen::Vector3d& w)
{
	Matrix3d What;
	What << 0.0, -w(2), w(1), 
		    w(2), 0.0, -w(0),
			-w(1), w(0), 0.0;
	return What;
}

//convert skew-symmetric matrix to 3-vector
EXPORT_SYM	
Eigen::Vector3d vee3(const Eigen::Matrix3d& What)
{
	Vector3d w (3);
	w(0) = (What(2,1) - What(1,2))/2.0;
	w(1) = (What(0,2) - What(2,0))/2.0;
	w(2) = (What(1,0) - What(0,1))/2.0;
	return w;
}

//A coefficient in exponential formula
static
double expmSO3_A(double theta)
{
	//cutoff at 1e-3, because taylor series truncation will be 6th order 1e-3, ~1e-18
	if (fabs(theta) > 1e-3) //return sin(theta)/theta
		return sin(theta)/theta; 
	else //return the truncated taylor series expansion
		return 1.0 - theta*theta/6.0 + theta*theta*theta*theta/120.0;
}

//B coefficient in exponential formula
static 
double expmSO3_B(double theta)
{
	if (fabs(theta) > 2e-3) //return B factor
		return (1.0-cos(theta))/theta/theta;
	else //return the truncated taylor series for the B factor
		return 1.0/2.0 - theta*theta/24.0 + theta*theta*theta*theta/720.0;
}

//Compute the matrix exponential of the matrix
//What, which is assumed to be a skew-symmetric matrix
EXPORT_SYM	
Eigen::Matrix3d expm(const Eigen::Matrix3d& What)
{
	double theta = 1.0/sqrt(2.0)*What.norm();
	return Matrix3d::Identity() + expmSO3_A(theta)*What + expmSO3_B(theta)*What*What;
}

//Protected ArcCos function
//so that it can't blow up in our face
double ProtectedAcos(double arg)
{
	if (arg >= 1.0)
		return 0.0;
	else if (arg <= -1.0)
		return M_PI;
	else
		return acos( arg );
}

//coefficient A in log expression
double logSO3_A(double theta)
{
	if ( fabs(theta) < 1e-3 ) //this is good down to about 1e-18 even at theta = 1e-3
		return 1.0/2.0 + theta*theta/12.0 + 7.0/720.0*theta*theta*theta*theta;
	else
		return theta/sin(theta)/2.0;
}

//Returns the cosine of the angle of rotation
//for a rotation matrix
double Rotation_CosineAngle(const Eigen::Matrix3d& R)
{
	return (R.trace() - 1.0)/2.0;
}

//Returns the sine of the angle of rotation
//for a rotation matrix
double Rotation_SineAngle(const Eigen::Matrix3d& R)
{
	return ( (R-R.transpose()).norm()/2.0/sqrt(2.0) );
}

//Compute the angle of rotation for 
//a rotation matrix R.  Use atan2 to avoid
//problems with either acos or asin.  
EXPORT_SYM	
double RotationAngle(const Eigen::Matrix3d& R)
{
	return atan2(Rotation_SineAngle(R), Rotation_CosineAngle(R));		//it is impossible for both arguments to be zero, so no chance of a blow-up
}

//This handles the case when the rotation
//has an angle near pi
// This method works by computing the logarithm of R^2.  Since
// it will only be called when the angle is within ~0.5 radian of
// pi, we know that R^2 will be a rotation which has the opposite
// axis and an angle of (pi-RotationAngle(R))*2.  This information
// is used to compute the log of R. This turns out to be more
// accurate than computing the log of R by more complicated means
// (at least in 64-bit IEEE floating point math).
Eigen::Matrix3d logSpecial(const Eigen::Matrix3d& R)
{
	Eigen::Matrix3d wwHat = SO3::log(R*R); //compute the log of R^2, which can be done very accurately
	double theta = wwHat.norm()/sqrt(2.0);
	return -wwHat/theta*(M_PI-theta/2.0);
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
	//logSpecial will work anywhere near theta = pi, whereas
	//the inaccuracy in logSO3_A becomes quite bad near theta = pi
	if ( fabs(theta - M_PI) > 0.5 ) //if the angle is within 0.5 radian of pi, use logSpecial 
		return logSO3_A(theta)*( R - R.transpose() );
	else
		return logSpecial(R); 
}

//Form the rotation matrix which represents a rotation
//about the axis given by the angle given. The axis
//may not be zero.
EXPORT_SYM  
Eigen::Matrix3d AxisAngle(const Eigen::Vector3d& axis, double angle)
{
	assert( axis.dot(axis) != 0 );
	return SO3::expm( angle * SO3::hat3( Normalize3d(axis) ) );
}



//Compute the axis of rotation for the rotation matrix R
//Note that this becomes very ill-conditioned when the rotation
//axis is small, whereas the logarithm doesn't suffer this
//problem. 
EXPORT_SYM  
Eigen::Vector3d RotationAxis(const Eigen::Matrix3d& R)
{
	if (RotationAngle(R) == 0.0) {
		return Eigen::Vector3d::Zero();
	} else {
		return Normalize3d(SO3::vee3( SO3::log( R ) )); //log is accurate, so use it
	}
}

EXPORT_SYM  
Eigen::Matrix3d RotationX(double angle)
{
	Eigen::Matrix3d R;
	R << 1,      0,       0,
		 0, cos(angle), -sin(angle),
		 0, sin(angle), cos(angle);
	return R;
}

EXPORT_SYM  
Eigen::Matrix3d RotationY(double angle)
{
	Eigen::Matrix3d R;
	R << cos(angle), 0, sin(angle),
		 0,          1, 0,
		 -sin(angle),0, cos(angle);
	return R;
}

EXPORT_SYM  
Eigen::Matrix3d RotationZ(double angle)
{
	Eigen::Matrix3d R;
	R << cos(angle), -sin(angle), 0,
		 sin(angle), cos(angle), 0,
		 0,          0,          1;
	return R;
}

}