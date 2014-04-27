#include <Eigen/Dense>
#include "SO3.h"

#if ( defined(WIN32) || defined( _WIN64 ) ) && defined( BUILD_SHARED_LIB )
#define EXPORT_SYM	_declspec(dllexport)
#else
#define EXPORT_SYM
#endif

namespace SE3 {
	typedef Eigen::Matrix<double,6,1> Vector6d;
	
	EXPORT_SYM Eigen::Matrix4d hat6(const Vector6d& xi);
	EXPORT_SYM Vector6d vee6(const Eigen::Matrix4d& XiHat);
	EXPORT_SYM Eigen::Matrix4d expm(const Eigen::Matrix4d& XiHat);
	EXPORT_SYM Eigen::Matrix4d log(const Eigen::Matrix4d& g);
	EXPORT_SYM Eigen::Matrix<double,6,6> Adjoint(const Eigen::Matrix4d& g);

	//template <class Derived>
	//EXPORT_SYM const ::Eigen::Block<const Derived,3,3>	GetRotation( const ::Eigen::MatrixBase<Derived>& T )
	//{
	//	return T.derived().block<3,3>(0,0);
	//}
}