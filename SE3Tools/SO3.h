#include <Eigen/Dense>

#if ( defined(WIN32) || defined( _WIN64 ) ) && defined( BUILD_SHARED_LIB )
#define EXPORT_SYM	_declspec(dllexport)
#else
#define EXPORT_SYM
#endif

namespace SO3 {
	EXPORT_SYM	Eigen::Vector3d vee3(const Eigen::Matrix3d& What);
	EXPORT_SYM	Eigen::Matrix3d hat3(const Eigen::Vector3d& w);
	EXPORT_SYM	Eigen::Matrix3d log(const Eigen::Matrix3d& R);
	EXPORT_SYM	Eigen::Matrix3d expm(const Eigen::Matrix3d& What);
	EXPORT_SYM  double RotationAngle(const Eigen::Matrix3d& R);
	EXPORT_SYM  Eigen::Vector3d RotationAxis(const Eigen::Matrix3d& R);
	EXPORT_SYM  Eigen::Matrix3d AxisAngle(const Eigen::Vector3d& axis, double angle);
}