#include <Eigen/Dense>

#if ( defined(WIN32) || defined( _WIN64 ) ) && defined( BUILD_SHARED_LIB )
#define EXPORT_SYM	_declspec(dllexport)
#else
#define EXPORT_SYM
#endif

/*!
 * In this documentation, group elements of SO(3) will be called R. The functions work
 * with the representation of SO(3) in the space of 3 x 3  real rotation matrices. These rotation
 * matrices satisfy the following two properties
 * \f[
 *    R R^T = I
 * \f]
 * and
 * \f[
 *    \det R = +1
 * \f]
 * The elements of so(3), the Lie algebra of SO(3), or the tangent space of SO(3) at the identity element, are represented
 * by skew symmetric real matrices. These "cross product matrices" are detailed in SO3::hat3.
 */
namespace SO3 {
	/*!
	 * \brief Convert a 3 x 3 skew symmetric real matrix \f$ \hat{\omega} \f$ to its representation as a 3-vector.
	 * 
	 * This is the inverse operation from hat3.
	 *
	 * \param What		The 3 x 3 skew symmetric matrix.
	 * \return w    	The 3-vector
	 */
	EXPORT_SYM	Eigen::Vector3d vee3(const Eigen::Matrix3d& What);

	/*!
	 * \brief Convert a 3d vector \f$ \omega \f$ to its skew-symmetric cross product matrix
	 * 
	 * This function converts a vector \f$ \omega = \begin{bmatrix} \omega_x & \omega_y & \omega_z \end{bmatrix}^T \f$
	 * into the matrix representing the linear transformation \f$ A(x) = \omega \times x \f$. This matrix is given by
	 * \f[
	 *    \hat{\omega} = \begin{bmatrix} 0 & -\omega_z & \omega_y \\
	                                    \omega_z & 0 & -\omega_x \\
										-\omega_y & \omega_x & 0 \end{bmatrix}
	 * \f]
	 *
	 * \param w		The 3d vector which will be converted to a cross product matrix.
	 * \return w^	The 3 x 3 cross product matrix. 
	 */
	EXPORT_SYM	Eigen::Matrix3d hat3(const Eigen::Vector3d& w);

	/*!
	 * \brief	The logarithm map on SO(3).
	 * 
	 * This function implements a single branch of the logarithm on SO(3), which returns
	 * a 3 x 3 skew symmetric matrix \f$ \hat{\omega} = \log R \f$. The logarithm function is 
	 * multi-valued, but this function returns only the branch for which the rotation angle
	 * \f$ \theta = || \omega || \f$ is less than or equal to \f$ \pi \f$. When \f$ R = I \f$
	 * is input the result is the zero matrix. The result of this function is uniquely defined,
	 * and the same answer is given for any particular input on every function call. Note that
	 * the function is not continuous in the sense that small changes in R can result in large
	 * changes in the output WHat. It is safe (and the result is accurate) to call this function
	 * with rotations that have an angle near or exactly 0. With angles near 
	 * \f$ \pi \f$ the function may execute more slowly as more work is required to accurately 
	 * compute the logarithm in these cases. Note that near \f$ \theta= \pi \f$ the SO3::RotationAngle
	 * function becomes less accurate, and since this function calls that function to compute
	 * the angle the results near \f$ \theta = \pi \f$ may not be as accurate as elsewhere in
	 * the space of rotations.
	 *
	 * \param R		The rotation matrix for which the logarithm is desired
	 * \return WHat	A skew-symmetric matrix 
	 */
	EXPORT_SYM	Eigen::Matrix3d log(const Eigen::Matrix3d& R);

	/*!
	 * \brief	The exponential map from se(3) to SO(3).
	 * 
	 * This function implements the exponential map which takes 3 x 3
	 * skew symmetric matrices and maps them onto the rotation matrices,
	 * using \f$ R = \exp \hat{\omega} \f$. The vector \f$ \omega \f$ may be
	 * decomposed into the axis of rotation and the rotation angle with
	 * \f$ \theta = ||\omega|| \f$ and \f$ a = \omega / \theta \f$, where 
	 * \f$ a \f$ is the unit vector pointing in the direction of rotation in
	 * the sense of the right-hand rule. The function is implemented with
	 * Rodriguez' formula
	 * \f[
	 *		R = I + \frac{\sin(\theta)}{\theta} \hat{\omega} + \frac{(1 - \cos(\theta))}{\theta^2} \hat{\omega}^2
	 * \f]
	 * 
	 * Note that \f$ \hat{\omega} \f$ is used to encode both the rotation axis and
	 * the angle, and \f$ \theta \f$ is computed from the input What. Near angles
	 * of zero the truncated taylor expansions of the coefficients are used instead of their
	 * trigonometric representations to retain accuracy (the limits as \f$ \theta \rightarrow 0 \f$ are finite).
	 * Angles larger than \f$ \pi \f$ are easily handled by Rodriguez' formula, but note that 
	 * precision loss will occur as the uncertainty in the angle grows with increasingly large
	 * angles, simply due to floating accuracy.
	 *
	 * \param  What	The skew symmetric matrix.
	 * \return R	The rotation matrix.
	 */
	EXPORT_SYM	Eigen::Matrix3d expm(const Eigen::Matrix3d& What);

	/*!
	 * \brief  Compute the angle of rotation (modulo \f$ 2\pi \f$) for a given rotation matrix.
	 * 
	 * This function computes the angle of rotation using the formula
	 * \f[
	 *		\theta = \cos^{-1} \left( \frac{\text{trace}(R)-1}{2} \right)
	 * \f]
	 *
	 * For rotation matrices which represent angles near zero, the above function
	 * is not used. For these cases, the argument to \f$ \cos^{-1} \f$ is nearly 1
	 * and instead the off-diagonal elements of \f$ R \f$ are used
	 * to compute the rotation angle, since \f$ R \approx I + \hat{\omega} \f$.
	 * The cutoff is currently hardcoded for when \f$ |\text{trace}(R) - 3| <= 1 \times 10^{-6} \f$.
	 * For angles near \f$ \pi \f$ the result may have some error as no special handling of this case
	 * has been implemented. 
	 *
	 * \param  R	The rotation matrix.
	 * \return theta	The angle of rotation.
	 */
	EXPORT_SYM  double RotationAngle(const Eigen::Matrix3d& R);

	/*!
	 * \brief  Extract the Axis of Rotation from a rotation matrix \f$ R \f$.
	 *
	 * This function returns the axis of rotation for a rotation matrix. Note that
	 * this function has a singularity at the Identity rotation, where the axis
	 * is undefined. In this case, the function returns the zero vector. Otherwise,
	 * the normalized axis of rotation is computed as \f$n = \frac{ \text{vee3}( \log(R) )}{|| \text{vee3}( \log(R) ) ||} \f$
	 *
	 * \param  R	The rotation matrix.
	 * \return n	The normalized axis of rotation
	 */
	EXPORT_SYM  Eigen::Vector3d RotationAxis(const Eigen::Matrix3d& R);

	/*!
	 * \brief  Compute a rotation matrix from the axis angle representation.
	 *
	 * This function returns the rotation matrix represented by a rotation of a given angle
	 * about a given axis. The axis is normalized prior to computation of the rotation matrix,
	 * so changes in scaling of the parameter axis don't change the result returned by this function.
	 *
	 * \param  axis	 The rotation axis (Will be normalized prior to use)
	 * \param  angle The rotation angle.
	 * \return R	 The rotation matrix.
	 */
	EXPORT_SYM  Eigen::Matrix3d AxisAngle(const Eigen::Vector3d& axis, double angle);

	/*!
	 * \brief Compute the canonical rotation about the \f$ e_1 \f$-axis, or \f$x\f$-axis.
	 *
	 * This function computes the rotation about the first canonical basis vector, with
	 * the angle given by the user. This matrix is given by
	 *
	 * \f[
	 *		R = \begin{bmatrix} 1 & 0 & 0 \\ 0 & \cos(\theta) & -\sin(\theta) \\ 0 & \sin(\theta) & \cos(\theta) \end{bmatrix}
	 * \f]
	 *
	 * \param  angle	 The rotation angle.
	 * \return R	     The rotation matrix.
	 */
	EXPORT_SYM  Eigen::Matrix3d RotationX(double angle);

	/*!
	 * \brief Compute the canonical rotation about the \f$ e_2 \f$-axis, or \f$y\f$-axis.
	 *
	 * This function computes the rotation about the second canonical basis vector, with
	 * the angle given by the user. This matrix is given by
	 *
	 * \f[
	 *		R = \begin{bmatrix} \cos(\theta) & 0 & \sin(\theta) \\ 0 & 1 & 0 \\ -\sin(\theta) & 0 & cos(\theta) \end{bmatrix}
	 * \f]
	 *
	 * \param  angle	 The rotation angle.
	 * \return R	     The rotation matrix.
	 */
	EXPORT_SYM  Eigen::Matrix3d RotationY(double angle);

	/*!
	 * \brief Compute the canonical rotation about the \f$ e_3 \f$-axis, or \f$z\f$-axis.
	 *
	 * This function computes the rotation about the third canonical basis vector, with
	 * the angle given by the user. This matrix is given by
	 *
	 * \f[
	 *		R = \begin{bmatrix} \cos(\theta) & -\sin(\theta) & 0 \\ \sin(\theta) & \cos(\theta) & 0 \\ 0 & 0 & 1 \end{bmatrix}
	 * \f]
	 *
	 * \param  angle	 The rotation angle.
	 * \return R	     The rotation matrix.
	 */
	EXPORT_SYM  Eigen::Matrix3d RotationZ(double angle);
}