#ifndef JPL_H_
#define JPL_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>

namespace ORB_SLAM2 {

/// Wrapper to convert the Hamilton quaternion to JPL quaternion
/// @anchor Doom

/// Returns a matrix with angular velocities used for quaternion derivatives/
// integration with the JPL notation.
/**
 The quaternion to be multiplied with this matrix has to be in the order x y z w !!!
 \param vec 3D vector with angular velocities.
 \return 4x4 matrix for multiplication with the quaternion.
 */
template<class Derived>
inline Eigen::Matrix<typename Derived::Scalar, 4, 4> OmegaMatJPL(
    const Eigen::MatrixBase<Derived> & vec) {
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);
  return (Eigen::Matrix<typename Derived::Scalar, 4, 4>() << 0, vec[2], -vec[1],
      vec[0], -vec[2], 0, vec[0], vec[1], vec[1], -vec[0], 0, vec[2], -vec[0],
      -vec[1], -vec[2], 0)
      .finished();
}

}

#endif //JPH.h
