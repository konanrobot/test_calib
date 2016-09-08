#ifndef _IMU_BIAS_
#define _IMU_BIAS_
#include <iostream>
#include <Eigen/Core>

namespace ORB_SLAM2 {
class ConstantBias{
private:
    Eigen::Vector3d biasGyro_;
public:
    ConstantBias() :
       biasGyro_(0.0, 0.0, 0.0) {
    }

    ConstantBias(const Eigen::Vector3d& biasGyro) :
        biasGyro_(biasGyro) {
    }

    Eigen::Vector3d vector() const {
      Vector3 v;
      v << biasGyro_;
      return v;
    }

    /** get gyroscope bias */
    const Vector3& gyroscope() const {
      return biasGyro_;
    }

    /** addition of vector on right */
    ConstantBias operator+(const Eigen::Vector3d& v) const {
      return ConstantBias(biasGyro_ + v.tail<3>());
    }

    /** addition */
    ConstantBias operator+(const ConstantBias& b) const {
      return ConstantBias(biasGyro_ + b.biasGyro_);
    }

    /** subtraction */
    ConstantBias operator-(const ConstantBias& b) const {
      return ConstantBias(biasGyro_ - b.biasGyro_);
    }
};
}
#endif
