#ifndef GYRO_PREINT_H
#define GYRO_PREINT_H

#include "JPL.hpp"

namespace g2o {
    using namespace Eigen;
    class GyroPreint
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    protected:
        Quaterniond _estimate;
//        Matrix3d _covariance;
//        Matrix3d _Jacobian;
//        Vector3d _bias;
    public:
        GyroPreint(){
          _estimate.setIdentity();
//          _covariance.setZero();
//          _Jacobian.setZero();
//          _bias.setZero();
        }

        inline const Quaterniond& estimation() const {return _estimate;}

        inline void setEstimation(const Quaterniond& estimate_) {_estimate = estimate_;}

//        inline const Matrix3d& covariance() const {return _covariance;}

//        inline void setCovariance(const Matrix3d& covariance_) {_covariance=covariance_;}

//        inline const Matrix3d& jacobian() const {return _Jacobian;}

//        inline void setJacobian(const Matrix3d& Jacobian_) {_Jacobian=Jacobian_;}

//        inline void setBias(const Vector3d& bias_) {_bias = bias_;}

//        inline const Vector3d& bias() const {return _bias;}

    };
}

#endif //gyropreint.h
