#ifndef PREINTEGRATOR_H_
#define PREINTEGRATOR_H_

#include <iostream>

#include "Thirdparty/g2o/g2o/types/JPL.hpp"
#include "Thirdparty/g2o/g2o/types/gyropreint.h"

using namespace std;

namespace ORB_SLAM2
{

class Preintegrator
{
public:
    ///
    /// \brief getPreintResult get the preintegration result
    /// \param gyro gyro datas between two keyframes
    /// \param q_est JPL quaternion estimation of the rotation between two keyframes
    /// \param info the infomation matrix(covariance.inverse)
    /// \param bias_jac bias jacobian
    ///
    static void getPreintResult(const vector<pair<double, double>> &gyro,
                                const Eigen::Quaterniond &q_k,
                                const Eigen::Vector3d &bias_bar,
                                Eigen::Quaterniond &q_est,
                                Eigen::Matrix3d &info,
                                Eigen::Matrix3d &bias_jac);
};

}

#endif //Preintegrator.h
