#ifndef _PREINTEGRATION_PARAM_
#define _PREINTEGRATION_PARAM_
#include <iostream>
#include <Eigen/Core>

namespace ORB_SLAM2 {
    struct PreintegrationParams: PreintegratedRotationParams {
    //  Matrix3 accelerometerCovariance; ///< continuous-time "Covariance" of accelerometer
      Eigen::Matrix3d integrationCovariance; ///< continuous-time "Covariance" describing integration uncertainty
    //  bool use2ndOrderCoriolis; ///< Whether to use second order Coriolis integration
    //  Vector3 n_gravity; ///< Gravity vector in nav frame

    PreintegrationParams()
        :integrationCovariance(Eigen::Matrix3d::Identity()){}

    void setIntegrationCovariance(const Eigen::Matrix3d& cov)   { integrationCovariance = cov; }

    const Eigen::Matrix3d& getIntegrationCovariance()   const { return integrationCovariance; }

    protected:
      /// Default constructor for serialization only: uninitialized!
      PreintegrationParams() {}
};
}

#endif
