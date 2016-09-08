#pragma once

#include <iostream>
#include <Eigen/Core>

namespace ORB_SLAM2 {

/// Parameters for pre-integration:
/// Usage: Create just a single Params and pass a shared pointer to the constructor
struct PreintegratedRotationParams {
  Eigen::Matrix3d gyroscopeCovariance;  ///< continuous-time "Covariance" of gyroscope measurements
  Eigen::Matrix4d body_P_sensor;    ///< The pose of the sensor in the body frame

  PreintegratedRotationParams() : gyroscopeCovariance(Eigen::Matrix3d::Identity()) {}

  void setGyroscopeCovariance(const Eigen::Matrix3d& cov)   { gyroscopeCovariance = cov;  }
  void setBodyPSensor(const Eigen::Matrix4d& pose)            { body_P_sensor = pose;  }

  const Eigen::Matrix3d& getGyroscopeCovariance()     const { return gyroscopeCovariance; }
  Eigen::Matrix4d   getBodyPSensor()   const { return body_P_sensor; }
};

/**
 * PreintegratedRotation is the base class for all PreintegratedMeasurements
 * classes (in AHRSFactor, ImuFactor, and CombinedImuFactor).
 * It includes the definitions of the preintegrated rotation.
 */
class PreintegratedRotation {
 public:
  typedef PreintegratedRotationParams Params;

 protected:
  /// Parameters
  Params p_;

  double deltaTij_;           ///< Time interval from i to j
  Eigen::Matrix3d deltaRij_;             ///< Preintegrated relative orientation (in frame i)
  Eigen::Matrix3d delRdelBiasOmega_;  ///< Jacobian of preintegrated rotation w.r.t. angular rate bias

  /// Default constructor for serialization
  PreintegratedRotation() {}

 public:
  /// @name Constructors
  /// @{

  /// Default constructor, resets integration to zero
  explicit PreintegratedRotation(const Params& p) : p_(p) {
    resetIntegration();
  }

  /// Explicit initialization of all class members
  PreintegratedRotation(const Params& p,
                        double deltaTij, const Eigen::Matrix3d& deltaRij,
                        const Eigen::Matrix3d& delRdelBiasOmega)
      : p_(p), deltaTij_(deltaTij), deltaRij_(deltaRij), delRdelBiasOmega_(delRdelBiasOmega) {}

  /// @}

  /// @name Basic utilities
  /// @{

  /// Re-initialize PreintegratedMeasurements
  void resetIntegration();

  /// @name Access instance variables
  /// @{
  const Params& params() const {
    return p_;
  }
  const double& deltaTij() const {
    return deltaTij_;
  }
  const Eigen::Matrix3d& deltaRij() const {
    return deltaRij_;
  }
  const Eigen::Matrix3d& delRdelBiasOmega() const {
    return delRdelBiasOmega_;
  }
  /// @}

  /// @name Main functionality
  /// @{

  /// Take the gyro measurement, correct it using the (constant) bias estimate
  /// and possibly the sensor pose, and then integrate it forward in time to yield
  /// an incremental rotation.
  Eigen::Matrix3d incrementalRotation(const Eigen::Vector3d& measuredOmega, const Eigen::Vector3d& biasHat, double deltaT,
                           OptionalJacobian<3, 3> D_incrR_integratedOmega) const;

  /// Calculate an incremental rotation given the gyro measurement and a time interval,
  /// and update both deltaTij_ and deltaRij_.
  void integrateMeasurement(const Eigen::Vector3d& measuredOmega, const Eigen::Vector3d& biasHat, double deltaT,
                            OptionalJacobian<3, 3> D_incrR_integratedOmega = boost::none,
                            OptionalJacobian<3, 3> F = boost::none);

  /// Return a bias corrected version of the integrated rotation, with optional Jacobian
  Eigen::Matrix3d biascorrectedDeltaRij(const Eigen::Vector3d& biasOmegaIncr,
                             OptionalJacobian<3, 3> H = boost::none) const;


  /// @}

};

}  /// namespace gtsam
