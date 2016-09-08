#include "PreintegrateRotation.h"

using namespace std;

namespace ORB_SLAM2 {

void PreintegratedRotation::resetIntegration() {
  deltaTij_ = 0.0;
  deltaRij_ = Eigen::Matrix3d::Zero();
  delRdelBiasOmega_ = Eigen::Matrix3d::Zero();
}

Rot3 PreintegratedRotation::incrementalRotation(const Eigen::Vector3d& measuredOmega,
    const Eigen::Vector3d& biasHat, double deltaT,
    OptionalJacobian<3, 3> D_incrR_integratedOmega) const {

  // First we compensate the measurements for the bias
  Eigen::Vector3d correctedOmega = measuredOmega - biasHat;

  // Then compensate for sensor-body displacement: we express the quantities
  // (originally in the IMU frame) into the body frame
    Eigen::Matrix3d body_R_sensor = p_.body_P_sensor.block<3,3>(0,0);
    // rotation rate vector in the body frame
    correctedOmega = body_R_sensor * correctedOmega;

  // rotation vector describing rotation increment computed from the
  // current rotation rate measurement
  const Eigen::Vector3d integratedOmega = correctedOmega * deltaT;
  return Rot3::Expmap(integratedOmega, D_incrR_integratedOmega); // expensive !!
}

void PreintegratedRotation::integrateMeasurement(const Vector3& measuredOmega,
    const Vector3& biasHat, double deltaT,
    OptionalJacobian<3, 3> optional_D_incrR_integratedOmega,
    OptionalJacobian<3, 3> F) {
  Matrix3 D_incrR_integratedOmega;
  const Rot3 incrR = incrementalRotation(measuredOmega, biasHat, deltaT,
      D_incrR_integratedOmega);

  // If asked, pass first derivative as well
  if (optional_D_incrR_integratedOmega) {
    *optional_D_incrR_integratedOmega << D_incrR_integratedOmega;
  }

  // Update deltaTij and rotation
  deltaTij_ += deltaT;
  deltaRij_ = deltaRij_.compose(incrR, F);

  // Update Jacobian
  const Matrix3 incrRt = incrR.transpose();
  delRdelBiasOmega_ = incrRt * delRdelBiasOmega_
      - D_incrR_integratedOmega * deltaT;
}

Rot3 PreintegratedRotation::biascorrectedDeltaRij(const Vector3& biasOmegaIncr,
    OptionalJacobian<3, 3> H) const {
  const Vector3 biasInducedOmega = delRdelBiasOmega_ * biasOmegaIncr;
  const Rot3 deltaRij_biascorrected = deltaRij_.expmap(biasInducedOmega,
      boost::none, H);
  if (H)
    (*H) *= delRdelBiasOmega_;
  return deltaRij_biascorrected;
}

//Vector3 PreintegratedRotation::integrateCoriolis(const Rot3& rot_i) const {
//  if (!p_->omegaCoriolis)
//    return Vector3::Zero();
//  return rot_i.transpose() * (*p_->omegaCoriolis) * deltaTij_;
//}

} // namespace gtsam
