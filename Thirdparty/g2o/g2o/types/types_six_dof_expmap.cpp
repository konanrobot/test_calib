// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "types_six_dof_expmap.h"

#include "../core/factory.h"
#include "../stuff/macros.h"

namespace g2o {

using namespace std;


Vector2d project2d(const Vector3d& v)  {
  Vector2d res;
  res(0) = v(0)/v(2);
  res(1) = v(1)/v(2);
  return res;
}

Vector3d unproject2d(const Vector2d& v)  {
  Vector3d res;
  res(0) = v(0);
  res(1) = v(1);
  res(2) = 1;
  return res;
}

VertexSE3Expmap::VertexSE3Expmap() : BaseVertex<6, SE3Quat>() {
}

bool VertexSE3Expmap::read(std::istream& is) {
  Vector7d est;
  for (int i=0; i<7; i++)
    is  >> est[i];
  SE3Quat cam2world;
  cam2world.fromVector(est);
  setEstimate(cam2world.inverse());
  return true;
}

bool VertexSE3Expmap::write(std::ostream& os) const {
  SE3Quat cam2world(estimate().inverse());
  for (int i=0; i<7; i++)
    os << cam2world[i] << " ";
  return os.good();
}

VertexBias::VertexBias() : BaseVertex<3, Vector3d>() {
}

bool VertexBias::read(std::istream& is)
{
  (void) is;
  return false;
}

bool VertexBias::write(std::ostream& os) const
{
  (void) os;
  return false;
}

EdgeSE3ProjectXYZ::EdgeSE3ProjectXYZ() : BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexSE3Expmap>() {
}

bool EdgeSE3ProjectXYZ::read(std::istream& is){
  for (int i=0; i<2; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectXYZ::write(std::ostream& os) const {

  for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


void EdgeSE3ProjectXYZ::linearizeOplus() {
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
  SE3Quat T(vj->estimate());
  VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
  Vector3d xyz = vi->estimate();
  Vector3d xyz_trans = T.map(xyz);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double z = xyz_trans[2];
  double z_2 = z*z;

  Matrix<double,2,3> tmp;
  tmp(0,0) = fx;
  tmp(0,1) = 0;
  tmp(0,2) = -x/z*fx;

  tmp(1,0) = 0;
  tmp(1,1) = fy;
  tmp(1,2) = -y/z*fy;

  _jacobianOplusXi =  -1./z * tmp * T.rotation().toRotationMatrix();

  _jacobianOplusXj(0,0) =  x*y/z_2 *fx;
  _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *fx;
  _jacobianOplusXj(0,2) = y/z *fx;
  _jacobianOplusXj(0,3) = -1./z *fx;
  _jacobianOplusXj(0,4) = 0;
  _jacobianOplusXj(0,5) = x/z_2 *fx;

  _jacobianOplusXj(1,0) = (1+y*y/z_2) *fy;
  _jacobianOplusXj(1,1) = -x*y/z_2 *fy;
  _jacobianOplusXj(1,2) = -x/z *fy;
  _jacobianOplusXj(1,3) = 0;
  _jacobianOplusXj(1,4) = -1./z *fy;
  _jacobianOplusXj(1,5) = y/z_2 *fy;
}

Vector2d EdgeSE3ProjectXYZ::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}


Vector3d EdgeStereoSE3ProjectXYZ::cam_project(const Vector3d & trans_xyz, const float &bf) const{
  const float invz = 1.0f/trans_xyz[2];
  Vector3d res;
  res[0] = trans_xyz[0]*invz*fx + cx;
  res[1] = trans_xyz[1]*invz*fy + cy;
  res[2] = res[0] - bf*invz;
  return res;
}

EdgeStereoSE3ProjectXYZ::EdgeStereoSE3ProjectXYZ() : BaseBinaryEdge<3, Vector3d, VertexSBAPointXYZ, VertexSE3Expmap>() {
}

bool EdgeStereoSE3ProjectXYZ::read(std::istream& is){
  for (int i=0; i<=3; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeStereoSE3ProjectXYZ::write(std::ostream& os) const {

  for (int i=0; i<=3; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

void EdgeStereoSE3ProjectXYZ::linearizeOplus() {
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
  SE3Quat T(vj->estimate());
  VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
  Vector3d xyz = vi->estimate();
  Vector3d xyz_trans = T.map(xyz);

  const Matrix3d R =  T.rotation().toRotationMatrix();

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double z = xyz_trans[2];
  double z_2 = z*z;

  _jacobianOplusXi(0,0) = -fx*R(0,0)/z+fx*x*R(2,0)/z_2;
  _jacobianOplusXi(0,1) = -fx*R(0,1)/z+fx*x*R(2,1)/z_2;
  _jacobianOplusXi(0,2) = -fx*R(0,2)/z+fx*x*R(2,2)/z_2;

  _jacobianOplusXi(1,0) = -fy*R(1,0)/z+fy*y*R(2,0)/z_2;
  _jacobianOplusXi(1,1) = -fy*R(1,1)/z+fy*y*R(2,1)/z_2;
  _jacobianOplusXi(1,2) = -fy*R(1,2)/z+fy*y*R(2,2)/z_2;

  _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*R(2,0)/z_2;
  _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)-bf*R(2,1)/z_2;
  _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2)-bf*R(2,2)/z_2;

  _jacobianOplusXj(0,0) =  x*y/z_2 *fx;
  _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *fx;
  _jacobianOplusXj(0,2) = y/z *fx;
  _jacobianOplusXj(0,3) = -1./z *fx;
  _jacobianOplusXj(0,4) = 0;
  _jacobianOplusXj(0,5) = x/z_2 *fx;

  _jacobianOplusXj(1,0) = (1+y*y/z_2) *fy;
  _jacobianOplusXj(1,1) = -x*y/z_2 *fy;
  _jacobianOplusXj(1,2) = -x/z *fy;
  _jacobianOplusXj(1,3) = 0;
  _jacobianOplusXj(1,4) = -1./z *fy;
  _jacobianOplusXj(1,5) = y/z_2 *fy;

  _jacobianOplusXj(2,0) = _jacobianOplusXj(0,0)-bf*y/z_2;
  _jacobianOplusXj(2,1) = _jacobianOplusXj(0,1)+bf*x/z_2;
  _jacobianOplusXj(2,2) = _jacobianOplusXj(0,2);
  _jacobianOplusXj(2,3) = _jacobianOplusXj(0,3);
  _jacobianOplusXj(2,4) = 0;
  _jacobianOplusXj(2,5) = _jacobianOplusXj(0,5)-bf/z_2;
}


//Only Pose

bool EdgeSE3ProjectXYZOnlyPose::read(std::istream& is){
  for (int i=0; i<2; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectXYZOnlyPose::write(std::ostream& os) const {

  for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


void EdgeSE3ProjectXYZOnlyPose::linearizeOplus() {
  VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
  Vector3d xyz_trans = vi->estimate().map(Xw);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double invz = 1.0/xyz_trans[2];
  double invz_2 = invz*invz;

  _jacobianOplusXi(0,0) =  x*y*invz_2 *fx;
  _jacobianOplusXi(0,1) = -(1+(x*x*invz_2)) *fx;
  _jacobianOplusXi(0,2) = y*invz *fx;
  _jacobianOplusXi(0,3) = -invz *fx;
  _jacobianOplusXi(0,4) = 0;
  _jacobianOplusXi(0,5) = x*invz_2 *fx;

  _jacobianOplusXi(1,0) = (1+y*y*invz_2) *fy;
  _jacobianOplusXi(1,1) = -x*y*invz_2 *fy;
  _jacobianOplusXi(1,2) = -x*invz *fy;
  _jacobianOplusXi(1,3) = 0;
  _jacobianOplusXi(1,4) = -invz *fy;
  _jacobianOplusXi(1,5) = y*invz_2 *fy;
}

Vector2d EdgeSE3ProjectXYZOnlyPose::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}


Vector3d EdgeStereoSE3ProjectXYZOnlyPose::cam_project(const Vector3d & trans_xyz) const{
  const float invz = 1.0f/trans_xyz[2];
  Vector3d res;
  res[0] = trans_xyz[0]*invz*fx + cx;
  res[1] = trans_xyz[1]*invz*fy + cy;
  res[2] = res[0] - bf*invz;
  return res;
}


bool EdgeStereoSE3ProjectXYZOnlyPose::read(std::istream& is){
  for (int i=0; i<=3; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeStereoSE3ProjectXYZOnlyPose::write(std::ostream& os) const {

  for (int i=0; i<=3; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

void EdgeStereoSE3ProjectXYZOnlyPose::linearizeOplus() {
  VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
  Vector3d xyz_trans = vi->estimate().map(Xw);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double invz = 1.0/xyz_trans[2];
  double invz_2 = invz*invz;

  _jacobianOplusXi(0,0) =  x*y*invz_2 *fx;
  _jacobianOplusXi(0,1) = -(1+(x*x*invz_2)) *fx;
  _jacobianOplusXi(0,2) = y*invz *fx;
  _jacobianOplusXi(0,3) = -invz *fx;
  _jacobianOplusXi(0,4) = 0;
  _jacobianOplusXi(0,5) = x*invz_2 *fx;

  _jacobianOplusXi(1,0) = (1+y*y*invz_2) *fy;
  _jacobianOplusXi(1,1) = -x*y*invz_2 *fy;
  _jacobianOplusXi(1,2) = -x*invz *fy;
  _jacobianOplusXi(1,3) = 0;
  _jacobianOplusXi(1,4) = -invz *fy;
  _jacobianOplusXi(1,5) = y*invz_2 *fy;

  _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*y*invz_2;
  _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)+bf*x*invz_2;
  _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2);
  _jacobianOplusXi(2,3) = _jacobianOplusXi(0,3);
  _jacobianOplusXi(2,4) = 0;
  _jacobianOplusXi(2,5) = _jacobianOplusXi(0,5)-bf*invz_2;
}

EdgeSE3Calib::EdgeSE3Calib()
{
}

bool EdgeSE3Calib::read(std::istream& is)
{
  (void) is;
  return false;
}

bool EdgeSE3Calib::write(std::ostream& os) const
{
  (void) os;
  return false;
}

//void EdgeSE3Calib::computeError()
//{
//  // K_T_G previous camera pose
////  const VertexSE3Expmap* lastCamNode = static_cast<const VertexSE3Expmap*>(_vertices[0]);
//  // K+1_T_G current camera pose
//  const VertexSE3Expmap* curCamNode = static_cast<const VertexSE3Expmap*>(_vertices[0]);
//  // Calib offset between camera and odometry
//  const VertexSE3Expmap* offsetNode = static_cast<const VertexSE3Expmap*>(_vertices[1]);
//  // delta T between two camera
////  const SE3Quat& Tcl = curCamNode->estimate() * lastCamNode->estimate().inverse();

//  // measurement = K+1_T_k, transformation between last and current odometry frame.
//  // offsetNode is the estimation of calibration between stereo and odom.
//  SE3Quat temp = measurement().inverse()*(offsetNode->estimate().inverse()*curCamNode->estimate()*offsetNode->estimate());
//  _error = temp.log();
//}

void EdgeSE3Calib::computeError()
{
  // K_T_G previous camera pose
  const VertexSE3Expmap* lastCamNode = static_cast<const VertexSE3Expmap*>(_vertices[0]);
  // K+1_T_G current camera pose
  const VertexSE3Expmap* curCamNode = static_cast<const VertexSE3Expmap*>(_vertices[1]);
  // Calib offset between camera and odometry
  const VertexSE3Expmap* offsetNode = static_cast<const VertexSE3Expmap*>(_vertices[2]);

//  Eigen::Quaterniond lastQ = lastCamNode->estimate().rotation();
//  Eigen::Quaterniond curQ = curCamNode->estimate().rotation();
//  Eigen::Quaterniond offsetQ = offsetNode->estimate().rotation();
//  Eigen::Quaterniond odomQ = measurement().rotation();

//  Eigen::Quaterniond temp1 = LeftTimes(inverseJPL(offsetQ),curQ);
//  Eigen::Quaterniond temp2 = LeftTimes(temp1,inverseJPL(lastQ));
//  Eigen::Quaterniond temp3 = LeftTimes(temp2,offsetQ);
//  Eigen::Quaterniond result = LeftTimes(temp3,inverseJPL(odomQ));

//  Eigen::Vector3d delta_t = getVectorJPL(result);

//  Eigen::Matrix3d cc = toRotationMatJPL(curQ)*toRotationMatJPL(lastQ).transpose();
//  Eigen::Matrix3d co = toRotationMatJPL(offsetQ);
//  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
//  Eigen::Vector3d pco = offsetNode->estimate().translation();
//  Eigen::Vector3d pc1 = curCamNode->estimate().translation();
//  Eigen::Vector3d pc0 = lastCamNode->estimate().translation();
//  Eigen::Vector3d poo = measurement().translation();
//  Eigen::Vector3d vresult = co.transpose()*((cc-I)*pco+pc1-cc*pc0)-poo;

//  g2o::Vector6d temp;
//  for(int i = 0; i < 3; i++)
//      temp[i] = delta_t[i];
//  for(int i = 0; i < 3; i++)
//      temp[i+3] = vresult[i];

//  _error = temp;

  // delta T between two camera
  const SE3Quat& Tcl = curCamNode->estimate() * lastCamNode->estimate().inverse();

//   measurement = K+1_T_k, transformation between last and current odometry frame.
//   offsetNode is the estimation of calibration between stereo and odom.
  SE3Quat temp = measurement().inverse()*(offsetNode->estimate().inverse()*Tcl*offsetNode->estimate());
  _error = temp.log();
}

EdgePreint::EdgePreint()
{
}

bool EdgePreint::read(std::istream& is)
{
  (void) is;
  return false;
}

bool EdgePreint::write(std::ostream& os) const
{
  (void) os;
  return false;
}

void EdgePreint::computeError()
{
  // last camera pose
  const VertexSE3Expmap* lastCamNode = static_cast<const VertexSE3Expmap*>(_vertices[0]);
  // current camera pose
  const VertexSE3Expmap* curCamNode = static_cast<const VertexSE3Expmap*>(_vertices[1]);
  // bias node
  const VertexBias* biasNode = static_cast<const VertexBias*>(_vertices[2]);

//  GyroPreint meas = measurement();
//  Eigen::Quaterniond q_est = meas.estimation();
//  Eigen::Matrix3d jaccobian = meas.jacobian();
//  bias_bar = meas.bias();

  Eigen::Quaterniond k_1_q = curCamNode->estimate().rotation();
  Eigen::Quaterniond k_q = lastCamNode->estimate().rotation();
  Eigen::Quaterniond k_q_T = k_q.inverse();
  Eigen::Vector3d bias = biasNode->estimate();

  Eigen::Quaterniond q1 = LeftTimes(k_1_q, k_q_T);
  Eigen::Quaterniond q2 = LeftTimes(q1, convertJPL(measurement()).inverse());
  Eigen::Vector3d biasv = 0.5*jaccobian*(bias-bias_bar);
  Eigen::Quaterniond q3;
  q3.x() = biasv[0]; q3.y() = biasv[1]; q3.z() = biasv[2];
  q3.w() = 1;
  Eigen::Quaterniond q4 = LeftTimes(q2, q3);
  Eigen::Vector3d delta_theta = getVectorJPL(q4);
  Eigen::Vector3d e;
  e << delta_theta;

  _error = e;
}

//void EdgePreint::linearizeOplus()
//{
//    VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
//    VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
//    VertexBias * vk = static_cast<VertexBias *>(_vertices[2]);

//    Eigen::Quaterniond k_1_q = vj->estimate().rotation();
//    Eigen::Quaterniond k_q = vi->estimate().rotation();
//    Eigen::Vector3d bias = vk->estimate();

//    Eigen::Quaterniond q_est = convertJPL(measurement());

//    Eigen::Quaterniond q1 = LeftTimes(k_1_q, k_q.inverse());
//    Eigen::Quaterniond q2 = LeftTimes(q1, q_est.inverse());
//    Eigen::Vector3d biasv = 0.5*jaccobian*(bias-bias_bar);
//    Eigen::Quaterniond q3;
//    q3.x() = biasv[0]; q3.y() = biasv[1]; q3.z() = biasv[2];
//    q3.w() = 1;
//    Eigen::Quaterniond q4 = LeftTimes(q2, q3);

//    qr = q2;
//    qrb = q4;
//    qn = q1;
//    qmb = measurement().estimation();

//    Eigen::Matrix3d tmp = qrb.w()*Eigen::Matrix3d::Identity()
//            + g2o::skewJPL(g2o::getVectorJPL(qrb));

//    _jacobianOplus[0](0,0) = tmp(0,0);
//    _jacobianOplus[0](0,1) = tmp(0,1);
//    _jacobianOplus[0](0,2) = tmp(0,2);
//    _jacobianOplus[0](0,3) = 0;
//    _jacobianOplus[0](0,4) = 0;
//    _jacobianOplus[0](0,5) = 0;

//    _jacobianOplus[0](1,0) = tmp(1,0);
//    _jacobianOplus[0](1,1) = tmp(1,1);
//    _jacobianOplus[0](1,2) = tmp(1,2);
//    _jacobianOplus[0](1,3) = 0;
//    _jacobianOplus[0](1,4) = 0;
//    _jacobianOplus[0](1,5) = 0;

//    _jacobianOplus[0](2,0) = tmp(2,0);
//    _jacobianOplus[0](2,1) = tmp(2,1);
//    _jacobianOplus[0](2,2) = tmp(2,2);
//    _jacobianOplus[0](2,3) = 0;
//    _jacobianOplus[0](2,4) = 0;
//    _jacobianOplus[0](2,5) = 0;

//    Eigen::Matrix3d tmp1 = (qn.w()*Eigen::Matrix3d::Identity()
//            - g2o::skewJPL(g2o::getVectorJPL(qn)))
//            *(qmb.w()*Eigen::Matrix3d::Identity()
//              - g2o::skewJPL(g2o::getVectorJPL(qmb)))
//            +g2o::getVectorJPL(qn)*g2o::getVectorJPL(qmb).transpose();

//    _jacobianOplus[1](0,0) = tmp1(0,0);
//    _jacobianOplus[1](0,1) = tmp1(0,1);
//    _jacobianOplus[1](0,2) = tmp1(0,2);
//    _jacobianOplus[1](0,3) = 0;
//    _jacobianOplus[1](0,4) = 0;
//    _jacobianOplus[1](0,5) = 0;

//    _jacobianOplus[1](1,0) = tmp1(1,0);
//    _jacobianOplus[1](1,1) = tmp1(1,1);
//    _jacobianOplus[1](1,2) = tmp1(1,2);
//    _jacobianOplus[1](1,3) = 0;
//    _jacobianOplus[1](1,4) = 0;
//    _jacobianOplus[1](1,5) = 0;

//    _jacobianOplus[1](2,0) = tmp1(2,0);
//    _jacobianOplus[1](2,1) = tmp1(2,1);
//    _jacobianOplus[1](2,2) = tmp1(2,2);
//    _jacobianOplus[1](2,3) = 0;
//    _jacobianOplus[1](2,4) = 0;
//    _jacobianOplus[1](2,5) = 0;

//    Eigen::Matrix3d temp = (qr.w()*Eigen::Matrix3d::Identity()
//            - g2o::skewJPL(g2o::getVectorJPL(qr)))*measurement().jacobian();

//    _jacobianOplus[2](0,0) = temp(0,0);
//    _jacobianOplus[2](0,1) = temp(0,1);
//    _jacobianOplus[2](0,2) = temp(0,2);

//    _jacobianOplus[2](1,0) = temp(1,0);
//    _jacobianOplus[2](1,1) = temp(1,1);
//    _jacobianOplus[2](1,2) = temp(1,2);

//    _jacobianOplus[2](2,0) = temp(2,0);
//    _jacobianOplus[2](2,1) = temp(2,1);
//    _jacobianOplus[2](2,2) = temp(2,2);
//}

//ParameterSE3Offset::ParameterSE3Offset(){
//  setOffset();
//}

//void ParameterSE3Offset::setOffset(const SE3Quat& offset_){
//  _offset = offset_;
//  _inverseOffset = _offset.inverse();
//}

//bool ParameterSE3Offset::read(std::istream& is) {
//  return false;
//}

//bool ParameterSE3Offset::write(std::ostream& os) const {
//  return false;
//}

//CacheSE3Offset::CacheSE3Offset() :
//  Cache(),
//  _offsetParam(0)
//{
//}

//bool CacheSE3Offset::resolveDependancies(){
//  _offsetParam = dynamic_cast <ParameterSE3Offset*> (_parameters[0]);
//  return _offsetParam != 0;
//}

//void CacheSE3Offset::updateImpl(){
//  const VertexSE3Expmap* v = static_cast<const VertexSE3Expmap*>(vertex());
//  _n2w = v->estimate() * _offsetParam->offset();
//  _w2n = _n2w.inverse();
//  _w2l = v->estimate().inverse();
//}

//void CacheSE3Offset::setOffsetParam(ParameterSE3Offset* offsetParam)
//{
//  _offsetParam = offsetParam;
//}

//// point to camera projection, monocular
//EdgeSE3Prior::EdgeSE3Prior() : BaseUnaryEdge<6, SE3Quat, VertexSE3Expmap>() {
//  setMeasurement(SE3Quat(Quaterniond::Identity(), Eigen::Vector3d()));
//  information().setIdentity();
//  _cache = 0;
//  _offsetParam = 0;
//  resizeParameters(1);
//  installParameter(_offsetParam, 0);
//}


//bool EdgeSE3Prior::resolveCaches(){
//  assert(_offsetParam);
//  ParameterVector pv(1);
//  pv[0]=_offsetParam;
//  resolveCache(_cache, (OptimizableGraph::Vertex*)_vertices[0],"CACHE_SE3_OFFSET",pv);
//  return _cache != 0;
//}



//bool EdgeSE3Prior::read(std::istream& is) {
//  return false;
//}

//bool EdgeSE3Prior::write(std::ostream& os) const {
//  return false;
//}


//void EdgeSE3Prior::computeError() {
//  SE3Quat delta=_inverseMeasurement * _cache->n2w();
//  _error = delta.log();
//}

//template <typename Derived>
//void computeEdgeSE3PriorGradient(SE3Quat& E,
//                                 const Eigen::MatrixBase<Derived>& JConstRef,
//                                 const SE3Quat& Z,
//                                 const SE3Quat& X,
//                                 const SE3Quat& P=SE3Quat())
//{
//  Eigen::MatrixBase<Derived>& J = const_cast<Eigen::MatrixBase<Derived>&>(JConstRef);
//  J.derived().resize(6,6);
//  // compute the error at the linearization point
//  const SE3Quat A = Z.inverse()*X;
//  const SE3Quat& B = P;
//  Eigen::Matrix3d Ra = A.rotation().toRotationMatrix();
//  Eigen::Matrix3d Rb = B.rotation().toRotationMatrix();
//  Eigen::Vector3d tb = B.translation();
//  E = A*B;
//  Eigen::Matrix3d Re = E.rotation().toRotationMatrix();

//  Eigen::Matrix<double, 3, 9, Eigen::ColMajor> dq_dR;
//  compute_dq_dR (dq_dR,
//      Re(0,0),Re(1,0),Re(2,0),
//      Re(0,1),Re(1,1),Re(2,1),
//      Re(0,2),Re(1,2),Re(2,2));

//  J.setZero();

//  // dte/dt
//  J.template block<3,3>(0,0)=Ra;

//  // dte/dq =0
//  // dte/dqj
//  {
//    Eigen::Matrix3d S;
//    skew(S,tb);
//    J.template block<3,3>(0,3)=Ra*S;
//  }

//  // dre/dt =0

//  // dre/dq
//  {
//    double buf[27];
//    Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::ColMajor> > M(buf);
//    Eigen::Matrix3d Sx,Sy,Sz;
//    skew(Sx,Sy,Sz,Rb);
//#ifdef __clang__
//    Matrix3D temp = Ra * Sx;
//    Eigen::Map<Matrix3D> M2(temp.data());
//    Eigen::Map<Matrix3D> Mx(buf);    Mx = M2;
//    temp = Ra*Sy;
//    Eigen::Map<Matrix3D> My(buf+9);  My = M2;
//    temp = Ra*Sz;
//    Eigen::Map<Matrix3D> Mz(buf+18); Mz = M2;
//#else
//    Eigen::Map<Eigen::Matrix3d> Mx(buf);    Mx = Ra*Sx;
//    Eigen::Map<Eigen::Matrix3d> My(buf+9);  My = Ra*Sy;
//    Eigen::Map<Eigen::Matrix3d> Mz(buf+18); Mz = Ra*Sz;
//#endif
//    J.template block<3,3>(3,3) = dq_dR * M;
//  }

//}

//void EdgeSE3Prior::linearizeOplus(){
//  VertexSE3Expmap *from = static_cast<VertexSE3Expmap*>(_vertices[0]);
//  SE3Quat E;
//  SE3Quat Z, X, P;
//  X=from->estimate();
//  P=_cache->offsetParam()->offset();
//  Z=_measurement;
//  computeEdgeSE3PriorGradient(E, _jacobianOplusXi, Z, X, P);
//}


//bool EdgeSE3Prior::setMeasurementFromState(){
//  setMeasurement(_cache->n2w());
//  return true;
//}


//void EdgeSE3Prior::initialEstimate(const OptimizableGraph::VertexSet& /*from_*/, OptimizableGraph::Vertex* /*to_*/) {
//  VertexSE3Expmap *v = static_cast<VertexSE3Expmap*>(_vertices[0]);
//  assert(v && "Vertex for the Prior edge is not set");

//  SE3Quat newEstimate = _offsetParam->offset().inverse() * measurement();
//  if (_information.block<3,3>(0,0).array().abs().sum() == 0){ // do not set translation, as that part of the information is all zero
//    newEstimate.setTranslation(v->estimate().translation());
//  }
//  if (_information.block<3,3>(3,3).array().abs().sum() == 0){ // do not set rotation, as that part of the information is all zero
//    newEstimate.setRotation(v->estimate().rotation());
//  }
//  v->setEstimate(newEstimate);
//}

} // end namespace
