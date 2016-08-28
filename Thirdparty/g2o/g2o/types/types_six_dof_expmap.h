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

// Modified by Raúl Mur Artal (2014)
// Added EdgeSE3ProjectXYZ (project using focal_length in x,y directions)
// Modified by Raúl Mur Artal (2016)
// Added EdgeStereoSE3ProjectXYZ (project using focal_length in x,y directions)
// Added EdgeSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)
// Added EdgeStereoSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)

#ifndef G2O_SIX_DOF_TYPES_EXPMAP
#define G2O_SIX_DOF_TYPES_EXPMAP

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_multi_edge.h"
#include "../core/base_unary_edge.h"
#include "../core/cache.h"
#include "se3_ops.h"
#include "se3quat.h"
#include "dquat2mat.h"
#include "types_sba.h"
#include "gyropreint.h"
#include <Eigen/Geometry>

namespace g2o {
namespace types_six_dof_expmap {
void init();
}

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;

/**
 * \brief SE3 Vertex parameterized internally with a transformation matrix
 and externally with its exponential map
 */
class  VertexSE3Expmap : public BaseVertex<6, SE3Quat>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  VertexSE3Expmap();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  virtual void setToOriginImpl() {
    _estimate = SE3Quat();
  }

  virtual void oplusImpl(const double* update_)  {
    Eigen::Map<const Vector6d> update(update_);
//    Vector6d update;
//    update[1] = update_[1]; update[3] = update_[3]; update[5] = update_[5];
    setEstimate(SE3Quat::exp(update)*estimate());
  }
};

///
/// \brief The bias vertex class
///
class  VertexBias : public BaseVertex<3, Vector3d>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  VertexBias();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  virtual void setToOriginImpl() {
    _estimate = Vector3d();
  }

  virtual void oplusImpl(const double* update_)  {
    Eigen::Map<const Vector3d> update(update_);
    setEstimate(estimate()+update);
  }
};

class  EdgeSE3ProjectXYZ: public  BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3ProjectXYZ();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    Vector2d obs(_measurement);
    _error = obs-cam_project(v1->estimate().map(v2->estimate()));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    return (v1->estimate().map(v2->estimate()))(2)>0.0;
  }
    

  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  double fx, fy, cx, cy;
};


class  EdgeStereoSE3ProjectXYZ: public  BaseBinaryEdge<3, Vector3d, VertexSBAPointXYZ, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeStereoSE3ProjectXYZ();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    Vector3d obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(v2->estimate()),bf);
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    return (v1->estimate().map(v2->estimate()))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector3d cam_project(const Vector3d & trans_xyz, const float &bf) const;

  double fx, fy, cx, cy, bf;
};

class  EdgeSE3ProjectXYZOnlyPose: public  BaseUnaryEdge<2, Vector2d, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3ProjectXYZOnlyPose(){}

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    Vector2d obs(_measurement);
    _error = obs-cam_project(v1->estimate().map(Xw));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    return (v1->estimate().map(Xw))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  Vector3d Xw;
  double fx, fy, cx, cy;
};


class  EdgeStereoSE3ProjectXYZOnlyPose: public  BaseUnaryEdge<3, Vector3d, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeStereoSE3ProjectXYZOnlyPose(){}

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    Vector3d obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(Xw));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    return (v1->estimate().map(Xw))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector3d cam_project(const Vector3d & trans_xyz) const;

  Vector3d Xw;
  double fx, fy, cx, cy, bf;
};

//class EdgeSE3Calib : public BaseBinaryEdge<6, SE3Quat, VertexSE3Expmap, VertexSE3Expmap>
//{
//  public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//    EdgeSE3Calib();

//    void computeError();

////    void setMeasurement(const SE3Quat& m){
////        _measurement = m;
////    }

//    virtual bool read(std::istream& is);
//    virtual bool write(std::ostream& os) const;
//};

class EdgeSE3Calib : public BaseMultiEdge<6, SE3Quat>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgeSE3Calib();

    void computeError();

    void setMeasurement(const SE3Quat& m){
        _measurement = m;
    }

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
};

class EdgePreint : public BaseMultiEdge<3, Eigen::Vector3d>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgePreint();

    void computeError();

    void setMeasurement(const Eigen::Vector3d& m){
        _measurement = m;
    }

//    virtual void linearizeOplus();

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

    Eigen::Matrix3d jaccobian;
    Eigen::Vector3d bias_bar;

 private:
    Eigen::Quaterniond qr;
    Eigen::Quaterniond qrb;
    Eigen::Quaterniond qn;
    Eigen::Quaterniond qmb;
};

//class ParameterSE3Offset: public Parameter
//{
//  public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
//    ParameterSE3Offset();

//    virtual bool read(std::istream& is);
//    virtual bool write(std::ostream& os) const;

//    /**
//     * update the offset to a new value.
//     * re-calculates the different representations, e.g., the rotation matrix
//     */
//    void setOffset(const SE3Quat& offset_=SE3Quat(Quaterniond::Identity(),Eigen::Vector3d()));

//    //! rotation of the offset as 3x3 rotation matrix
//    const SE3Quat& offset() const { return _offset;}

//    //! rotation of the inverse offset as 3x3 rotation matrix
//    const SE3Quat& inverseOffset() const { return _inverseOffset;}

//  protected:
//    SE3Quat _offset;
//    SE3Quat _inverseOffset;
//};

///**
// * \brief caching the offset related to a vertex
// */
//class CacheSE3Offset: public Cache {
//  public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
//    CacheSE3Offset();
//    virtual void updateImpl();

//    const ParameterSE3Offset* offsetParam() const { return _offsetParam;}
//    void setOffsetParam(ParameterSE3Offset* offsetParam);

//    const SE3Quat& w2n() const { return _w2n;}
//    const SE3Quat& n2w() const { return _n2w;}
//    const SE3Quat& w2l() const { return _w2l;}

//  protected:
//    ParameterSE3Offset* _offsetParam; ///< the parameter connected to the cache
//    SE3Quat _w2n;
//    SE3Quat _n2w;
//    SE3Quat _w2l;

//  protected:
//    virtual bool resolveDependancies();
//};

//class EdgeSE3Prior : public BaseUnaryEdge<6, SE3Quat, VertexSE3Expmap> {
//public:
//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//  EdgeSE3Prior();
//  virtual bool read(std::istream& is);
//  virtual bool write(std::ostream& os) const;

//  // return the error estimate as a 3-vector
//  void computeError();

//  // jacobian
//  virtual void linearizeOplus();

//  virtual void setMeasurement(const SE3Quat& m){
//    _measurement = m;
//    _inverseMeasurement = m.inverse();
//  }

//  virtual bool setMeasurementData(const double* d){
//    return false;
//  }

//  virtual bool getMeasurementData(double* d) const{
//    return false;
//  }

//  virtual int measurementDimension() const {return 7;}

//  virtual bool setMeasurementFromState() ;

//  virtual double initialEstimatePossible(const OptimizableGraph::VertexSet& /*from*/,
//           OptimizableGraph::Vertex* /*to*/) {
//    return 1.;
//  }

//  virtual void initialEstimate(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* to);
//protected:
//  SE3Quat _inverseMeasurement;
//  virtual bool resolveCaches();
//  ParameterSE3Offset* _offsetParam;
//  CacheSE3Offset* _cache;
//};

} // end namespace

#endif
