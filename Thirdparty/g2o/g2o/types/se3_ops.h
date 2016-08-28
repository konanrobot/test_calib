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

#ifndef G2O_MATH_STUFF
#define G2O_MATH_STUFF

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace g2o {
  using namespace Eigen;

  inline Matrix3d skew(const Vector3d&v);
  inline Vector3d deltaR(const Matrix3d& R);
  inline Vector2d project(const Vector3d&);
  inline Vector3d project(const Vector4d&);
  inline Vector3d unproject(const Vector2d&);
  inline Vector4d unproject(const Vector3d&);
  template <typename Derived, typename DerivedOther>
  inline void skew(Eigen::MatrixBase<Derived>& s, const Eigen::MatrixBase<DerivedOther>& v){
    const double x=2*v(0);
    const double y=2*v(1);
    const double z=2*v(2);
    s <<  0.,  z, -y, -z,  0,  x,  y, -x,  0;
  }
  template <typename Derived, typename DerivedOther>
  void skew(Eigen::MatrixBase<Derived>& Sx,
      Eigen::MatrixBase<Derived>& Sy,
      Eigen::MatrixBase<Derived>& Sz,
      const Eigen::MatrixBase<DerivedOther>& R){
    const double
      r11=2*R(0,0), r12=2*R(0,1), r13=2*R(0,2),
      r21=2*R(1,0), r22=2*R(1,1), r23=2*R(1,2),
      r31=2*R(2,0), r32=2*R(2,1), r33=2*R(2,2);
    Sx <<    0,    0,    0,  -r31, -r32, -r33,   r21,   r22,  r23;
    Sy <<  r31,  r32,  r33,     0,    0,    0,  -r11,  -r12, -r13;
    Sz << -r21, -r22, -r23,   r11,   r12, r13,     0,    0,    0;
  }

#include "se3_ops.hpp"

}

#endif //MATH_STUFF
