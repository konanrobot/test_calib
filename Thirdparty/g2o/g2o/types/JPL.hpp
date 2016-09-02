#ifndef JPL_H_
#define JPL_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>

/// inline JPL quaternion method using the Eigen quaternion (hamilton).
/// @author Doom
namespace g2o
{
    ///
    /// @brief get the JPL quaternion's vector
    /// @param _q JPL quaternion (use Eigen quaternion to store)
    ///
    inline Eigen::Matrix<double, 3, 1> getVectorJPL(const Eigen::Quaterniond &_q)
    {
        return (Eigen::Matrix<double, 3, 1>() << _q.x(), _q.y(),
                _q.z()).finished();
    }

    inline Eigen::Quaterniond inverseJPL(const Eigen::Quaterniond &_q)
    {
        Eigen::Quaterniond _r;
        for(int i =0; i < 3; i++)
        {
           _r.coeffs()(i) = -_q.coeffs()(i);
        }
        _r.w() = _q.w();
        return _r;
    }

    inline Eigen::Quaterniond smallAngle(const Eigen::Vector3d &_v)
    {
        Eigen::Quaterniond _r;
        for(int i =0; i < 3; i++)
        {
           _r.coeffs()(i)=0.5*_v[i];
        }
        _r.w() = 1;
        return _r;
    }

    ///
    /// @brief convert a vector to quaternion
    /// @param _v rotation vector
    ///
    inline Eigen::Quaterniond convertJPL(const Eigen::Vector3d &_v)
    {
        Eigen::Quaterniond _r;
        for(int i =0; i < 3; i++)
        {
           _r.coeffs()(i)=_v[i];
        }
        _r.w() = 0.; // recover the positive w
        if (_r.norm()>1.){
          _r.normalize();
        } else {
          double w2=1.-_r.squaredNorm();
          _r.w()= (w2<0.) ? 0. : sqrt(w2);
        }
        return _r;
    }

    ///
    /// @brief convert a small vector to quaternion
    /// @param _v rotation vector
    ///
    inline Eigen::Quaterniond convertSmallJPL(const Eigen::Vector3d &_v)
    {
        Eigen::Quaterniond _r;
        for(int i =0; i < 3; i++)
        {
           _r.coeffs()(i)=0.5*_v[i];
        }
        _r.w() = 1;
        return _r;
    }

    ///
    /// @brief get skew matrix with JPL format
    /// @param v 3*1 vector, output 3*3 matrix
    ///
    inline Eigen::Matrix<double, 3, 3> skewJPL(const Eigen::Vector3d & v)
    {
        return (Eigen::Matrix<double, 3, 3>() << 0, -v[2], v[1],
            v[2], 0, -v[0], -v[1], v[0], 0).finished();
    }

    ///
    /// @brief get Omega matrix of the JPL quaternion
    /// @param omega 3*1(1*3) rotation velocity vector
    ///
    inline Eigen::Matrix<double, 4, 4> OmegaJPL(
            const Eigen::Vector3d & omega)
    {
         return (Eigen::Matrix<double, 4, 4>() << 0, omega[2], -omega[1],
             omega[0], -omega[2], 0, omega[0], omega[1], omega[1], -omega[0], 0, omega[2], -omega[0],
             -omega[1], -omega[2], 0)
             .finished();
    }

    ///
    /// @brief left operator of JPL
    /// [q4*I3*3-skew(q)]
    /// [     -q^T      ]
    /// where quat = [q; q4]
    /// @param _q JPL quaternion, which use Eigen format
    inline Eigen::Matrix<double, 4, 3> LeftJPL(
        const Eigen::Quaterniond & _q) {
      return (Eigen::Matrix<double, 4, 3>() <<
      // This is the left matrix of JPL quaternion ---
          _q.w() * Eigen::Matrix<double, 3, 3>::Identity()
              - skewJPL(getVectorJPL(_q)), -getVectorJPL(_q).transpose()
      // ---
          ).finished();
    }

    ///
    /// @brief right operator of JPL
    /// [q4*I3*3+skew(q)]
    /// [     -q^T      ]
    /// where quat = [q; q4]
    /// @param _q JPL quaternion, which use Eigen format
    inline Eigen::Matrix<double, 4, 3> RightJPL(
        const Eigen::Quaterniond & _q) {
      return (Eigen::Matrix<double, 4, 3>() <<
      // This is the right matrix of JPL quaternion ---
          _q.w() * Eigen::Matrix<double, 3, 3>::Identity()
              + skewJPL(getVectorJPL(_q)), -getVectorJPL(_q).transpose()
      // ---
          ).finished();
    }

    ///
    /// @brief left times
    ///
    inline Eigen::Quaterniond LeftTimes(
        const Eigen::Quaterniond & q0,
        const Eigen::Quaterniond & q1)
    {
        Eigen::Matrix<double, 4, 4> left;
        Eigen::Matrix<double,4,1> q;
        q[0] = q0.x(); q[1] = q0.y();
        q[2] = q0.z(); q[3] = q0.w();
        left << LeftJPL(q0), q;
        Eigen::Matrix<double,4,1> right;
        right[0] = q1.x(); right[1] = q1.y();
        right[2] = q1.z(); right[3] = q1.w();
        Eigen::Matrix<double,4,1> result = left*right;
        Eigen::Quaterniond q_result;
        q_result.x() = result[0]; q_result.y() = result[1];
        q_result.z() = result[2]; q_result.w() = result[3];
        return q_result;
    }

    ///
    /// @brief right times
    ///
    inline Eigen::Quaterniond RightTimes(
        const Eigen::Quaterniond & q0,
        const Eigen::Quaterniond & q1)
    {
        Eigen::Matrix<double, 4, 4> right;
        Eigen::Matrix<double,4,1> q;
        q[0] = q1.x(); q[1] = q1.y();
        q[2] = q1.z(); q[3] = q1.w();
        right << RightJPL(q1), q;
        Eigen::Matrix<double,4,1> left;
        left[0] = q0.x(); left[1] = q0.y();
        left[2] = q0.z(); left[3] = q0.w();
        Eigen::Matrix<double,4,1> result = right*left;
        Eigen::Quaterniond q_result;
        q_result.x() = result[0]; q_result.y() = result[1];
        q_result.z() = result[2]; q_result.w() = result[3];
        return q_result;
    }

    ///
    /// @brief get the JPL quaternion's rotation matrix
    /// @param _q JPL quaternion
    ///
    inline Eigen::Matrix<double, 3, 3> toRotationMatJPL(
        const Eigen::Quaterniond & _q) {

      const double tx  = double(2)*_q.x();
      const double ty  = double(2)*_q.y();
      const double tz  = double(2)*_q.z();
      const double twx = tx*_q.w();
      const double twy = ty*_q.w();
      const double twz = tz*_q.w();
      const double txx = tx*_q.x();
      const double txy = ty*_q.x();
      const double txz = tz*_q.x();
      const double tyy = ty*_q.y();
      const double tyz = tz*_q.y();
      const double tzz = tz*_q.z();

      return (Eigen::Matrix<double, 3, 3>() <<
          double(1)-(tyy+tzz), txy+twz, txz-twy, txy-twz,
          double(1)-(txx+tzz), tyz+twx, txz+twy, tyz-twx,
          double(1)-(txx+tyy)).finished();
    }    

    ///
    /// @brief Right Jacobian of JPL quaternion
    /// @param phi angle vector
    ///
    inline Eigen::Matrix<double, 3, 3> getRightJacobian(
        const Eigen::Vector3d & phi) {
        return (Eigen::Matrix<double, 3, 3>() <<
              Eigen::Matrix<double, 3, 3>::Identity()
              - ((1-cos(phi.norm()))/std::pow(phi.norm(),2))*skewJPL(phi)
              + ((phi.norm() - sin(phi.norm()))/std::pow(phi.norm(),3))
                *skewJPL(phi)*skewJPL(phi)).finished();
    }

    inline Eigen::Quaterniond fromRotationMatJPL(
            Eigen::Matrix3d &mat)
    {
        Eigen::Quaterniond q;
        double t = mat.trace();
        if (t > 0.0)
        {
          t = sqrt(t + 1.0);
          q.w() = 0.5*t;
          t = 0.5/t;
          q.x() = (mat.coeff(1,2) - mat.coeff(2,1)) * t;
          q.y() = (mat.coeff(2,0) - mat.coeff(0,2)) * t;
          q.z() = (mat.coeff(0,1) - mat.coeff(1,0)) * t;
        }
        else
        {
          int i = 0;
          if (mat.coeff(1,1) > mat.coeff(0,0))
            i = 1;
          if (mat.coeff(2,2) > mat.coeff(i,i))
            i = 2;
          int j = (i+1)%3;
          int k = (j+1)%3;

          t = sqrt(mat.coeff(i,i)-mat.coeff(j,j)-mat.coeff(k,k) + 1.0);
          q.coeffs().coeffRef(i) = 0.5 * t;
          t = 0.5/t;
          q.w() = (mat.coeff(j,k)-mat.coeff(k,j))*t;
          q.coeffs().coeffRef(j) = (mat.coeff(j,i)+mat.coeff(i,j))*t;
          q.coeffs().coeffRef(k) = (mat.coeff(k,i)+mat.coeff(i,k))*t;
        }
        return q;
    }
}


#endif // JPL.h
