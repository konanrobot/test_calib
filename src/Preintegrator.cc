#include "Preintegrator.h"

namespace ORB_SLAM2 {
void Preintegrator::getPreintResult(const vector<pair<double, double>> &gyro,
                                    const Eigen::Quaterniond &q_k,
                                    const Eigen::Vector3d &bias_bar,
                                    Eigen::Quaterniond &q_est,
                                    Eigen::Matrix3d &info,
                                    Eigen::Matrix3d &bias_jac)
{
    double time_tau = gyro[0].first;
    double omega_m_tau = gyro[0].second;
    Eigen::Vector4d q_tau;
    q_tau[0] = q_k.x(); q_tau[1] = q_k.y(); q_tau[2] = q_k.z(); q_tau[3] = q_k.w();
    Eigen::Matrix3d P_tau; P_tau.setZero();
    Eigen::Matrix3d Jq_tau; Jq_tau.setZero();
    Eigen::Vector4d q_tau_1;
    Eigen::Matrix3d P_tau_1;
    Eigen::Matrix3d Jq_tau_1;
    for(int i = 1; i < gyro.size(); i++)
    {
        double time_tau_1 = gyro[i].first;
        double omega_m_tau_1 = gyro[i].second;
        Eigen::Vector3d omega_m; omega_m.setZero();
        omega_m[2] = omega_m_tau;
        Eigen::Vector3d omega = omega_m - bias_bar;
        double theta = omega.norm();
        double delta_t = time_tau_1 - time_tau;
        Eigen::Vector3d omega_m_1; omega_m_1.setZero();
        omega_m_1[2] = omega_m_tau_1;
        Eigen::Vector3d omega_1 = omega_m_1 - bias_bar;
        if(theta<0.00001)
        {
            Eigen::Matrix4d M;
            M = Eigen::Matrix4d::Identity()
                    + delta_t/2*g2o::OmegaJPL(omega);
            q_tau_1 = M*q_tau;
            Eigen::Matrix3d Phi_tau = Eigen::Matrix3d::Identity()
                    + delta_t*g2o::skewJPL(-omega)
                    + delta_t*delta_t*g2o::skewJPL(-omega)*g2o::skewJPL(-omega)/2;
            Eigen::Matrix3d Qd = delta_t*0.05*0.05*Eigen::Matrix3d::Identity();
            P_tau_1 = Phi_tau*P_tau*Phi_tau.transpose()+Qd;
            Eigen::Matrix3d hat_R_tau_tau_1 = Eigen::Matrix3d::Identity()
                    + g2o::skewJPL(omega) + g2o::skewJPL(omega)*g2o::skewJPL(omega);
            Jq_tau_1 = hat_R_tau_tau_1*Jq_tau
                    +g2o::getRightJacobian(omega_1)*delta_t;
            time_tau = time_tau_1;
            omega_m_tau = omega_m_tau_1;
            q_tau = q_tau_1;
            P_tau = P_tau_1;
            Jq_tau = Jq_tau_1;
        }
        else
        {
            Eigen::Matrix4d M;
            M = cos(theta*delta_t/2)*Eigen::Matrix4d::Identity()
                    + sin(theta*delta_t/2)*g2o::OmegaJPL(omega)/theta;
            q_tau_1 = M*q_tau;
            double theta1 = (-omega).norm();
            Eigen::Matrix3d Phi_tau = Eigen::Matrix3d::Identity()
                    + sin(theta1*delta_t)*g2o::skewJPL(-omega)/theta1
                    + (1-cos(theta1*delta_t))*g2o::skewJPL(-omega)*g2o::skewJPL(-omega)/(theta1*theta1);
            Eigen::Matrix3d Qd = delta_t*0.05*0.05*Eigen::Matrix3d::Identity();
            P_tau_1 = Phi_tau*P_tau*Phi_tau.transpose()+Qd;
            Eigen::Matrix3d hat_R_tau_tau_1 = Eigen::Matrix3d::Identity()
                    + sin(theta)*g2o::skewJPL(omega)/theta
                    + (1-cos(theta))*g2o::skewJPL(omega)*g2o::skewJPL(omega)/(theta*theta);
            Jq_tau_1 = hat_R_tau_tau_1*Jq_tau
                    +g2o::getRightJacobian(omega_1)*delta_t;
            time_tau = time_tau_1;
            omega_m_tau = omega_m_tau_1;
            q_tau = q_tau_1;
            P_tau = P_tau_1;
            Jq_tau = Jq_tau_1;
        }
    }
    q_est.x() = q_tau_1[0]; q_est.y() = q_tau_1[1]; q_est.z() = q_tau_1[2];
    q_est.w() = q_tau_1[3];
    info = P_tau_1.inverse();
    bias_jac = Jq_tau_1;
}
}
