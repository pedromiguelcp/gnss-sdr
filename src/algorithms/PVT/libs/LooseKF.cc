/*
 * LooseKF.h
 * 21 state Kalman Filter Data Holder for Loosely Coupled Integration
 *  Created on: Sept 12, 2018
 *      Author: Aaron Boda
 */

#include "LooseKF.h"
#include "inertial_navigator.h"

const double Om = 7.2921155e-5;

// Develop transition matrix
void LooseKF::Transition(double dT, arma::mat Ne, arma::mat Fe, arma::mat Ceb)
{
    _F = arma::zeros(15, 15);
    arma::mat dfdx = arma::zeros(3, 3);
    arma::vec3 Oeie = arma::zeros(3);
    Oeie(2) = Om;

    // --- Derivatives wrt position model ---
    // 1) for position
    dfdx = arma::eye(3, 3);
    _F.submat(0, 0, 2, 2) = dfdx;
    // 2) for velocity
    dfdx = arma::eye(3, 3);
    dfdx = dfdx * dT;
    _F.submat(0, 3, 2, 5) = dfdx;

    // --- Derivatives wrt velocity model ---
    // 1) for position
    _F.submat(3, 0, 5, 2) = Ne;
    // 2) for velocity
    dfdx = arma::eye(3, 3) - 2 * dT * SkewMat(Oeie);
    _F.submat(3, 3, 5, 5) = dfdx;
    // 3) for attitude
    dfdx = -dT * Fe;
    _F.submat(3, 6, 5, 8) = dfdx;
    // 5) for specific force
    dfdx = dT * Ceb;
    _F.submat(3, 12, 5, 14) = dfdx;

    // --- Derivatives wrt Attitude model ---
    // 3) for attitude
    dfdx = arma::eye(3, 3) - dT * SkewMat(Oeie);
    _F.submat(6, 6, 8, 8) = dfdx;
    // 4) for angular rate
    dfdx = dT * Ceb;
    _F.submat(6, 9, 8, 11) = dfdx;

    // --- Derivatives wrt Angular Rate Model ---
    // 4) for angular rate
    arma::mat Dg(3, 3);
    Dg.fill(0.0007);
    dfdx = arma::eye(3, 3) + dT * Dg;
    _F.submat(9, 9, 11, 11) = dfdx;

    // --- Derivatives wrt Specific Force model ---
    // 5) for specific force
    arma::mat Da(3, 3);
    Da.fill(0.003);
    dfdx = arma::eye(3, 3) + dT * Da;
    _F.submat(12, 12, 14, 14) = dfdx;
}

// Develop process noise coefficient matrix
void LooseKF::ProcessNoiseCoeff(double dT, arma::mat Ceb)
{
    _G.zeros(15, 12);  // [gyro_noise(3), acc_noise(3), gyro_bias_rw(3), acc_bias_rw(3)]

    // attitude is driven by gyro noise  (δθ̇ ≈ Ceb*ng)
    _G.submat(6, 0, 8, 2) = Ceb;

    // velocity is driven by accel noise (δv̇ ≈ Ceb*na)
    _G.submat(3, 3, 5, 5) = Ceb;

    // gyro-bias random walk: ḃg = n_bg
    _G.submat(9, 6, 11, 8) = arma::eye(3, 3);

    // accel-bias random walk: ḃa = n_ba
    _G.submat(12, 9, 14, 11) = arma::eye(3, 3);

    // continuous-time noise PSDs (tune to your IMU; examples)
    // ng [rad/s/√Hz], na [m/s²/√Hz], bg_rw [rad/s√Hz], ba_rw [m/s²√Hz]
    arma::vec q = {
        0.002, 0.002, 0.002,  // gyro white
        0.03, 0.03, 0.03,     // accel white
        5e-6, 5e-6, 5e-6,     // gyro bias RW
        2e-4, 2e-4, 2e-4      // accel bias RW
    };
    arma::mat Qw = arma::diagmat(arma::square(q));

    // Discretize (first order): Q ≈ G Qw Gᵀ Δt
    _Q = _G * Qw * _G.t() * dT;
}


// Sets observation vector of IMU
void LooseKF::SetObs(Inertial_Navigator& imu)
{
    _Zobs = arma::zeros(6);
    _Robs = arma::zeros(6, 6);

    // innovation z = INS - truth  (pos[0..2], vel[3..5])
    for (int i = 0; i < 3; ++i)
        {
            _Zobs(i) = imu.pos_ant_ecef(i) - imu.obs_pva.pos_ecef[i];
            _Zobs(i + 3) = imu.vel_ant_ecef(i) - imu.obs_pva.vel_ecef[i];
        }

    // measurement covariance (tune; small if you trust iTrace)
    _Robs.submat(0, 0, 2, 2) = arma::eye(3, 3) * 0.05 * 0.05;  // 5 cm (ECEF)
    _Robs.submat(3, 3, 5, 5) = arma::eye(3, 3) * 0.05 * 0.05;  // 5 cm/s
}


// Kalman Filter Algorithm
void LooseKF::Filter(Inertial_Navigator& imu)
{
    arma::vec::fixed<6> I = arma::zeros(6);
    arma::mat K = arma::zeros(15, 6);
    arma::mat H = arma::eye(6, 15);

    // State Covariance Prediction
    _Xupd = arma::zeros(15);
    _Ppre = _F * _Pupd * _F.t() + _Q;
    // State Prediction
    _Xpre = _F * _Xupd;
    // Compute Innovation
    I = _Zobs;
    // Gain
    K = _Ppre * H.t() * arma::pinv(H * _Ppre * H.t() + _Robs);
    // State Estimate (Updated)
    _Xupd = _Xpre + (K * I);
    // State Covaraince (Updated)
    _Pupd = (arma::eye(15, 15) - K * H) * _Ppre;

    // Update States
    sol.posXYZ = imu.pos_ant_ecef - _Xupd.subvec(0, 2);
    sol.velXYZ = imu.vel_ant_ecef - _Xupd.subvec(3, 5);
    sol.attXYZ = imu.att_rpy - _Xupd.subvec(6, 8);
    sol.dw = imu.GYRbias_b + _Xupd.subvec(9, 11);
    sol.df = imu.ACCbias_b + _Xupd.subvec(12, 14);

    imu.pos_ant_ecef -= _Xupd.subvec(0, 2);
    imu.vel_ant_ecef -= _Xupd.subvec(3, 5);
    imu.att_rpy -= _Xupd.subvec(6, 8);
    imu.GYRbias_b += _Xupd.subvec(9, 11);
    imu.ACCbias_b += _Xupd.subvec(12, 14);
}

// CONSTRUCTOR AND DESTRUCTOR DEFINITIONS
LooseKF::LooseKF()
{
    // Number of states
    size_t u = 15;
    // Initialize matrices
    _Xpre = arma::zeros(u);
    _Xupd = arma::zeros(u);
    _F = arma::zeros(u, u);
    _Q = arma::zeros(u, u);
    _Ppre = arma::zeros(u, u);
    _Pupd = arma::zeros(u, u);
    _G = arma::zeros(u, 6);
    _Qw = arma::zeros(6, 6);
    _scaleFactor = 1.0;

    // Process Noise Covariance
    _Qw(0, 0) = pow(0.013, 2);
    _Qw(1, 1) = pow(0.013, 2);
    _Qw(2, 2) = pow(0.013, 2);
    _Qw(3, 3) = pow(0.160, 2);
    _Qw(4, 4) = pow(0.160, 2);
    _Qw(5, 5) = pow(0.160, 2);

    // State variance
    _Pupd(0, 0) = pow(0.250, 2);
    _Pupd(1, 1) = pow(0.250, 2);
    _Pupd(2, 2) = pow(0.250, 2);
    _Pupd(3, 3) = pow(0.050, 2);
    _Pupd(4, 4) = pow(0.050, 2);
    _Pupd(5, 5) = pow(0.050, 2);
    _Pupd(6, 6) = pow(0.001, 2);
    _Pupd(7, 7) = pow(0.001, 2);
    _Pupd(8, 8) = pow(0.001, 2);
    _Pupd(9, 9) = pow(0.0007, 2);
    _Pupd(10, 10) = pow(0.0007, 2);
    _Pupd(11, 11) = pow(0.0007, 2);
    _Pupd(12, 12) = pow(0.004, 2);
    _Pupd(13, 13) = pow(0.003, 2);
    _Pupd(14, 14) = pow(0.003, 2);
}
LooseKF::~LooseKF() {}

void LooseKF::clearKF()
{
    _Xpre.resize(0);
    _Xupd.resize(0);
    _F.resize(0, 0);
    _G.resize(0, 0);
    _Q.resize(0, 0);
    _Qw.resize(0, 0);
    _Ppre.resize(0, 0);
    _Pupd.resize(0, 0);
}

// A function to build skew symmetric matrix
arma::mat LooseKF::SkewMat(const arma::vec3& Vec)
{
    arma::mat Skew = arma::zeros(3, 3);
    Skew(0, 1) = -Vec(2);
    Skew(0, 2) = Vec(1);
    Skew(1, 0) = Vec(2);
    Skew(1, 2) = -Vec(0);
    Skew(2, 0) = -Vec(1);
    Skew(2, 1) = Vec(0);
    return Skew;
}