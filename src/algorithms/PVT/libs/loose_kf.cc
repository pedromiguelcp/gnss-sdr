/*!
 * \file loose_kf.cc
 * \brief Class that implements a Loosely Coupled integration
 * \author Pedro Pereira, 2025. pereirapedrocp@gmail.com
 *
 * -----------------------------------------------------------------------------
 *
 * GNSS-SDR is a Global Navigation Satellite System software-defined receiver.
 * This file is part of GNSS-SDR.
 *
 * Copyright (C) 2010-2022  (see AUTHORS file for a list of contributors)
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * -----------------------------------------------------------------------------
 */

#include "loose_kf.h"
#include "inertial_navigator.h"

const double Om = 7.2921155e-5;

// Develop transition matrix
void Loose_Kf::Transition(double dT, arma::mat Ne, arma::mat Fe, arma::mat Ceb)
{
    _F = arma::zeros(15, 15);
    arma::vec3 Oeie = arma::zeros(3);
    Oeie(2) = Om;

    // eq 14.50
    // --- Derivatives wrt position model ---
    // position
    _F.submat(0, 0, 2, 2) = arma::eye(3, 3);
    // velocity
    _F.submat(0, 3, 2, 5) = arma::eye(3, 3) * dT;

    // --- Derivatives wrt velocity model ---
    // position (maps a position error into a velocity error - gravity changes with position.)
    _F.submat(3, 0, 5, 2) = Ne * dT;
    // velocity
    _F.submat(3, 3, 5, 5) = arma::eye(3, 3) - 2 * dT * SkewMat(Oeie);
    // attitude
    _F.submat(3, 6, 5, 8) = -dT * Fe;
    // specific force
    _F.submat(3, 12, 5, 14) = dT * Ceb;

    // --- Derivatives wrt Attitude model ---
    // attitude
    _F.submat(6, 6, 8, 8) = arma::eye(3, 3) - dT * SkewMat(Oeie);
    // angular rate
    _F.submat(6, 9, 8, 11) = dT * Ceb;

    // --- Derivatives wrt Angular Rate Model ---
    // angular rate
    _F.submat(9, 9, 11, 11) = arma::eye(3, 3);

    // --- Derivatives wrt Specific Force model ---
    // specific force
    _F.submat(12, 12, 14, 14) = arma::eye(3, 3);
}

// Develop process noise coefficient matrix
void Loose_Kf::ProcessNoiseCoeff(double dT, arma::mat Ceb)
{
    _G.zeros(15, 12);

    // attitude is driven by gyro noise
    _G.submat(6, 0, 8, 2) = Ceb;

    // velocity is driven by accel noise
    _G.submat(3, 3, 5, 5) = Ceb;

    // gyro-bias random walk
    _G.submat(9, 6, 11, 8) = arma::eye(3, 3);

    // accel-bias random walk
    _G.submat(12, 9, 14, 11) = arma::eye(3, 3);

    // continuous-time noise PSDs
    arma::vec q = {
        0.002, 0.002, 0.002,  // gyro white
        0.03, 0.03, 0.03,     // accel white
        5e-6, 5e-6, 5e-6,     // gyro bias RW
        2e-4, 2e-4, 2e-4      // accel bias RW
    };
    arma::mat Qw = arma::diagmat(arma::square(q));

    // Discretize (first order)
    _Q = _G * Qw * _G.t() * dT;
}


// Sets observation vector of IMU
void Loose_Kf::SetObs(Inertial_Navigator& imu)
{
    _Zobs = arma::zeros(6);
    _Robs = arma::zeros(6, 6);

    for (int i = 0; i < 3; ++i)
        {
            _Zobs(i) = imu.obs_pva.pos_ecef[i] - imu.pos_ant_ecef(i);
            _Zobs(i + 3) = imu.obs_pva.vel_ecef[i] - imu.vel_ant_ecef(i);
        }

    // measurement covariance (small - high trust in iTrace)
    _Robs.submat(0, 0, 2, 2) = arma::eye(3, 3) * 0.05 * 0.05;  // 5 cm (ECEF)
    _Robs.submat(3, 3, 5, 5) = arma::eye(3, 3) * 0.05 * 0.05;  // 5 cm/s
}


// Kalman Filter Algorithm
void Loose_Kf::Filter(Inertial_Navigator& imu)
{
    // H selects position and velocity
    arma::mat H = arma::zeros(6, 15);
    H.submat(0, 0, 2, 2) = -arma::eye(3, 3);
    H.submat(3, 3, 5, 5) = -arma::eye(3, 3);
    // position/velocity errors at antenna depen on level arm | eq 14.111
    arma::vec3 h_r = imu.Ceb * imu.Lxyz;
    arma::vec3 om_b = imu.obs.Gyr - imu.GYRbias_b;
    arma::vec3 h_v = imu.Ceb * arma::cross(om_b, imu.Lxyz);
    arma::mat H_gyr = imu.Ceb * SkewMat(imu.Lxyz);
    H.submat(0, 6, 2, 8) = SkewMat(h_r);
    H.submat(3, 6, 5, 8) = SkewMat(h_v);
    H.submat(3, 9, 5, 11) = H_gyr;

    // predict
    _Xupd = arma::zeros(15);  // closed-loop
    _Xpre = _F * _Xupd;
    _Ppre = _F * _Pupd * _F.t() + _Q;

    // innovation (GNSS − INS)
    arma::vec y = _Zobs - H * _Xpre;

    // Kalman gain
    arma::mat S = H * _Ppre * H.t() + _Robs;
    arma::mat Sinv = arma::inv_sympd(S);
    arma::mat K = _Ppre * H.t() * Sinv;

    // state update
    _Xupd = _Xpre + K * y;

    // covariance update (Joseph form)
    arma::mat I15 = arma::eye(15, 15);
    arma::mat IKH = I15 - K * H;
    _Pupd = IKH * _Ppre * IKH.t() + K * _Robs * K.t();
    _Pupd = 0.5 * (_Pupd + _Pupd.t());  // enforce symmetry


    // Update States
    sol.posXYZ = imu.pos_ant_ecef - _Xupd.subvec(0, 2);
    sol.velXYZ = imu.vel_ant_ecef - _Xupd.subvec(3, 5);
    sol.attXYZ = imu.att_rpy - _Xupd.subvec(6, 8);
    sol.dw = imu.GYRbias_b + _Xupd.subvec(9, 11);
    sol.df = imu.ACCbias_b + _Xupd.subvec(12, 14);

    imu.pos_ant_ecef -= _Xupd.subvec(0, 2);
    imu.vel_ant_ecef -= _Xupd.subvec(3, 5);
    // imu.att_rpy -= _Xupd.subvec(6, 8);
    arma::vec3 dpsi = _Xupd.subvec(6, 8);
    imu.Ceb = (arma::eye(3, 3) - SkewMat(dpsi)) * imu.Ceb;
    imu.GYRbias_b += _Xupd.subvec(9, 11);
    imu.ACCbias_b += _Xupd.subvec(12, 14);
}

// CONSTRUCTOR AND DESTRUCTOR DEFINITIONS
Loose_Kf::Loose_Kf()
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

    // State variance
    // position (m^2), velocity (m/s)^2
    double qr = std::pow(2, 2);
    double qv = std::pow(0.2, 2);
    _Pupd(0, 0) = qr;
    _Pupd(1, 1) = qr;
    _Pupd(2, 2) = qr;
    _Pupd(3, 3) = qv;
    _Pupd(4, 4) = qv;
    _Pupd(5, 5) = qv;

    // attitude (rad^2) – allow ~1° initial misalignment
    double qatt = std::pow(0.017, 2);  // ≈ 1 deg
    _Pupd(6, 6) = qatt;
    _Pupd(7, 7) = qatt;
    _Pupd(8, 8) = qatt;

    // gyro bias (rad/s)^2
    double qbg = std::pow(0.02, 2);  // ≈ 72 deg/h
    _Pupd(9, 9) = qbg;
    _Pupd(10, 10) = qbg;
    _Pupd(11, 11) = qbg;

    // accel bias (m/s^2)^2 – generous (20 mg)
    double qba = std::pow(0.2, 2);
    _Pupd(12, 12) = qba;
    _Pupd(13, 13) = qba;
    _Pupd(14, 14) = qba;
}

Loose_Kf::~Loose_Kf() {}

void Loose_Kf::clearKF()
{
    _Xpre.resize(0);
    _Xupd.resize(0);
    _F.resize(0, 0);
    _G.resize(0, 0);
    _Q.resize(0, 0);
    _Ppre.resize(0, 0);
    _Pupd.resize(0, 0);
}

// A function to build skew symmetric matrix
arma::mat Loose_Kf::SkewMat(const arma::vec3& Vec)
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