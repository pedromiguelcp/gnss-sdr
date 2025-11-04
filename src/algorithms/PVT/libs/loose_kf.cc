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
#include "dcm.h"
#include "rtklib_rtkcmn.h"
#include "inertial_navigator.h"


// Develop transition matrix
void Loose_Kf::Transition(double rx_dT, const Inertial_Navigator& imu)
{
    dT = rx_dT;
    _F = arma::zeros(15, 15);
    const arma::vec3 om_ei_e = {0.0, 0.0, GNSS_OMEGA_EARTH_DOT};

    // eq 14.50
    // Derivatives wrt position model
    // position
    _F.submat(0, 0, 2, 2) = arma::eye(3, 3);
    // velocity
    _F.submat(0, 3, 2, 5) = arma::eye(3, 3) * dT;

    // Derivatives wrt velocity model
    // position (maps a position error into a velocity error - gravity changes with position.)
    _F.submat(3, 0, 5, 2) = dT * imu.Ne;
    // velocity
    _F.submat(3, 3, 5, 5) = arma::eye(3, 3) - 2 * dT * SkewMat(om_ei_e);
    // attitude
    _F.submat(3, 6, 5, 8) = -dT * imu.Fe;
    // specific force
    _F.submat(3, 12, 5, 14) = dT * imu.Ceb;

    // Derivatives wrt Attitude model
    // attitude
    _F.submat(6, 6, 8, 8) = arma::eye(3, 3) - dT * SkewMat(om_ei_e);
    // angular rate
    _F.submat(6, 9, 8, 11) = dT * imu.Ceb;

    // Derivatives wrt Angular Rate Model
    // angular rate
    _F.submat(9, 9, 11, 11) = arma::eye(3, 3);

    // Derivatives wrt Specific Force model
    // specific force
    _F.submat(12, 12, 14, 14) = arma::eye(3, 3);
}

// Develop process noise coefficient matrix
void Loose_Kf::ProcessNoiseCoeff(const Inertial_Navigator& imu)
{
    _G.zeros(15, 12);

    // attitude is driven by gyro noise
    _G.submat(6, 0, 8, 2) = imu.Ceb;

    // velocity is driven by accel noise
    _G.submat(3, 3, 5, 5) = imu.Ceb;

    // gyro-bias random walk
    _G.submat(9, 6, 11, 8) = arma::eye(3, 3);

    // accel-bias random walk
    _G.submat(12, 9, 14, 11) = arma::eye(3, 3);

    // continuous-time noise PSDs
    arma::vec q = {
        4.3633e-4, 4.3633e-4, 4.3633e-4,  // gyro white  (rad/s)/sqrt(Hz)
        4.9033e-3, 4.9033e-3, 4.9033e-3,  // accel white m/s^2/sqrt(Hz)
        1.6160e-7, 1.6160e-7, 1.6160e-7,  // gyro-bias RW  (rad/s)/sqrt(s)
        1.6344e-4, 1.6344e-4, 1.6344e-4   // accel-bias RW m/s^2/sqrt(s)
    };
    arma::mat Qw = arma::diagmat(arma::square(q));

    // Discretize (first order)
    _Q = _G * Qw * _G.t() * dT;
}


// Sets observation vector of IMU
void Loose_Kf::SetObs(const Inertial_Navigator& imu, const arma::vec3& GNSS_Pxyz, const arma::vec3& GNSS_Vxyz, const float* qr)
{
    // measurement = observed - computed
    _Zobs = arma::zeros(6);
    _Zobs.subvec(0, 2)  = GNSS_Pxyz - imu.pos_ant_ecef; // imu.obs_pva.pos_ecef - imu.pos_ant_ecef;
    _Zobs.subvec(3, 5)  = GNSS_Vxyz - imu.vel_ant_ecef; // imu.obs_pva.vel_ecef - imu.vel_ant_ecef;

    // measurement covariance
    _Robs.zeros(6, 6);
    // position block (ECEF)
    arma::mat Rpos(3, 3, arma::fill::zeros);
    Rpos(0, 0) = qr[0];
    Rpos(1, 1) = qr[1];
    Rpos(2, 2) = qr[2];
    Rpos(0, 1) = Rpos(1, 0) = qr[3];
    Rpos(1, 2) = Rpos(2, 1) = qr[4];
    Rpos(2, 0) = Rpos(0, 2) = qr[5];
    _Robs.submat(0, 0, 2, 2) = Rpos;
    // velocity block (ECEF)
    _Robs.submat(3, 3, 5, 5) = 0.5 * Rpos;

    // weight from innovation
    _Robs.diag() += 0.01 * arma::square(_Zobs);
}


// Kalman Filter Algorithm
void Loose_Kf::Filter(Inertial_Navigator& imu)
{
    // H selects position and velocity
    arma::mat H = arma::zeros(6, 15);
    H.submat(0, 0, 2, 2) = -arma::eye(3, 3);
    H.submat(3, 3, 5, 5) = -arma::eye(3, 3);
    // position/velocity errors at antenna depend on level arm | eq 14.111
    arma::vec3 h_r = imu.Ceb * imu.Lxyz;
    arma::vec3 om_b = imu.obs.Gyr - imu.GYRbias_b;
    const arma::vec3 om_ei_e = {0.0, 0.0, GNSS_OMEGA_EARTH_DOT};
    arma::vec3 h_v = imu.Ceb * arma::cross(om_b, imu.Lxyz) - arma::cross(om_ei_e, h_r);
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
    arma::mat PHt = _Ppre * H.t();
    arma::mat S = H * PHt + _Robs;
    // enforce symmetry + tiny diagonal load
    S = 0.5 * (S + S.t());
    S.diag() += 1e-12;
    arma::mat K = arma::solve(S, PHt.t()).t();

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
    _G = arma::zeros(u, 12);

    // State variance
    _Pupd(0, 0) = pow(5, 2);
    _Pupd(1, 1) = pow(5, 2);
    _Pupd(2, 2) = pow(7, 2);
    _Pupd(3, 3) = pow(2, 2);
    _Pupd(4, 4) = pow(2, 2);
    _Pupd(5, 5) = pow(3, 2);
    _Pupd(6, 6) = pow(0.1, 2);
    _Pupd(7, 7) = pow(0.1, 2);
    _Pupd(8, 8) = pow(0.1, 2);
    _Pupd(9, 9) = pow(0.001, 2);
    _Pupd(10, 10) = pow(0.001, 2);
    _Pupd(11, 11) = pow(0.001, 2);
    _Pupd(12, 12) = pow(0.001, 2);
    _Pupd(13, 13) = pow(0.001, 2);
    _Pupd(14, 14) = pow(0.001, 2);
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