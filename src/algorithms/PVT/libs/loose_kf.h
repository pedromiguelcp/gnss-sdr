/*!
 * \file loose_kf.h
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

#ifndef GNSS_SDR_LOOSE_KF_H_
#define GNSS_SDR_LOOSE_KF_H_

#include "inertial_navigator.h"

/** \addtogroup PVT
 * \{ */
/** \addtogroup PVT_libs
 * \{ */


class Loose_Kf
{
public:
    // CONSTRUCTOR
    Loose_Kf();
    // DESTRUCTOR
    ~Loose_Kf();

    arma::vec::fixed<15> _Xpre;  // Predicted State Vector
    arma::vec::fixed<15> _Xupd;  // Measurement Updated State Vector
    arma::mat _F;                // State Transition Matrix
    arma::mat _G;                // Process Noise Coefficient Matrix
    arma::mat _Q;                // Process Noise Covariance
    arma::mat S;
    arma::mat _Ppre;             // State Covariance (Predicted)
    arma::mat _Pupd;             // State Covariance (Updated)
    arma::vec::fixed<6> _Zobs;   // IMU Observation Vector
    arma::mat _Robs;             // IMU Observation Variance Matrix
    double dT;

    arma::vec3 nis_pos_i;
    arma::vec3 nis_vel_i;
    double nis_pos;
    double nis_vel;
    arma::vec3 y_pos;
    arma::vec3 y_vel;
    bool pos_ok;
    bool vel_ok;


    // Functions
    void Transition(const double rx_dT, const Inertial_Navigator& imu);
    void ProcessNoiseCoeff(const Inertial_Navigator& imu);
    void SetObs(const Inertial_Navigator& imu, const arma::vec3& GNSS_Pxyz, const arma::vec3& GNSS_Vxyz, const float* qr);
    void Filter(Inertial_Navigator& imu);
    void clearKF();
    void PredictOnly();
};

/** \} */
/** \} */
#endif  // GNSS_SDR_LOOSE_KF_H