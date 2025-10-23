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

    struct KFupd
    {
        arma::vec3 posXYZ;
        arma::vec3 velXYZ;
        arma::vec3 attXYZ;
        arma::vec3 df;
        arma::vec3 dw;
    };

    // Attributes
    KFupd sol;
    arma::vec::fixed<15> _Xpre;  // Predicted State Vector
    arma::vec::fixed<15> _Xupd;  // Measurement Updated State Vector
    arma::mat _F;                // State Transition Matrix
    arma::mat _G;                // Process Noise Coefficient Matrix
    arma::mat _Q;                // Process Noise Covariance
    arma::mat _Ppre;             // State Covariance (Predicted)
    arma::mat _Pupd;             // State Covariance (Updated)
    arma::mat _Zobs;             // IMU Observation Vector
    arma::mat _Robs;             // IMU Observation Variance Matrix

    // Functions
    void Transition(double dT, arma::mat Ne, arma::mat Fe, arma::mat Ceb);
    void ProcessNoiseCoeff(double dT, arma::mat Ceb);
    void SetObs(Inertial_Navigator& IMU_NAV);
    void Filter(Inertial_Navigator& IMU_NAV);
    void clearKF();

private:
    // Functions

    arma::mat SkewMat(const arma::vec3& Vec);
};

/** \} */
/** \} */
#endif  // GNSS_SDR_LOOSE_KF_H