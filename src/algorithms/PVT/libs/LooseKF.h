/*
 * LooseKF.h
 * 21 state Kalman Filter for GNSS-IMU Loosely Coupled Integration
 *  Created on: Sept 12, 2018
 *      Author: Aaron Boda
 */

#ifndef LOOSEKF_H_
#define LOOSEKF_H_

#include "inertial_navigator.h"

class LooseKF
{
public:
    // CONSTRUCTOR
    LooseKF();
    // DESTRUCTOR
    ~LooseKF();

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
    arma::mat _Qw;               // Process Noise Covariance
    arma::mat _Q;                // Process Noise Covariance
    arma::mat _Ppre;             // State Covariance (Predicted)
    arma::mat _Pupd;             // State Covariance (Updated)
    arma::mat _Zobs;             // IMU Observation Vector
    arma::mat _Robs;             // IMU Observation Variance Matrix
    double _scaleFactor;         // Scale Factor of covariance matrix

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

#endif /* LOOSEKF_H_ */