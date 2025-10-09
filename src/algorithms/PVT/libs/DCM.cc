/*
 * DCM.cpp
 * Direction Cosine Matrices
 *  Created on: Aug 07, 2018
 *      Author: Aaron Boda
 */

#include "DCM.h"
#include <armadillo>
#include <algorithm>
#include <cmath>
#include <vector>


// Rotation Matrix Rx
arma::mat RotateX(double theta)
{
    arma::mat R(3, 3);
    R.zeros(3, 3);
    R(0, 0) = 1;
    R(1, 1) = cos(theta);
    R(1, 2) = sin(theta);
    R(2, 1) = -sin(theta);
    R(2, 2) = cos(theta);
    return R;
}

// Rotation Matrix Ry
arma::mat RotateY(double theta)
{
    arma::mat R(3, 3);
    R.zeros(3, 3);
    R(0, 0) = cos(theta);
    R(0, 2) = -sin(theta);
    R(1, 1) = 1;
    R(2, 0) = sin(theta);
    R(2, 2) = cos(theta);
    return R;
}

// Rotation Matrix Rz
arma::mat RotateZ(double theta)
{
    arma::mat R(3, 3);
    R.zeros(3, 3);
    R(0, 0) = cos(theta);
    R(0, 1) = sin(theta);
    R(1, 0) = -sin(theta);
    R(1, 1) = cos(theta);
    R(2, 2) = 1;
    return R;
}

// ECEF to LLF DCM
arma::mat e2llfDCM(double lat, double lon)
{
    arma::mat Cne = arma::zeros(3, 3);
    Cne(0, 0) = -cos(lon) * sin(lat);
    Cne(0, 1) = -sin(lat) * sin(lon);
    Cne(0, 2) = cos(lat);
    Cne(1, 0) = -sin(lon);
    Cne(1, 1) = cos(lon);
    Cne(1, 2) = 0;
    Cne(2, 0) = -cos(lat) * cos(lon);
    Cne(2, 1) = -cos(lat) * sin(lon);
    Cne(2, 2) = -sin(lat);
    return Cne;
}

// Body to LLF DCM
arma::mat b2llfDCM(double r, double p, double y)
{
    arma::mat Cbn = arma::zeros(3, 3);
    arma::mat Cnb = arma::zeros(3, 3);
    Cbn = (RotateX(r) * RotateY(p)) * RotateZ(y);
    Cnb = Cbn.t();
    return Cnb;
}

// Body to ECEF DCM
arma::mat b2eDCM(double lat, double lon, double r, double p, double y)
{
    arma::mat Cne = arma::zeros(3, 3);
    arma::mat Cnb = arma::zeros(3, 3);
    arma::mat Ceb = arma::zeros(3, 3);
    Cne = e2llfDCM(lat, lon);
    Cnb = b2llfDCM(r, p, y);
    Ceb = Cne.t() * Cnb;
    return Ceb;
}

// Find RPY, given enu2b1 DCM
arma::vec3 dcm2euler(arma::mat Cbn)
{
    const double roll = std::atan2(Cbn(1, 2), Cbn(2, 2));
    const double s = std::clamp(-Cbn(0, 2), -1.0, 1.0);
    const double pitch = std::asin(s);
    const double yaw = std::atan2(Cbn(0, 1), Cbn(0, 0));
    return {roll, pitch, yaw};
}