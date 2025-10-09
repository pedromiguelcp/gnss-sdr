#pragma once
/*
 * DCM.h
 * Direction Cosine Matrices
 *  Created on: Aug 07, 2018
 *      Author: Aaron Boda
 */

#include <armadillo>

#ifndef DCM_H_
#define DCM_H_

arma::mat RotateX(double ang);
arma::mat RotateY(double ang);
arma::mat RotateZ(double ang);
arma::mat e2llfDCM(double lat, double lon);
arma::mat b2llfDCM(double roll, double pitch, double yaw);
arma::mat b2eDCM(double lat, double lon, double r, double p, double y);
arma::vec3 dcm2euler(arma::mat Cnb);

#endif /* DCM_H_ */
