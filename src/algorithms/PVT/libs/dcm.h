/*!
 * \file dcm.h
 * \brief Class that implements direction cosine matrices (DMC)
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

#ifndef GNSS_SDR_DCM_H_
#define GNSS_SDR_DCM_H_

#include <armadillo>

/** \addtogroup PVT
 * \{ */
/** \addtogroup PVT_libs
 * \{ */


arma::mat RotateX(double ang);
arma::mat RotateY(double ang);
arma::mat RotateZ(double ang);
arma::mat e2llfDCM(double lat, double lon);
arma::mat b2llfDCM(double roll, double pitch, double yaw);
arma::mat b2eDCM(double lat, double lon, double r, double p, double y);
arma::vec3 dcm2euler(arma::mat Cnb);


/** \} */
/** \} */
#endif  // GNSS_SDR_DCM_H