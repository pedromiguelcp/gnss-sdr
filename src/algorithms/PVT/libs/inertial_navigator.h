/*!
 * \file inertial_navigator.h
 * \brief Class that implements an inertial navigation mechanizator
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

#ifndef GNSS_SDR_INERTIAL_NAVIGATOR_H_
#define GNSS_SDR_INERTIAL_NAVIGATOR_H_

#include <armadillo>
#include <fstream>

/** \addtogroup PVT
 * \{ */
/** \addtogroup PVT_libs
 * \{ */


class Inertial_Navigator
{
public:
    Inertial_Navigator();
    ~Inertial_Navigator() = default;

    // -------- Data structures --------
    struct IMUEpochInfo
    {
        double imuTime{};
        arma::vec3 Acc;  // ax, ay, az (m/s^2)
        arma::vec3 Gyr;  // gx, gy, gz (rad/s)
    };

    struct IMUPVAEpochInfo
    {
        double imuPVATime{};
        arma::vec3 pos_llh;  // lat, long, h (deg, deg, m)
        arma::vec3 pos_ecef;
        arma::vec3 vel_ned;  // v_n, v_e, v_d (m/s)
        arma::vec3 vel_ecef;
        arma::vec3 att_rpy;  // roll, pitch, yaw (rad)
    };

    // -------- File reading --------
    // Reads one epoch from a Septentrio IMU log. If readHeader=true, skips header line.
    void readIMU(std::ifstream& fin_raw_imu_file, bool readHeader);

    // -------- ECEF mechanization --------
    void initializeMechanizer(std::ifstream& fin_raw_imu, double EndTime, const arma::vec3& iniPOS_ecef,
        const arma::vec3& iniVEL_ecef);

    // One propagation step (dT seconds)
    void stepMechanizer(std::ifstream& fin_raw_imu);

    // Predict antenna ECEF position/velocity from IMU-origin state at this epoch
    void predictAntennaECEF();
    void predictIMUECEF();

    // Align yaw from a GNSS PVT (position+velocity in ECEF).
    // yaw_boresight_rad lets you account for any IMU-to-vehicle yaw offset.
    // If snap_state_to_gnss=true, we also overwrite pos_ecef/vel_ecef with the GNSS values.
    void alignYawFromGnss(const arma::vec3& pos_ecef_gnss,
        const arma::vec3& vel_ecef_gnss,
        double yaw_boresight_rad = 0.0);

    bool yawAligned() const { return yaw_aligned_; }


    // -------- Public state --------
    IMUEpochInfo obs{};         // latest observation read
    IMUPVAEpochInfo obs_pva{};  // latest observation read

    // Init results
    arma::vec3 GYRbias_b;
    arma::vec3 ACCbias_b;

    // Mech state
    arma::vec3 pos_ecef;      // x,y,z (m)
    arma::vec3 pos_llh;       // lat, long, height (rad,rad,m)
    arma::vec3 vel_ecef;      // vx,vy,vz (m/s)
    arma::vec3 pos_ant_ecef;  // x,y,z (m)
    arma::vec3 vel_ant_ecef;  // vx,vy,vz (m/s)
    arma::vec3 att_rpy;       // roll,pitch,yaw (rad)
    arma::mat Ceb;            // body->ECEF
    arma::mat Ne;             // gravity gradient
    arma::mat Fe;
    arma::vec3 Lxyz;  // lever arm (m)

private:
    bool yaw_aligned_ = false;

    double normalise(const double value, const double start, const double end);
    arma::mat GravGrad(const arma::vec3& Vec);

    // GNSS velocity EMA in NED for horizontal-only correction
    arma::vec3 v_ned_gnss_filt;
    // Smoothing gains per GNSS epoch
    double ema_alpha = 0.2;
};

/** \} */
/** \} */
#endif  // GNSS_SDR_INERTIAL_NAVIGATOR_H