/*!
 * \file inertial_navigator.cc
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

#include "inertial_navigator.h"
#include "geofunctions.h"
#include "rtklib_rtkcmn.h"


Inertial_Navigator::Inertial_Navigator()
{
    pos_ecef = arma::zeros(3);
    pos_llh = arma::zeros(3);
    vel_ecef = arma::zeros(3);
    att_rpy = arma::zeros(3);
    GYRbias_b = arma::zeros(3);
    ACCbias_b = arma::zeros(3);
    Ceb = arma::zeros(3, 3);
    Fe = arma::zeros(3, 3);
    Ne = arma::zeros(3, 3);
    Lxyz = arma::zeros(3);
    v_ned_gnss_filt = arma::zeros(3);
}


void Inertial_Navigator::readIMU(std::ifstream& fin_raw_imu_file, bool readHeader)
{
    std::string line;
    if (readHeader)
        {
            getline(fin_raw_imu_file, line);
        }

    getline(fin_raw_imu_file, line);
    if (!line.empty())
        {
            std::istringstream iss(line);
            std::vector<std::string> words{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
            obs.imuTime = std::stod(words[17]);

            // raw data from file
            const double ax_s = std::stod(words[19]);
            const double ay_s = std::stod(words[20]);
            const double az_s = std::stod(words[21]);
            const double gx_s = std::stod(words[22]);
            const double gy_s = std::stod(words[23]);
            const double gz_s = std::stod(words[24]);

            obs.Acc(0) = ay_s;   // y-axis forward
            obs.Acc(1) = ax_s;   // x-axis right
            obs.Acc(2) = -az_s;  // +g when static (sensor Z up)
            obs.Gyr(0) = gy_s;
            obs.Gyr(1) = gx_s;
            obs.Gyr(2) = -gz_s;
        }
}


void Inertial_Navigator::readIMUPVA(std::ifstream& fin_pva_imu_file, bool readHeader)
{
    std::string line;
    if (readHeader)
        {
            getline(fin_pva_imu_file, line);
        }

    getline(fin_pva_imu_file, line);
    if (!line.empty())
        {
            std::istringstream iss(line);
            std::vector<std::string> words{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
            obs_pva.imuPVATime = std::stod(words[8]);

            // raw data from file
            const double pos_lat = std::stod(words[9]);
            const double pos_lon = std::stod(words[10]);
            const double pos_height = std::stod(words[11]);
            const double vel_north = std::stod(words[12]);
            const double vel_east = std::stod(words[13]);
            const double vel_up = std::stod(words[14]);
            const double roll = std::stod(words[15]);
            const double pitch = std::stod(words[16]);
            const double yaw = std::stod(words[17]);

            obs_pva.pos_llh(0) = pos_lat * M_PI / 180;
            obs_pva.pos_llh(1) = pos_lon * M_PI / 180;
            obs_pva.pos_llh(2) = pos_height;
            obs_pva.vel_ned(0) = vel_north;
            obs_pva.vel_ned(1) = vel_east;
            obs_pva.vel_ned(2) = -vel_up;
            obs_pva.att_rpy(0) = roll * M_PI / 180;
            obs_pva.att_rpy(1) = pitch * M_PI / 180;
            obs_pva.att_rpy(2) = -yaw * M_PI / 180;

            pos2ecef(obs_pva.pos_llh.memptr(), obs_pva.pos_ecef.memptr());
            arma::mat Cne = e2llfDCM(obs_pva.pos_llh(0), obs_pva.pos_llh(1));
            obs_pva.vel_ecef = Cne.t() * obs_pva.vel_ned;
        }
}


void Inertial_Navigator::correctVelRPY(std::ifstream& fin_pva_imu_file, double EndTime)
{
    while ((!fin_pva_imu_file.eof()) && (obs_pva.imuPVATime < EndTime))
        {
            readIMUPVA(fin_pva_imu_file, false);
        }
}


void Inertial_Navigator::initializeMechanizer(std::ifstream& fin_raw_imu, std::ifstream& fin_pva_imu, double EndTime, const arma::vec3& iniPOS_ecef,
    const arma::vec3& iniVEL_ecef)
{
    // Lever Arm
    // Lxyz = {0, 0, 0};
    Lxyz = {1.184, -0.004, -1.152};  // antenna w.r.t. IMU in body frame

    pos_ant_ecef = iniPOS_ecef;
    vel_ant_ecef = iniVEL_ecef;
    predictIMUECEF();

    // read first line (header)
    readIMU(fin_raw_imu, true);
    readIMUPVA(fin_pva_imu, true);

    double Ax_sum = 0, Ay_sum = 0, Az_sum = 0;
    double Gx_sum = 0, Gy_sum = 0, Gz_sum = 0;
    int count = 0;

    while ((!fin_raw_imu.eof()) && (obs.imuTime < EndTime))
        {
            readIMU(fin_raw_imu, false);
            if ((EndTime - obs.imuTime) < 600)
                {
                    Ax_sum += obs.Acc(0);
                    Ay_sum += obs.Acc(1);
                    Az_sum += obs.Acc(2);
                    Gx_sum += obs.Gyr(0);
                    Gy_sum += obs.Gyr(1);
                    Gz_sum += obs.Gyr(2);
                    count++;
                }
        }
    if (count == 0) return;

    // sensor means
    arma::vec3 acc_mean = {Ax_sum / count, Ay_sum / count, Az_sum / count};
    arma::vec3 gyr_mean = {Gx_sum / count, Gy_sum / count, Gz_sum / count};

    // rough r/p from raw acc
    double roll0 = std::atan2(-acc_mean(1), -acc_mean(2));
    double pitch0 = std::atan2(acc_mean(0), std::sqrt(acc_mean(1) * acc_mean(1) + acc_mean(2) * acc_mean(2)));
    double yaw0 = 0;

    // acc bias using gravity
    ecef2pos(pos_ecef.memptr(), pos_llh.memptr());
    Ceb = b2eDCM(pos_llh(0), pos_llh(1), roll0, pitch0, yaw0);
    arma::vec3 g_e = Gravity_ECEF(pos_ecef);
    arma::vec3 g_b = Ceb.t() * g_e;
    ACCbias_b = acc_mean + g_b;

    // gyro bias using earth rate
    const arma::vec3 om_ei_e = {0.0, 0.0, GNSS_OMEGA_EARTH_DOT};
    const arma::vec3 om_ei_b = Ceb.t() * om_ei_e;
    GYRbias_b = gyr_mean - om_ei_b;

    // final RPY (Groves-style) from bias-compensated means
    arma::vec3 acc_mean_bc = acc_mean - ACCbias_b;
    arma::vec3 gyr_mean_bc = gyr_mean - GYRbias_b;
    double roll = std::atan2(-acc_mean_bc(1), -acc_mean_bc(2));
    double pitch = std::atan2(acc_mean_bc(0), std::sqrt(acc_mean_bc(1) * acc_mean_bc(1) + acc_mean_bc(2) * acc_mean_bc(2)));
    double num = -gyr_mean_bc(1) * std::cos(roll) + gyr_mean_bc(2) * std::sin(roll);
    double den = gyr_mean_bc(0) * std::cos(pitch) + gyr_mean_bc(1) * std::sin(pitch) * std::sin(roll) + gyr_mean_bc(2) * std::cos(roll) * std::sin(pitch);
    double yaw = std::atan2(num, den);
    att_rpy = {roll, pitch, yaw};

    Ceb = b2eDCM(pos_llh(0), pos_llh(1), att_rpy(0), att_rpy(1), att_rpy(2));
}


void Inertial_Navigator::stepMechanizer(std::ifstream& fin_raw_imu)
{
    double imuPrevTime = obs.imuTime;
    readIMU(fin_raw_imu, false);
    double dT = obs.imuTime - imuPrevTime;
    if (!(dT > 0.0 && std::isfinite(dT))) return;

    // constants
    const arma::vec3 om_ei_e = {0.0, 0.0, GNSS_OMEGA_EARTH_DOT};

    // subtract biases
    arma::vec3 sf_b_ib = obs.Acc - ACCbias_b;
    arma::vec3 om_b_ib = obs.Gyr - GYRbias_b;

    // attitude update (first-order) | eq 5.27
    arma::mat Ceb_updt = Ceb * (arma::eye(3, 3) + (SkewMat(om_b_ib) * dT)) - SkewMat(om_ei_e) * Ceb * dT;
    // keep the DCM orthonormal (Newton–Schulz)
    Ceb_updt = 1.5 * Ceb_updt - 0.5 * Ceb_updt * (Ceb_updt.t() * Ceb_updt);

    // specific force | eq 5.28
    arma::vec3 sf_e_ib = (0.5 * (Ceb + Ceb_updt)) * sf_b_ib;

    // velocity update (with Coriolis) | eq 5.36
    arma::vec3 pos_ecef_mid = pos_ecef + vel_ecef * 0.5 * dT;
    arma::vec3 g_e = Gravity_ECEF(pos_ecef_mid);
    arma::vec3 acc_e = sf_e_ib + g_e - 2 * SkewMat(om_ei_e) * vel_ecef;
    arma::vec3 vel_ecef_updt = vel_ecef + acc_e * dT;

    // position update | eq 5.38
    arma::vec3 pos_ecef_updt = pos_ecef + (vel_ecef + vel_ecef_updt) * 0.5 * dT;

    // RPY update
    ecef2pos(pos_ecef_updt.memptr(), pos_llh.memptr());
    arma::mat Cne = e2llfDCM(pos_llh(0), pos_llh(1));
    arma::mat Cnb = Cne * Ceb_updt;
    arma::vec3 att_rpy_updt = dcm2euler(Cnb.t());
    att_rpy_updt(2) = normalise(att_rpy_updt(2), -M_PI, M_PI);

    // complementary filter to stabilize roll pitch
    double lambda = 0.95;
    double acc_roll = std::atan2(-sf_b_ib(1), -sf_b_ib(2));
    double acc_pitch = std::atan2(sf_b_ib(0), std::sqrt(sf_b_ib(1) * sf_b_ib(1) + sf_b_ib(2) * sf_b_ib(2)));
    att_rpy_updt(0) = lambda * att_rpy_updt(0) + (1 - lambda) * acc_roll;
    att_rpy_updt(1) = lambda * att_rpy_updt(1) + (1 - lambda) * acc_pitch;

    Ceb = Ceb_updt;
    pos_ecef = pos_ecef_updt;
    vel_ecef = vel_ecef_updt;
    att_rpy = att_rpy_updt;
    predictAntennaECEF();
    Fe = SkewMat(sf_e_ib);  // derivative of a rotated vector wrt a small rotation is a cross product (skew)
    Ne = GravGrad(pos_ecef);
}


void Inertial_Navigator::alignYawFromGnss(const arma::vec3& pos_ecef_gnss,
    const arma::vec3& vel_ecef_gnss,
    double yaw_boresight_rad)
{
    arma::vec3 p_llh_gnss;
    ecef2pos(pos_ecef_gnss.memptr(), p_llh_gnss.memptr());
    arma::mat Cne = e2llfDCM(p_llh_gnss(0), p_llh_gnss(1));
    arma::vec3 v_ned_gnss = Cne * vel_ecef_gnss;

    // smoothing
    v_ned_gnss_filt = (1.0 - ema_alpha) * v_ned_gnss_filt + ema_alpha * v_ned_gnss;

    const double vN = v_ned_gnss_filt(0);
    const double vE = v_ned_gnss_filt(1);
    const double speed_h = std::hypot(vN, vE);
    if (speed_h < 3) return;  // not enough horizontal motion

    pos_llh(0) = p_llh_gnss(0);
    pos_llh(1) = p_llh_gnss(1);
    pos2ecef(pos_llh.memptr(), pos_ant_ecef.memptr());  // snap horizontally only

    arma::vec3 v_ned_imu = Cne * vel_ant_ecef;
    v_ned_imu(0) = v_ned_gnss_filt(0);
    v_ned_imu(1) = v_ned_gnss_filt(1);
    vel_ant_ecef = Cne.t() * v_ned_imu;  // snap horizontally only
    predictIMUECEF();

    // yaw from course-over-ground
    const double yaw_course = std::atan2(vE, vN);  // heading from North (rad)
    att_rpy(2) = yaw_course + yaw_boresight_rad;
    att_rpy(2) = normalise(att_rpy(2), -M_PI, M_PI);

    // update attitude DCM
    ecef2pos(pos_ecef.memptr(), pos_llh.memptr());
    Ceb = b2eDCM(pos_llh(0), pos_llh(1), att_rpy(0), att_rpy(1), att_rpy(2));

    yaw_aligned_ = true;
}


void Inertial_Navigator::predictAntennaECEF()
{
    const arma::vec3 om_ei_e = {0.0, 0.0, GNSS_OMEGA_EARTH_DOT};
    arma::vec3 om_b_ib = obs.Gyr - GYRbias_b;

    arma::vec3 L_e = Ceb * Lxyz;
    arma::vec3 L_e_dot = Ceb * arma::cross(om_b_ib, Lxyz) - arma::cross(om_ei_e, L_e);

    // IMU origin -> Antenna
    pos_ant_ecef = pos_ecef + L_e;
    vel_ant_ecef = vel_ecef + L_e_dot;
}


void Inertial_Navigator::predictIMUECEF()
{
    const arma::vec3 om_ei_e = {0.0, 0.0, GNSS_OMEGA_EARTH_DOT};
    arma::vec3 om_b_ib = obs.Gyr - GYRbias_b;

    arma::vec3 L_e = Ceb * Lxyz;
    arma::vec3 L_e_dot = Ceb * arma::cross(om_b_ib, Lxyz) - arma::cross(om_ei_e, L_e);

    // Antenna -> IMU origin
    pos_ecef = pos_ant_ecef - L_e;
    vel_ecef = vel_ant_ecef - L_e_dot;
}


arma::mat Inertial_Navigator::SkewMat(const arma::vec3& Vec)
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


// Normalizes any number to an arbitrary range
// by assuming the range wraps around when going below min or above max
double Inertial_Navigator::normalise(const double value, const double start, const double end)
{
    const double width = end - start;
    const double offsetValue = value - start;

    return (offsetValue - (floor(offsetValue / width) * width)) + start;
}


// Numerical Jacobian of gravity (ECEF)
arma::mat Inertial_Navigator::GravGrad(const arma::vec3& r_e)
{
    arma::mat33 _Ne;
    _Ne.zeros();

    const double h = 1.0;  // meters
    const arma::vec3 ex{h, 0.0, 0.0};
    const arma::vec3 ey{0.0, h, 0.0};
    const arma::vec3 ez{0.0, 0.0, h};

    // Gravity at +/- steps along each axis
    const arma::vec3 g_xp = Gravity_ECEF(r_e + ex);
    const arma::vec3 g_xm = Gravity_ECEF(r_e - ex);
    const arma::vec3 g_yp = Gravity_ECEF(r_e + ey);
    const arma::vec3 g_ym = Gravity_ECEF(r_e - ey);
    const arma::vec3 g_zp = Gravity_ECEF(r_e + ez);
    const arma::vec3 g_zm = Gravity_ECEF(r_e - ez);

    // Columns: dg/dx, dg/dy, dg/dz
    _Ne.col(0) = (g_xp - g_xm) / (2.0 * h);
    _Ne.col(1) = (g_yp - g_ym) / (2.0 * h);
    _Ne.col(2) = (g_zp - g_zm) / (2.0 * h);

    return _Ne;
}
