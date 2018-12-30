#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  // Initialization flag:
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // Initializing R & H matrices:
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // Measurement covariance matrix (Laser):
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // Measurement covariance matrix (Radar):
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  
  // Initializing the measurement function matrix H (Laser):
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Initializing the state covariance matrix P:
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
  
  // Initializing the state transition matrix F:
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // Initializing the process covariance matrix Q:
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             1, 0, 1, 0,
             0, 1, 0, 1;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    
    // Initializing the state vector x:
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      cout << "EKF's first measurement (RADAR):" << endl;

      // Extracting the range, bearing and radial velocity from the first measurement:
      double rho = measurement_pack.raw_measurements_[0];     // Range
      double phi = measurement_pack.raw_measurements_[1];     // Bearing (angle)
      double rho_dot = measurement_pack.raw_measurements_[2]; // Radial velocity

      // Polar -> Cartesian transformation:
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);

      // Storing the catesian coordinates in our state vector:
      ekf_.x_ << x, y, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      cout << "EKF's first measurement (LIDAR):" << endl;
      
      // Storing the Laser's Px & Py, alongside V = 0 in our state vector:
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // First measurement's timestamp:
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << ekf_.x_ << endl;

    cout << "INITIALIZING COMPLETE!" << endl;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // Compute the elapsed time between current and previous measurements (in seconds):
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Updating the state transition matrix F:
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // Setting the acceleration (x, y) variance:
  double noise_ax = 9;
  double noise_ay = 9;

  // Cutting corners (for calculating Q):
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;
  double dt_3_2 = dt_3 / 2;
  double dt_4_4 = dt_4 / 4;
  
  // Updating the process covariance matrix Q:
  ekf_.Q_ << dt_4_4*noise_ax, 0, dt_3_2*noise_ax, 0,
             0, dt_4_4*noise_ay, 0, dt_3_2*noise_ay,
             dt_3_2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3_2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    // Upadting the H and R matrices to be the radar's:
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // TODO: Laser updates
    // Upadting the H and R matrices to be the laser's:
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
