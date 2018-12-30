#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */

  // The prediction steps as learned in the lessons:
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  
  // KF Measurement update step:
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  UpdatePostY(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  // Recover state parameters:
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  
  // Pre-compute a set of terms to avoid repeated calculations:
  double rho = sqrt(px*px + py*py);        // Range
  double phi = atan2(py, px);              // Bearing (angle)
  double rhodot = (px*vx + py*vy) / rho;   // Radial velocity

  //// EKF Measurement update step:
  
  // Creating the h(x') vector:
  VectorXd h_x = VectorXd(3);
  h_x << rho, phi, rhodot;

  // Calculating the error vector y, while keeping the bearing in the [-pi, pi] range:
  VectorXd y = z - h_x;
  while ( y(1) > M_PI || y(1) < -M_PI) {
    if ( y(1) > M_PI ) {
      y(1) -= M_PI;
    }
    else if ( y(1) < -M_PI ) {
      y(1) += M_PI;
    }
  }

  UpdatePostY(y);
}

void KalmanFilter::UpdatePostY(const VectorXd &y) {
  /**
   * Updates the state post calculating the y (error) vector.
   */

  // The measurement update steps (for laser & radar) as learned in the lessons:
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // New estimate:
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}