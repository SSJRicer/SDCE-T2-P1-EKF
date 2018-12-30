#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::cout;
using std::endl;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */

  // Initializing the RMSE vector:
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  unsigned int est_size = estimations.size();

  // Check the validity of the inputs:
  //  - Estimation vector size shouldn't be 0.
  //  - Estimation vector size should be equal to ground truth vector size.
  if ( est_size == 0 || est_size != ground_truth.size() ) {
     cout << "Invalid estimation or ground_truth data!" << endl;
     return rmse;
  }

  // Accumulate squared residuals:
  for (unsigned int i = 0; i < est_size; ++i) {
     VectorXd residual = estimations[i] - ground_truth[i];
     residual = residual.array() * residual.array();
     rmse += residual;
  }

  // Calculate the mean:
  rmse /= est_size;

  // Calculate the squared root:
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  
  // Initializing the H Jacobian matrix:
  MatrixXd Hj(3, 4);
  Hj << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

  // Recover state parameters:
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // Pre-compute a set of terms to avoid repeated calculations:
  double r_2 = px*px + py*py;  // Squared position radius
  double r = sqrt(r_2);        // Position radius
  double r_3_2 = r * r_2;      // Radius^(3/2)


  // Check division by 0 (or close to 0):
  if ( fabs(r_2) < 0.0001) {
     cout << "CalculateJacobian() - Error - Division by Zero!" << endl;
     return Hj;
  }

  // Compute the Jacobian matrix:
  Hj << (px / r), (py / r), 0, 0,
        -(py / r_2), (px / r_2), 0, 0,
        py*(vx*py-vy*px)/r_3_2, px*(px*vy-py*vx)/r_3_2, px/r, py/r;
   
  return Hj;
}
