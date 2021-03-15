#include "tools.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  double c1 = pow(px,2) + pow(py,2);
  // check division by zero
  if (fabs(c1) < 0.0001) {
//    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
//    return Hj;
    px += 0.001;
    py += 0.001;
    c1 = pow(px,2) + pow(py,2);
  }
  double c2 = sqrt(c1);
  double c3 = c1*c2;

  // compute the Jacobian matrix
  Hj << px/c2, py/c2, 0, 0,
      -py/c1, px/c1, 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}

VectorXd Tools::xy2polar(const VectorXd& x_state){
// calculate nonlinear h(x)
// convert from x-y coordinate to polar coordinate
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  double rho = sqrt(pow(px,2) + pow(py,2));
  double angle = atan2(py, px);
  double rho_derv = (px*vx + py*vy) / rho;

  VectorXd states_polar = VectorXd(3);
  states_polar << rho, angle, rho_derv;
  return states_polar;  
}

double Tools::normalAngle(const double angle_raw){
// Normalize angle in polar axis -- output in range [-pi, pi)
  double angle_offsetpi = std::remainder(angle_raw + M_PI, 2*M_PI);
  double angle = (angle_offsetpi>=0) ? (angle_offsetpi-M_PI) : (angle_offsetpi+M_PI);
  return angle;
}

VectorXd Tools::polar2xy(const VectorXd& x_state){
// convert from polar coordinate to x-y
  double rho = x_state(0);
  double angle = x_state(1);
  double rho_derv = x_state(2);
  
  double sin_angle = sin(angle);
  double cos_angle = cos(angle);
  
  double px = cos_angle * rho;
  double py = sin_angle * rho;
  double vx = cos_angle * rho_derv;
  double vy = sin_angle * rho_derv;

  VectorXd states_xy = VectorXd(4);
  states_xy << px, py, vx, vy;
  return states_xy;
}
