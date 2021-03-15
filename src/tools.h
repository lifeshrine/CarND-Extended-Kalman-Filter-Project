#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {
 public:
  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * A helper method to calculate RMSE.
   */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, 
                                const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);
  
  /**
   * A helper method to convert from xy coordinate to polar.
   */
  Eigen::VectorXd xy2polar(const Eigen::VectorXd& x_state);
  
  /**
   * A helper method to standardize polar axis angle.
   */
  double normalAngle(const double angle_raw);
  
  /**
   * A helper method to convert from polar coordinate to xy.
   */
  Eigen::VectorXd polar2xy(const Eigen::VectorXd& x_state);
};

#endif  // TOOLS_H_
