#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
  const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;
    if (estimations.size() != ground_truth.size() || estimations.size() == 0)
    return rmse;
    
    for(int i = 0; i < estimations.size(); i++) {
      VectorXd diff = estimations[i] - ground_truth[i];
      diff = diff.array() * diff.array();
      rmse += diff;
    }
    
    rmse /= estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
  }
  
  MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    float px = x_state[0];
    float py = x_state[1];
    float vx = x_state[2];
    float vy = x_state[3];
    
    float c1 = px * px + py * py;
    float c2 = sqrt(c1);
    float c3 = (c1 * c2);
    float c4 = vx * py - vy * px;
    
    MatrixXd Hj(3, 4);
    
    if((px == 0.0f && py == 0.0f) || c1 < 0.0001f)
    return Hj;
    
    Hj << px / c2,      py / c2,       0,       0,
          -py / c1,     px / c1,       0,       0,
          py * c4 / c3, -px * c4 / c3, px / c2, py / c2;
    return Hj;
  }
  