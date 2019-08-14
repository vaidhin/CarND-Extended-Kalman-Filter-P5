#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  KF_basic(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
      y = z - h(x)

  */
    float px1 = x_(0);
    float py1 = x_(1);
    float vx1 = x_(2);
    float vy1 = x_(3);

    /*
    MISTAKE :: used px*vy + vy*px instead of px*vx + py*vy
    */

    float rec1 = sqrt(px1*px1 + py1*py1);
    float rec2 = atan2(py1, px1);
    float rec3 = (px1*vx1 + py1*vy1) / rec1;



    VectorXd z_pred = VectorXd(3);
    z_pred << rec1, rec2, rec3;
    VectorXd y = z - z_pred;
   


    while(y(1) > M_PI || y(1) < -M_PI) {
      if(y(1) > M_PI) {
        y(1) -= 2 * M_PI;
      }
      if(rec2 < -M_PI) {
        y(1) += 2 * M_PI;
      }
    }
    KF_basic(y);
}

void KalmanFilter::KF_basic(const VectorXd &y) {
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
