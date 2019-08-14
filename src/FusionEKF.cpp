#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  ekf_.P_  = MatrixXd(4, 4);
  Hj_      = MatrixXd(3, 4);
  ekf_.F_  = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        	  0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
       		  0, 0.0009, 0,
       		  0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
  			  0, 1, 0, 0;



  //H_radar is jacobian, will be calculated later below

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    ekf_.P_ << 1, 0, 0, 0,
        	   0, 1, 0, 0,
        	   0, 0, 1000, 0,
        	   0, 0, 0, 1000;

  	ekf_.F_ << 1, 0, 1, 0,
       		   0, 1, 0, 1,
 			   0, 0, 1, 0,
 			   0, 0, 0, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
  		float ro     = measurement_pack.raw_measurements_[0];
  		float theta  = measurement_pack.raw_measurements_[1];
  		float ro_dot = measurement_pack.raw_measurements_[2];

    	float px_r = ro * cos(theta);
    	float py_r = ro * sin(theta);
    	float vx_r = ro_dot * cos(theta);
    	float vy_r = ro_dot * sin(theta);

    	ekf_.x_ << px_r, py_r, vx_r, vy_r;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    	//Get x, y positions from first measurement, set velocity co-ordinates to 0's since LIDAR cannot track velocity
    	ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0;

    }

   	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
 	long long current_timestamp_ = measurement_pack.timestamp_;
 	float dt = (current_timestamp_ - previous_timestamp_) / 1000000.0;
 	previous_timestamp_ = current_timestamp_;

 	/**
 		Set Predict Matrices - State transition function F & Process covariance Q
 		These will be same irrespective of LIDAR or RADAR
 	*/
 	

 	ekf_.F_(0, 2) = dt;
 	ekf_.F_(1, 3) = dt;
 	
 	float dt2 = dt * dt;
	float dt3 = dt2 * dt;
	float dt4 = dt3 * dt;
	float noise_ax = 9;
	float noise_ay = 9;

	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << dt4/4*noise_ax, 0, dt3/2*noise_ax, 0,
	         0, dt4/4*noise_ay, 0, dt3/2*noise_ay,
	         dt3/2*noise_ax, 0, dt2*noise_ax, 0,
	         0, dt3/2*noise_ay, 0, dt2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    
  	ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
  	ekf_.H_ = H_laser_;
  	ekf_.R_ = R_laser_;
  	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
