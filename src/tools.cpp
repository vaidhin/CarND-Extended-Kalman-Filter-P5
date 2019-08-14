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
	rmse << 0,0,0,0;

	if(estimations.size() != ground_truth.size()
		|| estimations.size() == 0){
	cout << "Invalid estimation or ground_truth data" << endl;
	return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  	MatrixXd Hj(3,4);
	
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero
	if((px == 0) & (py == 0)) {
	    cout << "CalcuateJacobian() - Error division by zero" << endl;
	} else {
	    //Calculating repeated elements
	    float xysq = (px * px) + (py * py);
	    float sqtxysq = sqrt(xysq); //Square root
	    float xy23 = xysq * sqtxysq; // power (3/2) = x * sqrt(x)
	    

	    Hj << px/sqtxysq,   py/sqtxysq, 0, 0,
	          -py/xysq,  px/xysq, 0, 0,
	          (py * (vx*py - vy*px))/xy23, (px * (vy*px - vx*py))/xy23, px/sqtxysq, py/sqtxysq;
	}
	
	return Hj;
}
