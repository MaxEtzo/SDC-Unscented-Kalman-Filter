#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	/**
	 * TODO:
	 * Calculate the RMSE here.
	 */

	VectorXd rmse(4);
	rmse << 0,0,0,0;
	unsigned int est_size = estimations.size();	
	if (est_size == 0 || est_size != ground_truth.size())
	{
		cout << "Invalid size of estimations vector" << endl;
		return rmse;
	}   

	//accumulate squared residuals
	for(unsigned int i=0; i < est_size; ++i){
		VectorXd residual = (estimations[i] - ground_truth[i]).array().square();
		rmse += residual;
	}

	//calculate the mean
	rmse /= est_size;
	//calculate the squared root
	rmse = rmse.array().sqrt();
	//return the result
	return rmse;
}
