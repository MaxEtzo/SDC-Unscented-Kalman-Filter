#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

    lidar_NIS.open("../lidar_NIS");
    radar_NIS.open("../radar_NIS");

	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	/// State dimension
	n_x_ = 5;

	/// Augmented state dimension
	n_aug_ = n_x_ + 2;

	/// Number of sigma points
	n_sig_ = 2 * n_aug_ + 1;

	/// Sigme point spreading parameter
	lambda_ = 3 - n_aug_;

	// initial state vector
	x_ = VectorXd::Zero(n_x_);

	// initial covariance matrix
	P_ = MatrixXd::Zero(n_x_, n_x_);
	for (int i = 0; i < n_x_; i++)
		P_(i,i) = 1.0;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 1.5;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = M_PI / 6;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

    // Radar measurement covariance matrix
    R_radar_ = MatrixXd::Zero(3,3);
    R_radar_(0,0) = std_radr_ * std_radr_;
    R_radar_(1,1) = std_radphi_ * std_radphi_;
    R_radar_(2,2) = std_radrd_ * std_radrd_;

    // Radar measurement covariance matrix
    R_laser_ = MatrixXd::Zero(2,2);
    R_laser_(0,0) = std_laspx_ * std_laspx_;
    R_laser_(1,1) = std_laspy_ * std_laspy_; 

	/// initialize weights
	weights_ = VectorXd(n_sig_);
	weights_I = MatrixXd::Zero(n_sig_, n_sig_);
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	weights_I(0,0) = weights_(0);
	for (int i = 1; i < 2 * n_aug_ + 1; i++){
		weights_(i) = 0.5 / (lambda_ + n_aug_);
		weights_I(i,i) = weights_(i);
	}
	
	/// Predicted sigma points initialization 
	Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig_);
    
    /// Predicted sigma points error
    Xsig_err_ = MatrixXd::Zero(n_x_, n_sig_);
}

UKF::~UKF() {
    lidar_NIS.close();
    radar_NIS.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    // first measurement 
	if(!is_initialized_){
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
			/**
			  Convert radar from polar to cartesian coordinates and initialize state.
			  */
			float cos_phi = cos(meas_package.raw_measurements_(1));
			float sin_phi = sin(meas_package.raw_measurements_(1));
			// estimate cartesian state
			x_(0) = meas_package.raw_measurements_(0) * cos_phi;
			x_(1) = meas_package.raw_measurements_(0) * sin_phi;
			P_(0,0) = R_radar_(0,0) + R_radar_(1,1);
			P_(1,1) = P_(0,0);
			time_us_ = meas_package.timestamp_;
			is_initialized_ = true;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
			/**
			  Initialize state.
			  */
			x_(0) = meas_package.raw_measurements_(0);
			x_(1) = meas_package.raw_measurements_(1);
			P_(0,0) = R_laser_(0,0); 
			P_(1,1) = R_laser_(1,1);
			time_us_ = meas_package.timestamp_;
			is_initialized_ = true;
		}
		return;
	}
    bool is_radar = meas_package.sensor_type_ == MeasurementPackage::RADAR;
	if ((is_radar && use_radar_) || (!is_radar && use_laser_))
	{
		double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
		time_us_ = meas_package.timestamp_;
		Prediction(delta_t);
        if (is_radar)
            UpdateRadar(meas_package);
        else
            UpdateLidar(meas_package);
	}
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    
	MatrixXd X_sig;
	AugmentedSigmaPoints(&X_sig);
	PredictSigmaPoints(delta_t, X_sig);

	//create vector for predicted state
	x_ = Xsig_pred_ * weights_;

	//create covariance matrix for prediction
	Xsig_err_ = Xsig_pred_.colwise() - x_;

	for (int i = 0; i < n_sig_; i++)
        Xsig_err_(3,i) = NormalizeAngle(Xsig_err_(3,i));

	P_ = Xsig_err_ * weights_I * Xsig_err_.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    int n_z = 2;
    MatrixXd Z_pred = Xsig_pred_.block(0,0,n_z,n_sig_);
    UpdateCommon(meas_package, Z_pred);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	int n_z = 3;
	//create matrix for sigma points in measurement space
	MatrixXd Z_pred = MatrixXd::Zero(n_z, n_sig_);
	for (int i = 0; i < n_sig_; i++)
	{
		double px = Xsig_pred_(0,i);
		double py = Xsig_pred_(1,i);
		double v  = Xsig_pred_(2,i);
		double yaw = Xsig_pred_(3,i);
		// polar coordinates (radar) state vector
		Z_pred(0,i) = sqrt(px*px + py*py);
		Z_pred(1,i) = atan2(py, px);
		if (fabs(Z_pred(0,i)) > 0.001)
            Z_pred(2,i) = (px*cos(yaw) + py*sin(yaw))*v / Z_pred(0,i);
        else
            Z_pred(2,i) = 0;
	}
    UpdateCommon(meas_package, Z_pred);
}

/**
 * UpdateCommon common update steps for both radar and lidar
 * @param meas_package measurement package (either radar or lidar)
 * @param Z_pred predicted state 
 */
void UKF::UpdateCommon(MeasurementPackage meas_package, MatrixXd Z_pred)
{
    bool is_radar = meas_package.sensor_type_ == MeasurementPackage::RADAR;

	// mean predicted measurement
	VectorXd z_mean = Z_pred * weights_;
	MatrixXd Z_err = Z_pred.colwise() - z_mean;
    if (is_radar)
        for (int i = 0; i < n_sig_; i++)
            Z_err(1,i) = NormalizeAngle(Z_err(1,i));

    // measurement covariance matrix S
    MatrixXd S = Z_err * weights_I * Z_err.transpose();
    S = S + (is_radar? R_radar_ : R_laser_);	
    // cross correlation matrix Tc
    MatrixXd Tc = Xsig_err_ * weights_I * Z_err.transpose();
    // Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    // update state mean and covariance matrix
    VectorXd z_diff = meas_package.raw_measurements_ - z_mean;
    if (is_radar)
        z_diff(1) = NormalizeAngle(z_diff(1));
    x_ = x_ + K * z_diff; 
    P_ = P_ - K * S * K.transpose();

    // calculate and write NIS
    double eps = z_diff.transpose() * S.inverse() * z_diff;
    if (is_radar)
        radar_NIS << eps << endl;
    else lidar_NIS << eps << endl;
}

/**
 * AugmentedSigmaPoints returns augmented sigma points
 * @param Xsig_out a pointer for result augmented sigma points
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out)
{
    //create augmented mean vector
    VectorXd x_aug = VectorXd::Zero(n_aug_);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

    //create augmented sigma points
    MatrixXd X_sig = MatrixXd::Zero(n_aug_, n_sig_); 

    //create augmented mean state
    x_aug.head(n_x_) = x_;

    //create augmented covariance matrix
    P_aug.topLeftCorner(n_x_,n_x_) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;

    //create square root matrix
    MatrixXd A = P_aug.llt().matrixL();

    //create augmented sigma points
    X_sig.col(0)  = x_aug;
    double scale_fact = sqrt(lambda_ + n_aug_);
    for (int i = 0; i < n_aug_; i++)
    {
        X_sig.col(i+1)     	= x_aug + scale_fact * A.col(i);
        X_sig.col(i+1+n_aug_) = x_aug - scale_fact * A.col(i);
    }
    *Xsig_out = X_sig;
}

/**
 * SigmaPointsPrediction assigns predicted sigma points
 * @param delta_t Time between k and k+1 in s
 * @param X_sig augmented sigma points 
 */
void UKF::PredictSigmaPoints(double delta_t, MatrixXd X_sig) {

    Xsig_pred_ = X_sig.block(0,0,n_x_,n_sig_); 
    double dt = delta_t;
    double dt2 = dt * dt / 2;
    for (int i = 0; i < n_sig_; i++)
    {
        // reading necessary data
        double v   = X_sig(2,i);
        double yaw = X_sig(3,i);
        double yawd = X_sig(4,i);
        double nu_a = X_sig(5,i);
        double nu_yawdd = X_sig(6,i);
        double cos_yaw = cos(yaw);
        double sin_yaw = sin(yaw);
        // error term
        VectorXd x_err = VectorXd(5);
        x_err << dt2 * cos_yaw * nu_a, dt2 * sin_yaw * nu_a, dt * nu_a, dt2 * nu_yawdd, dt * nu_yawdd;
        Xsig_pred_.col(i) += x_err;
        if (fabs(yawd) > 0.001)
        {
            Xsig_pred_(0,i) += v / yawd * (sin(yaw + yawd * dt) - sin_yaw);
            Xsig_pred_(1,i) += v / yawd * (cos_yaw - cos(yaw + yawd * dt));
            Xsig_pred_(3,i) += yawd * dt;
        }
        else
        {
            Xsig_pred_(0,i) += v * cos_yaw * dt;
            Xsig_pred_(1,i) += v * sin_yaw * dt;
        }
    }
}

double UKF::NormalizeAngle(double A) {
    double A_n = A;
    while(A_n >= M_PI) A_n -= 2 * M_PI;
    while(A_n <= -M_PI) A_n += 2 * M_PI;
    return A_n;
}

