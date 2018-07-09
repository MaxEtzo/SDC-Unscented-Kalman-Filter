#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
	private:
		/// Weights in Identity format
		MatrixXd weights_I;
		/**
		 * AugmentedSigmaPoints generates augmented sigma points
		 * @param Xsig_out a pointer for result augmented sigma points
		 */
		void AugmentedSigmaPoints(MatrixXd* Xsig_out); 

		/**
		 * SigmaPointsPrediction assigns predicted sigma points
		 * @param delta_t Time between k and k+1 in s
		 * @param Xsig_aug augmented sigma points 
		 * @param Xsig_out a pointer for result predicted sigma points
		 */
		void SigmaPointsPrediction(double delta_t, MatrixXd Xsig_aug, MatrixXd* Xsig_out);

		/**
		 * PredictRadarMeasurement predicts radar measurement
		 * @param z_out pointer for predicted radar state
		 * @param Zsig_out pointer for sigma points in radar space
		 * @param Zerr_out pointer for error matrix in radar space
		 * @param S_out pointer for measurement covariance matrix 
		 */
		void PredictRadarMeasurement(VectorXd* z_out, MatrixXd* Zsig_out, MatrixXd* Zerr_out, MatrixXd* S_out);
	public:

		///* initially set to false, set to true in first call of ProcessMeasurement
		bool is_initialized_;

		///* if this is false, laser measurements will be ignored (except for init)
		bool use_laser_;

		///* if this is false, radar measurements will be ignored (except for init)
		bool use_radar_;

		///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
		VectorXd x_;

		///* state covariance matrix
		MatrixXd P_;

		///* predicted sigma points matrix
		MatrixXd Xsig_pred_;

		///* time when the state is true, in us
		long long time_us_;

		///* Process noise standard deviation longitudinal acceleration in m/s^2
		double std_a_;

		///* Process noise standard deviation yaw acceleration in rad/s^2
		double std_yawdd_;

		///* Laser measurement noise standard deviation position1 in m
		double std_laspx_;

		///* Laser measurement noise standard deviation position2 in m
		double std_laspy_;

		///* Radar measurement noise standard deviation radius in m
		double std_radr_;

		///* Radar measurement noise standard deviation angle in rad
		double std_radphi_;

		///* Radar measurement noise standard deviation radius change in m/s
		double std_radrd_ ;

		///* Weights of sigma points
		VectorXd weights_;

		///* State dimension
		int n_x_;
		
		///* Augmented state dimension
		int n_aug_;

		///* Sigma point spreading parameter
		double lambda_;

		/**
		 * Constructor
		 */
		UKF();

		/**
		 * Destructor
		 */
		virtual ~UKF();

		/**
		 * ProcessMeasurement
		 * @param meas_package The latest measurement data of either radar or laser
		 */
		void ProcessMeasurement(MeasurementPackage meas_package);

		/**
		 * Prediction Predicts sigma points, the state, and the state covariance
		 * matrix
		 * @param delta_t Time between k and k+1 in s
		 */
		void Prediction(double delta_t);

		/**
		 * Updates the state and the state covariance matrix using a laser measurement
		 * @param meas_package The measurement at k+1
		 */
		void UpdateLidar(MeasurementPackage meas_package);

		/**
		 * Updates the state and the state covariance matrix using a radar measurement
		 * @param meas_package The measurement at k+1
		 */
		void UpdateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
