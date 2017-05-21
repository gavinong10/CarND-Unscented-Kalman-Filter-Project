#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
//  std_a_ = 30;
  std_a_ = 4;
//  std_a_ = 0.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
//  std_yawdd_ = 30;
    std_yawdd_ = 1;
//  std_yawdd_ = 0.5;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */


//  ///* State dimension
  n_x_ = 5;

//  ///* Augmented state dimension
  n_aug_ = 7;

//  ///* Sigma point spreading parameter
  lambda_ = 3 - n_x_;

//  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
//    x_ = VectorXd(5);
//    P_ = MatrixXd(5, 5);

    switch (meas_package.sensor_type_) {
      case MeasurementPackage::LASER: {
        if (!use_laser_) {
          return;
        }
      }
      case MeasurementPackage::RADAR: {
        if (!use_radar_) {
          return;
        }
      }
    }


    cout << "Initializing" << endl;

    switch (meas_package.sensor_type_) {

      // X is a vector of 5:
      // px
      // py
      // v
      // yaw
      // yawd


      case MeasurementPackage::LASER: {
        // Laser has vector of (px, py)
        auto px = meas_package.raw_measurements_[0];
        auto py = meas_package.raw_measurements_[1];

        x_ << px, py, 0, 0, 0;
        break;
      }
      case MeasurementPackage::RADAR: {
        // Radar has vector of (rho, phi, rhodot)
        auto rho = meas_package.raw_measurements_[0];
        auto phi = meas_package.raw_measurements_[1];
        auto rhodot = meas_package.raw_measurements_[2];

        x_ << rho * cos(phi), rho * sin(phi), rhodot, phi, 0;
        break;
      }
      default: {
        throw "Unknown Measurement Package for initialization!";
      }

    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Prediction
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  cout << "Making Prediction" << endl;
  Prediction(dt);

  // UpdateLidar or UpdateRadar
  switch (meas_package.sensor_type_) {
    case MeasurementPackage::LASER: {
      cout << "Updating with Lidar measurement" << endl;
      UpdateLidar(meas_package);
      break;
    }
    case MeasurementPackage::RADAR: {
      cout << "Updating with Radar measurement" << endl;
      UpdateRadar(meas_package);
      break;
    }
    default: {
      throw "Unknown Measurement Package for update!";
    }

  }
  time_us_ = meas_package.timestamp_;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Generate sigma points
  auto n_sigma = 2 * n_aug_ + 1;

#ifdef DEBUG
  cout << "n_sigma" << n_sigma << endl;
#endif
  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);

  x_aug.head(n_x_) = x_;

#ifdef DEBUG
  cout << "Creating P_aug" << endl;
#endif
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;

#ifdef DEBUG
  cout << "Creating Q" << endl;
#endif
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_ * std_a_, 0,
      0, std_yawdd_ * std_yawdd_;

  P_aug.block(P_.rows(), P_.cols(), Q.rows(), Q.cols()) = Q;

#ifdef DEBUG
  cout << "P_aug has been set" << endl;
#endif

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma);
  Xsig_aug.col(0) = x_aug;

  // Calculate sqrt of matrix P:
  MatrixXd A = P_aug.llt().matrixL();
  float c = sqrt(lambda_ + n_aug_);
  MatrixXd scaled_A = c * A;

#ifdef DEBUG
  cout << "About to populate Xsig_aug" << endl;
  cout << "Xsig_aug dims " << Xsig_aug.rows() << ", " << Xsig_aug.cols() << endl;
  cout << "scaled_A dims " << scaled_A.rows() << ", " << scaled_A.cols() << endl;

#endif
  for (int i = 0; i < scaled_A.cols(); i++) {
    Xsig_aug.col(i + 1) = x_aug + scaled_A.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - scaled_A.col(i);
  }

  // Make sigma point prediction
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
#ifdef DEBUG
  cout << "Xsig_pred_ defined" << endl;
#endif
  for (int i = 0; i < Xsig_aug.cols(); i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    //predicted state values
    double px_p, py_p;

    if (fabs(yawd) > 0.001) {
      px_p = px + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = py + v / yawd * (-cos(yaw + yawd * delta_t) + cos(yaw));
    } else {
      px_p = px + v * cos(yaw) * delta_t;
      py_p = py + v * sin(yaw) * delta_t;
    }

    // Add noise
    px_p += 0.5 * delta_t * delta_t * nu_a * cos(yaw);
    py_p += 0.5 * delta_t * delta_t * nu_a * sin(yaw);

    double v_p = v + delta_t * nu_a;
    double yaw_p = yaw + yawd * delta_t + 0.5 * delta_t * delta_t * nu_yawdd;
    double yawd_p = yawd + delta_t * nu_yawdd;

    //write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // predict mean and covariance
  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  //Normalize angle
  while (x_(3) > M_PI) x_(3)  -= 2. * M_PI;
  while (x_(3) < M_PI) x_(3)  += 2. * M_PI;

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

#ifdef DEBUG
  cout << "Finished Prediction" << endl;
#endif
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //set measurement dimension, lidar can measure x and y
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;                                 //phi
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;
  S = S + R;


  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //Normalize angle
    while (x_diff(3) > M_PI) x_diff(3)  -= 2. * M_PI;
    while (x_diff(3) < M_PI) x_diff(3)  += 2. * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  // Normalize angle
  while (x_(3) > M_PI) x_(3)  -= 2. * M_PI;
  while (x_(3) < M_PI) x_(3)  += 2. * M_PI;

  P_ = P_ - K * S * K.transpose();

  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

#ifdef DEBUG
  cout << "Finished Lidar update" << z_pred << endl;
#endif

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  cout << "Xsig_pred_ dims " << Xsig_pred_.rows() << ", " << Xsig_pred_.cols() << endl;


  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                        //r
    Zsig(1, i) = atan2(p_y, p_x);                                 //phi
    Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  while (z_pred(1) > M_PI) z_pred(1) -= 2. * M_PI;
  while (z_pred(1) < -M_PI) z_pred(1) += 2. * M_PI;


#ifdef DEBUG
  cout << "z_pred defined: " << z_pred << endl;
#endif
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;
  S = S + R;


  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

#ifdef DEBUG
  cout << "cross correlation matrix calculated" << endl;
#endif

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  //Normalize angle
  while (x_(3) > M_PI) x_(3)  -= 2. * M_PI;
  while (x_(3) < M_PI) x_(3)  += 2. * M_PI;
  P_ = P_ - K * S * K.transpose();

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

#ifdef DEBUG
  cout << "finished UpdateRadar" << endl;
#endif
}
