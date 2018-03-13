#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.54;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/3.0;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
 
  n_x_ = 5;
  n_aug_ = 7;

  lambda_ = 3 - n_aug_;
  weights_ = VectorXd::Zero(2*n_aug_ + 1);
  /// set weights
  weights_(0) = double(lambda_/(lambda_ + double(n_aug_)));
  for (int i = 1; i < 2*n_aug_ + 1; ++i) {
    weights_(i) = double(1.0/(2.0*(lambda_ + double(n_aug_))));
  }
}

UKF::~UKF() {
  if (nis_file_.is_open()) {
    nis_file_.close();
  }
}

/** init_helper_
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::init_helper_(MeasurementPackage meas_package) {
  time_us_ = meas_package.timestamp_;
  P_(0, 0) = 0.5;
  P_(1, 1) = 0.5;
  P_(2, 2) = 0.5;
  P_(3, 3) = 0.5;
  P_(4, 4) = 0.5;
  //x_ << 1.0,1.0,0.15,0.1,M_PI/20.0;
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    x_(0) = meas_package.raw_measurements_(0);
    x_(1) = meas_package.raw_measurements_(1);
    is_initialized_ = true;
    return;
  }
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    double x = meas_package.raw_measurements_(0) * 
                             cos(meas_package.raw_measurements_(1));    
    double y = meas_package.raw_measurements_(0) *
                                   sin(meas_package.raw_measurements_(1));    

    x_(0) = x;
    x_(1) = y;
    is_initialized_ = true;
    return;
  }
  std::cout << "Error in input..." << std::endl;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) {
    return;
  }
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) {
    return;
  }
  if (!is_initialized_) {
    return init_helper_(meas_package);
  }
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else {
    std::cout << "Input error..." << std::endl;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  Xsig_aug_ = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);
  VectorXd x_aug = VectorXd::Zero(7);
  MatrixXd P_aug = MatrixXd::Zero(7, 7);
  // setting first 5 elements as x_
  x_aug.head(x_.size()) = x_;
  // setting next two elements as 0.
  x_aug(x_.size()) = 0;
  x_aug(x_.size() + 1) = 0;

  // Creating P_aug with
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  Xsig_aug_.col(0) = x_aug;

  MatrixXd sqmat = P_aug.llt().matrixL();

  for (int i = 1; i < n_aug_ + 1; i++) {
    Xsig_aug_.col(i) = x_aug + sqrt(lambda_ + n_aug_) * sqmat.col(i - 1);
    Xsig_aug_.col(n_aug_ + i) = x_aug - sqrt(lambda_ + n_aug_) * sqmat.col(i - 1);
  }
  // Now we have generated representative sigma points.
  // Next, predict the new sigma points.
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    double px_p, py_p;

    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // Based on the predicted sigma points, now generate the mean and co-var matrix
  /// Predict mean
  x_.fill(0.0);
  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    x_ = x_ + Xsig_pred_.col(i) * weights_(i);
  }
  P_ = MatrixXd::Zero(5, 5);
  /// Predict variance
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI){ x_diff(3)-=2.*M_PI;}
    while (x_diff(3)<-M_PI){ x_diff(3)+=2.*M_PI;}
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  MatrixXd Zsig = Xsig_pred_.block(0, 0, lid_dim_,
            2 * n_aug_ + 1);
  MatrixXd R = MatrixXd::Zero(lid_dim_, lid_dim_);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  update_helper_(Zsig, lid_dim_, meas_package, R);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  MatrixXd Zsig = MatrixXd::Zero(rad_dim_, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        
    Zsig(1,i) = atan2(p_y,p_x);                                 
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   
  }
  MatrixXd R = MatrixXd::Zero(rad_dim_, rad_dim_);
  R <<    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;
  update_helper_(Zsig, rad_dim_, meas_package, R);
}

void UKF::update_helper_(const MatrixXd &Zsig,
                         const int sensor_dim,
                         const MeasurementPackage &meas_package,
                         const MatrixXd &R) {
  VectorXd z_pred = VectorXd::Zero(sensor_dim);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  MatrixXd S = MatrixXd::Zero(sensor_dim, sensor_dim);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R;

  MatrixXd Tc = MatrixXd::Zero(n_x_, sensor_dim);
  VectorXd z = VectorXd::Zero(sensor_dim);
  for (int i = 0; i < sensor_dim; ++i) {
    z(i) = meas_package.raw_measurements_(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;


    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd Sinv = S.inverse();
  MatrixXd K = Tc * Sinv;
  VectorXd z_diff = z - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  double nis_val = z_diff.transpose() * Sinv * z_diff;
  nis_file_.open("nis_vals.txt", ios::app);
  nis_file_ << nis_val << "\n";
  nis_file_.close();
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
