#include "ukf.h"
#include "Eigen/Dense"
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
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .5;

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
  is_init = false;

  n_x_ = x_.size();

  n_aug_ = n_x_ + 2;
  
  n_sig_ = 2 * n_aug_ + 1;

  lambda_ = 3 - n_aug_;

  time_us_ = 0;

  weights_ = VectorXd(n_sig_);

  Xsig_pred = MatrixXd(n_x_, n_sig_);
}

UKF::~UKF() {}

double NormalizeAngle(double phi)
{
  return  atan2(sin(phi), cos(phi));
}
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
   if(!is_init){

    x_ << 0, 0, 0, 0, 0;

    // x_aug.head(5) = x_;
    // x_aug(5) = 0;
    // x_aug(6) = 0; 


    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    // P_aug.fill(0.0);
    // P_aug.topLeftCorner(n_x_, n_x_) = P_;
    // P_aug(5,5) = std_a_ * std_a_;
    // P_aug(6,6) = std_yawdd_ * std_yawdd_;



    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
        float rho = meas_package.raw_measurements_[0]; // range
        float phi = meas_package.raw_measurements_[1]; // bearing
        float rho_dot = meas_package.raw_measurements_[2]; // velocity of rho

        if (fabs(rho) < 0.001)
        {
          rho = 0.001;
        }
        // Coordinates convertion from polar to cartesian
        float px = rho * cos(phi); 
        float py = rho * sin(phi);
        float vx = rho_dot * cos(phi);
        float vy = rho_dot * sin(phi);
        float v  = sqrt(vx * vx + vy * vy);
        x_ << px, py, v, 0, 0;


      }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        /**
        Initialize state.
        */
        float px = meas_package.raw_measurements_[0];
        float py = meas_package.raw_measurements_[1];
        x_ << px, py, 0, 0, 0;
      
        if(fabs(x_(0)) < 0.001 and fabs(x_(1) < 0.001)){
          x_(0) = 0.001;
          x_(1) = 0.001;
          }
      }
  
    weights_.fill(0.5/(n_aug_ + lambda_));
    weights_(0) = lambda_/(lambda_ + n_aug_);


    time_us_ = meas_package.timestamp_;

    is_init = true;
    cout << "finish init"<<endl;
    cout << "x_" << endl << x_ << endl;
    cout << "P_" << endl << P_ << endl;
    return;
 }

  // prediction 
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  // cout << "predict:" << endl;
  // cout << "x_" << endl << x_ << endl;
  // cout << "P_" << endl << P_ << endl;
  // cout << "sig_pred" <<Xsig_pred << endl;
  // update 
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // cout << "Use Radar" << endl;
    UpdateRadar(meas_package);
    // cout << "Radar Update Complete " << endl;
    }

  if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    // cout << "Use laser" << endl;
    UpdateLidar(meas_package);
    // cout << "lidar Update Complete " << endl;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {float} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(float delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //predict sigma points
  
  float delta_t2 = delta_t * delta_t;

  VectorXd x_aug = VectorXd(n_aug_);

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  //set x_aug using x_ 
  //x_aug =[px, py, v(speed m/s), angle(rad), angle(rad/s) noise, noise ]
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    float p_x = Xsig_aug(0,i);
    float p_y = Xsig_aug(1,i);
    float v = Xsig_aug(2,i);
    float yaw = Xsig_aug(3,i);
    float yawd = Xsig_aug(4,i);
    float nu_a = Xsig_aug(5,i);
    float nu_yawdd = Xsig_aug(6,i);

    //avoid division by zero

    if (fabs(yawd) > 0.001){
      // prediction + noise
      Xsig_pred(0,i) = p_x + (v / yawd) * (sin(yaw + yawd*delta_t) - sin(yaw)) + 0.5 * nu_a * delta_t2 * cos(yaw);
      Xsig_pred(1,i) = p_y + (v / yawd) * (cos(yaw) - cos(yaw + yawd * delta_t)) + 0.5 * nu_a * delta_t2 * sin(yaw);
      Xsig_pred(2,i) = v + nu_a * delta_t;
      Xsig_pred(3,i) = yaw + yawd * delta_t + 0.5 * nu_yawdd * delta_t2;
      Xsig_pred(4,i) = yawd + nu_yawdd * delta_t;
    }
    else {
      Xsig_pred(0,i) = p_x + v * delta_t * cos(yaw) + 0.5*nu_a*delta_t2 * cos(yaw);
      Xsig_pred(1,i) = p_y + v * delta_t * sin(yaw) + 0.5*nu_a*delta_t2 * sin(yaw);
      Xsig_pred(2,i) = v + nu_a * delta_t;
      Xsig_pred(3,i) = yaw + 0.5 * nu_yawdd * delta_t2;
      Xsig_pred(4,i) = yawd + nu_yawdd * delta_t;
    }
  }
  //predicted state mean
  x_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    //angle normalization

    x_diff(3) = NormalizeAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;

  }
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
  int n_z = 2; 

  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  Zsig.fill(0.0);

  VectorXd z = meas_package.raw_measurements_;

  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  MatrixXd R = MatrixXd(n_z,n_z);

   //transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    // extract values for better readibility
    float p_x = Xsig_pred(0,i);
    float p_y = Xsig_pred(1,i);
    

    // measurement model
    Zsig(0,i) = p_x;                       
    Zsig(1,i) = p_y;                                 
  
  }

  //mean predicted measurement

  for (int i=0; i < n_sig_; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S

  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization

    z_diff(1) = NormalizeAngle(z_diff(1));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix

  R <<    std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;

  S = S + R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;


    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;


    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;



  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  // Calculate and update NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  cout << "NIS_laser_ = " << NIS_laser_ << endl;
  // cout << "x_ after Lidar update = " << endl << x_ << endl;
  // cout << "P_ after Lidar update = " << endl << P_ << endl;
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
  int n_z = 3; //rho, phi, rho_dot 

  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  Zsig.fill(0.0);

  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0.0);
  VectorXd z = meas_package.raw_measurements_;
   //transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    // extract values for better readibility
    const double p_x = Xsig_pred(0,i);
    const double p_y = Xsig_pred(1,i);
    const double v  = Xsig_pred(2,i);
    const double yaw = Xsig_pred(3,i);

    const double v1 = cos(yaw)*v;
    const double v2 = sin(yaw)*v;

    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //rho
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //rho_dot
  }
  //   if (fabs(p_x) < 0.001){
  //     p_x = 0.001;
  //   }

  //   if (fabs(p_y) < 0.001){
  //     p_y = 0.001;
  //   }
    
  //   float rho = sqrt(p_x*p_x + p_y*p_y);

  //   if (fabs(rho) < 0.001)
  //   {
  //     rho = 0.001;
  //   }
  //   // measurement model
  //   Zsig(0,i) = rho;                        //rho
  //   Zsig(1,i) = atan2(p_y,p_x);             //phi
  //   Zsig(2,i) = (p_x*v1 + p_y*v2 ) / rho;   //rho_dot
  // }

  //   if (p_x == 0 && p_y == 0)
  //   {
  //     Zsig(0,i) = 0;
  //     Zsig(1,i) = 0;
  //     Zsig(2,i) = 0;
  //   } else {
  //     Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //rho
  //     Zsig(1,i) = atan2(p_y,p_x);                                 //phi
  //     Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //rho_dot
  //   }
  // }

  //mean predicted measurement

  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = NormalizeAngle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix

  R <<    std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0,std_radrd_ * std_radrd_;
  S = S + R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    //angle normalization
    x_diff(3) = NormalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = NormalizeAngle(z_diff(1));
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  // Calculate and update NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  // cout << "NIS_radar_ = " << NIS_radar_ << endl;
  // cout << "x_ after radar update = " << endl << x_ << endl;
  // cout << "P_ after radar update = " << endl << P_ << endl;
}
