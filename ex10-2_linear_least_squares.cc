// Linear Least Squares parameter estimation. Example 10-2 of Vallado 4ed.
// Note: the determinant, inverse, and dot product functions are written
// specifically for this problem and can only take 2x2 and 1x2 vectors as input. 

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <tuple>

double determinant(std::vector<std::vector<double> > arr);
std::vector<std::vector<double> > inverse(std::vector< std::vector<double> > arr);
std::vector<double> dot_product(std::vector<std::vector<double> > a, std::vector<double> b); 
std::vector<double> propagate(std::vector<double> state_, std::vector<double> xobs_);
std::vector<double> residuals(std::vector<double> yobs_, std::vector<double> yexp_);
double root_mean_sq(std::vector<double> res_);
std::tuple< std::vector<std::vector<double> >, std::vector<double> > 
                        fill_ata_atb(std::vector<double> xobs_, std::vector<double> yobs_); 

int main() {

  int i, ii;

  // The state vector (alpha and beta)
  std::vector<double> state;
  // The inverse of the information matrix (A^T * A)^-1
  // i.e. the covariance matrix without weighting
  std::vector<std::vector<double> > inv;

  // The observed values of x  
  std::vector<double> xobs = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  // The observed values of y (the observation matrix, b)
  std::vector<double> yobs = {1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 7.0, 6.0};

  // The two components of the general solution (A^T * A, and A^T * b)
  std::vector<std::vector<double> > ata;
  std::vector<double> atb;

  // The expected values after propagation
  std::vector<double> yexp;
  // Residuals (observed - expected). 
  // Must be initialised to <= RMS for first loop 
  std::vector<double> res(xobs.size(), 0.0);
  // Absolute value of residuals. 
  std::vector<double> absres(xobs.size(), 0.0);
  // Root mean square of residuals.
  // Must be initialised  >= residuals for first loop
  double rms = 0.0;
  // Observational variance
  double obs_var;
  // Standard deviation in alpha, beta
  double stdev_a;
  double stdev_b;
 
  // Complete two iterations of the propagate and fit algorithm
  for ( i = 0; i < 2; i++ ) {
    std::cout << "\n********* Loop #" << i + 1 << " *********" <<std::endl;
    // Remove observed value with residuals outside 2*RMS
    // This will have no affect during the first iteration, if initialised 
    // RMS >= residuals.
    for ( ii = 0; ii < xobs.size(); ii++ ) {
      if (std::abs(res[ii]) > 2.0 * rms) {
        std::cout << "Dropping observation " << ii << ": (" << xobs[ii] << "," << yobs[ii] << ")" << std::endl;
        xobs.erase(xobs.begin() + ii);
        yobs.erase(yobs.begin() + ii);
        res.erase(res.begin() + ii);
        ii--;
      }
    }

    // Fill the matrices for the general solution
    tie(ata, atb) = fill_ata_atb(xobs, yobs);
    // Calculate the covariance matrix (A^T * A)*-1
    inv = inverse(ata);
    // Calculate the state vector
    state = dot_product(inv, atb);

    std::cout << "ATA matrix  " << ata[0][0] << "  " << ata[0][1] << std::endl;
    std::cout << "            " << ata[1][0] << "  " << ata[1][1] << std::endl;
    std::cout << "ATb matrix  " << atb[0] << "  " << atb[1] << std::endl; 
  
    // Propagate the x observations with the state to find the expected y values
    yexp = propagate(state, xobs);

    std::cout << "\nCalculated Y values" << std::endl;
    for ( ii = 0; ii < yexp.size(); ii++ ) {
      std::cout << yexp[ii] << " ";
    }
    std::cout << std::endl;

    // Calculate the residuals
    res = residuals(yobs, yexp);
    std::transform(res.begin(), res.end(), absres.begin(),
                         static_cast<float (*)(float)>(&std::abs));

    std::cout << "\nResidual values" << std::endl;
    for ( ii = 0; ii < res.size(); ii++ ) {
      std::cout << absres[ii] << " ";
    }
    std::cout << std::endl;

    // Calculate the RMS of the residuals
    rms = root_mean_sq(res);

    std::cout << "\nRMS: " << rms << std::endl;

    // Standard deviation in alpha and beta are taken from the covariance matrix
    stdev_a = std::sqrt(inv[0][0]);
    stdev_b = std::sqrt(inv[1][1]);

    // Without a weighting matrix, observational variance is estimated as:
    obs_var = std::sqrt(pow(rms, 2.0) * res.size() / (res.size() - 1.0)); 

    std::cout << "\nState:" << std::endl;
    std::cout << "alpha " << state[0] << " +/- " << obs_var * stdev_a << std::endl;
    std::cout << "beta  " << state[1] << " +/- " << obs_var * stdev_b << std::endl;
 
  }
}

std::tuple< std::vector<std::vector<double> >, std::vector<double> > 
            fill_ata_atb(std::vector<double> xobs_, std::vector<double> yobs_) {
  int i;
  // Fill the matrix components of the general solution, as shown in Vallado
  // Exercise 10-2.

  // The two components of the general solution (A^T * A, and A^T * B)
  std::vector<std::vector<double> > ata_ = {{0.0, 0.0},{0.0, 0.0}};
  std::vector<double> atb_ = {0.0, 0.0};
  // fill ata and atb
  for ( i = 0; i < xobs_.size(); i++ ) {
    ata_[0][0] += 1.0;
    ata_[0][1] += xobs_[i];
    ata_[1][0] += xobs_[i];
    ata_[1][1] += pow(xobs_[i], 2);

    atb_[0] += yobs_[i];
    atb_[1] += xobs_[i] * yobs_[i];
  }
  return std::make_tuple(ata_, atb_);
} 

double determinant(std::vector<std::vector<double> > arr) {
  // Calculate the determinant of a 2x2 matrix
  double det;

  det = (arr[0][0] * arr[1][1]) - (arr[0][1] * arr[1][0]);

  return det;
}

std::vector<std::vector<double> > inverse(std::vector<std::vector<double> > arr) {
  // Calculate the inverse of a 2x2 matrix

  std::vector<std::vector<double> > inv(arr[0].size(), std::vector<double>(arr[1].size()));
  
  double det = determinant(arr);

  inv[0][0] = arr[1][1]  / det;
  inv[1][1] = arr[0][0]  / det;
  inv[0][1] = -arr[1][0] / det;
  inv[1][0] = -arr[0][1] / det;

  return inv;
}

std::vector<double> dot_product(std::vector<std::vector<double> > a, 
                                  std::vector<double> b) {
  // calculate the dot production of a 2x2 matrix with a 1x2 vector
                          
  std::vector<double> result(2, 0.0);
  std::vector<double> a0 = {a[0][0], a[0][1]};
  std::vector<double> a1 = {a[1][0], a[1][1]};

  result[0] = std::inner_product(a0.begin(),a0.end(),b.begin(),0.0);
  result[1] = std::inner_product(a1.begin(),a1.end(),b.begin(),0.0);

  return result;
}

std::vector<double> propagate(std::vector<double> state_, 
                                   std::vector<double> xobs_) {
  // Propagate the observed x value using the state vector
  // to calculate the expected y values.

  int i;
  std::vector<double> yexp_(xobs_.size(), 0.0);

  for (i = 0; i < xobs_.size(); i++) {
    yexp_[i] = state_[0] + state_[1] * xobs_[i];   
  }

  return yexp_;
}

std::vector<double> residuals(std::vector<double> yobs_, 
                                 std::vector<double> yexp_) {
  // calculate the residuals (observed  - expected)

  // residuals vector
  std::vector<double> res_(yobs_.size(), 0.0);

  std::transform(yobs_.cbegin(), yobs_.cend(), 
                   yexp_.cbegin(), res_.begin(), std::minus<double>());

  return res_;  
}

double root_mean_sq(std::vector<double> res_) {
  // Calculate the square root of the mean expected squared value
  double rms_;
  std::vector<double> ressq(res_.size(), 0.0);
  std::transform(res_.cbegin(), res_.cend(),
                   res_.cbegin(), ressq.begin(), std::multiplies<double>());
  rms_ = std::accumulate(ressq.begin(), ressq.end(), 0.0) / res_.size();
  rms_ = std::sqrt(rms_);

  return rms_;
}
