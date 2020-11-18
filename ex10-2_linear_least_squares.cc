// Linear Least Squares parameter estimation. Example 10-2 of Vallado 4ed.
// Note: the determinant, inverse, and dot product functions can only take
// 2x2 and 1x2 vectors as input. For parameter estimation with two variables
// vectors of these dimensions are guaranteed.

// TO DO: 
// Makefile
// gitignore
// read/write : https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>

double determinant(std::vector<std::vector<double>> arr);
std::vector<std::vector<double>> inverse(std::vector< std::vector<double>> arr);
std::vector<double> dot_product(std::vector<std::vector<double>> a, std::vector<double> b); 
std::vector<double> propagate(std::vector<double> state_, std::vector<double> xobs_);
std::vector<double> residuals(std::vector<double> yobs_, std::vector<double> yexp_);

int main() {

  int i, j;

  // The state vector (alpha and beta)
  std::vector<double> state;
  // The inverse of the information matrix (A^T * A)^-1
  std::vector<std::vector<double>> inv;

  // The observed values of x and y
  std::vector<double> xobs = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  std::vector<double> yobs = {1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 7.0, 6.0};

  // The two components of the general solution (A^T * A, and A^T * B)
  std::vector<std::vector<double>> ata = {{0.0, 0.0},{0.0, 0.0}};
  std::vector<double> atb = {0.0, 0.0};

  // The expected values after propagation
  std::vector<double> yexp;

  // fill ata and atb
  for ( i = 0; i < 8; i++ ) {
    ata[0][0] += 1.0;
    ata[0][1] += xobs[i];
    ata[1][0] += xobs[i];
    ata[1][1] += pow(xobs[i], 2);

    atb[0] += yobs[i];
    atb[1] += xobs[i] * yobs[i];
  }

  inv = inverse(ata);
  state = dot_product(inv, atb);

  std::cout << "ATA matrix  " << ata[0][0] << "  " << ata[0][1] << std::endl;
  std::cout << "            " << ata[1][0] << "  " << ata[1][1] << std::endl;
  std::cout << "ATb matrix  " << atb[0] << "  " << atb[1] << std::endl; 
  
  std::cout << "\nState:" << std::endl;
  std::cout << "alpha " << state[0] << std::endl;
  std::cout << "beta  " << state[1] << std::endl;

  //yexp = propagate(state, xobs);

  //for ( i = 0; i < 8; i++ ) {
  //  std::cout << yexp[i] << std::endl;
  //}

}

double determinant(std::vector<std::vector<double> > arr) {
  // Calculate the determinant of a 2x2 matrix
  double det;

  det = (arr[0][0] * arr[1][1]) - (arr[0][1] * arr[1][0]);

  return det;
}

std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> arr) {
  // Calculate the inverse of a 2x2 matrix

  int i, j;
  // tmp variable to store matrix elements
  double tmp;  
  std::vector<std::vector<double>> inv = arr;
  
  double det = determinant(arr);

  inv[0][0] = arr[1][1] / det;
  inv[1][1] = arr[0][0] / det;
  inv[0][1] = -arr[1][0] / det;
  inv[1][0] = -arr[0][1]/ det;

  return inv;
}

std::vector<double> dot_product(std::vector<std::vector<double>> a, 
                                  std::vector<double> b) {
  // calculate the dot production of a 2x2 matrix with a 1x2 vector
                          
  std::vector<double> result = {0.0, 0.0};
  std::vector<double> a0 = {a[0][0], a[0][1]};
  std::vector<double> a1 = {a[1][0], a[1][1]};

  result[0] = std::inner_product(a0.begin(),a0.end(),b.begin(),0.0);
  result[1] = std::inner_product(a1.begin(),a1.end(),b.begin(),0.0);

  return result;
}

std::vector<double> propagate(std::vector<double> state_, 
                                   std::vector<double> xobs_) {

  int i;
  std::vector<double> yexp_;

  for (i = 0; i < xobs_.size(); i++) {
    yexp_[i] = state_[0] + state_[1] * xobs_[i];   
  }

  return yexp_;
}

//std::vector<double> residuals(std::vector<double> obs, 
//                                 std::vector<double> mdl) {
//  // calculate the residuals (observed  - modelled)
//
//  // residuals vector
//  std::vector<double> res;
//
//  
//}