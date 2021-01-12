#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

// C libraries
#include <float.h>		// DBL_MAX
#include <math.h> 		// sin, cos
#include <stdio.h>      // printf, fgets
#include <stdlib.h>     // atoi
#include <omp.h>		// OpenMP
#include <time.h>		// time
#include <mpi.h>		// MPI

// C++ libraries
#include <string>		// strings
using std::string;

#include <iostream>		// terminal output
#include <cmath>
using std::cout;
using std::endl;

#include <fstream>		// file io
using std::ifstream;
using std::getline;
using std::istringstream;

#include <algorithm>	// math stuff
using std::min;
using std::max;
using std::find;

#include <chrono>
#include <random>

#include <vector>  // standard vector
using std::vector;
 
#include <ctime>
using std::clock_t;

// #include<tr1/array>
// using std::tr1::array;

#include <iomanip>

#include </home/siavash/Downloads/eigen-eigen-323c052e1731/Eigen/Dense>
#include </home/siavash/Downloads/eigen-eigen-323c052e1731/Eigen/Core>
#include </home/siavash/Downloads/eigen-eigen-323c052e1731/Eigen/IterativeLinearSolvers>
#include </home/siavash/Downloads/eigen-eigen-323c052e1731/Eigen/Sparse>

using Eigen::Vector2d;
using Eigen::Vector2i;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix2d;
using Eigen::Matrix;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

typedef Matrix<double, 8, 1> Vector8d;
typedef Matrix<double, 6, 1> Vector6d;

#include </home/siavash/Downloads/eigen-eigen-323c052e1731/Eigen/StdVector>




#endif /* DEFINITIONS_H_ */
