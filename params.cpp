#include<iostream>
#include<cmath>
#include<Eigen/Dense>
#include<fstream>
#include<omp.h>


//using Eigen::MatrixXd;
//typedef Eigen::MatrixXd dmatrix;
using dmatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;   //the loops are accessed in row major so this is a dirty way to improve the performance without upending the code


//grid:
const double lx=1.0;
const double ly=1.0;
const int nx=128;
const int ny=128;


const double dx=lx/nx;
const double dy=ly/ny;

//pseudo speed of sound parameter 
const double delta=4.5;

//simulation parameters:
const double dt=0.0001;
const double Re=1000;
const double tol=0.0000001;

//data parameter:
const int print_interval=100000;
