#ifndef _FEM_FLUID_H_
#define _FEM_FLUID_H_


#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fluid_domain.h"
#include "allocation.h"
#include "Eigen/Core"
#include "Eigen/LU"

enum class DataType{
 VELOCITY=0,
 PRESSURE=1,
 PHI=2
};

using namespace Eigen;
using namespace std;

typedef Matrix<double, 4, 2> Matrix4_2d;

class FEMF :public fluidDomain{
  private:
  void inputDomainInfo();
  void inputBoundaryInfo();
  void inputImageInfo();
  void inputSolverInfo();
  void inputOutputInfo();
  void allocate();

  public:
  /*
  VectorXd un;
  VectorXd vn;
  Vector4d umid;
  Vector4d vmid;
  Vector4d p;
  VectorXd R;
  MatrixXd K;
  VectorXd R2;
  MatrixXd K2;
  VectorXd R3;
  MatrixXd K3;
  Vector4d N;
  Matrix4_2d dNdr;
  Matrix4_2d x_current;
  VectorXd U2;
  VectorXd U3;
  */
  double Renum;
  double nu;
  double rho;
  double mu;
  double dt;
  int x_size;
  int boundary_element_size;
  int s;
  int check;
  int t;
  int time;
  string outputDir,fileName;

  void export_vtu(const std::string &file,DataType dataType,VectorXd &U2,VectorXd &U3);
  void initialize();
  void main_NavierStokes();
  void umidMatrix(MatrixXd &un,MatrixXd &vn,MatrixXd &K,VectorXd &R,Vector4d &N,Matrix4_2d &dNdr,Matrix4_2d &x_current,const int numOfNodeInElm,const double weight,const int ic,const int t);
  void bc(MatrixXd &K,VectorXd &R,int x,double u_value,double v_value,int check);
  int solve_Matrix(MatrixXd &K,VectorXd &R,VectorXd &a,int s);
  void poissonMatrix(Vector4d &umid,Vector4d &vmid,MatrixXd &K2,VectorXd &R2,Vector4d &N,Matrix4_2d &dNdr,Matrix4_2d &x_current,const int numOfNodeInElm,const double weight,const int ic);
  void u_nextMatrix(Vector4d &umid,Vector4d &vmid,Vector4d &p,MatrixXd &K3,VectorXd &R3,Vector4d &N,Matrix4_2d &dNdr,Matrix4_2d &x_current,const int numOfNodeInElm,const double weight,const int ic);
};


#endif