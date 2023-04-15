#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/SparseLU"
#include "fem_fluid.h"



using namespace std;
using namespace Eigen;


typedef Matrix<double, 2, 2> Matrix2_2d;
typedef Matrix<double, 4, 4> Matrix4_4d;
typedef Matrix<double, 4, 2> Matrix4_2d;
typedef Matrix<double, 2, 4> Matrix2_4d;
typedef Matrix<double, 8, 8> Matrix8_8d;
typedef Matrix<double, 6, 2> Matrix6_2d;
typedef Matrix<double, 2, 6> Matrix2_6d;
typedef Matrix<double, 8, 1> Vector8d;
typedef Triplet<double> T;

void C2D4_dNdr(Matrix4_2d &dNdr,const double &g1,const double &g2);
void C2D4_N(Vector4d &N,const double &g1,const double &g2);

void calc_dxdr2D(Matrix2_2d &dxdr,Matrix4_2d &dNdr,Matrix4_2d &x1,const int &numOfNodeInElm);
void calc_dNdx2D(Matrix4_2d &dNdx,Matrix4_2d &dNdr,Matrix2_2d &dxdr ,const int &numOfNodeInElm);
void umidMatrix(VectorXd &un,VectorXd &vn,MatrixXd &K,VectorXd &R,Vector4d &N,Matrix4_2d &dNdr,Matrix4_2d &x_current,const int numOfNodeInElm,const double weight,const int ic);
void bc(MatrixXd &K,VectorXd &R,int x,double u_value,double v_value,int check);
int solve_Matrix(MatrixXd &K,VectorXd &R,VectorXd &a,int s);
void poissonMatrix(VectorXd &umid,VectorXd &vmid,MatrixXd &K2,VectorXd &R2,Vector4d &N,Matrix4_2d &dNdr,Matrix4_2d &x_current,const int numOfNodeInElm,const double weight,const int ic);
void u_nextMatrix(VectorXd &umid,VectorXd &vmid,Vector4d &p,MatrixXd &K3,VectorXd &R3,Vector4d &N,Matrix4_2d &dNdr,Matrix4_2d &x_current,const int numOfNodeInElm,const double weight,const int ic);

class Gauss{
public:
	Gauss(){}
	Gauss(const int num){
		switch(num){
			case 1:
	  		point[0] = -0.577350296189626; point[1] =  0.577350296189626;
      	weight[0] = 1e0;							 weight[1] = 1e0;
      	break;
      case 2:
      	point[0] = -0.774596669241483; point[1] = 0e0; point[2] = 0.774596669241483;
      	weight[0] = 0.555555555555555; weight[1] = 0.888888888888888; weight[2] = 0.555555555555555;
      	break;
      case 3:
      	point[0] = -0.861135311594053;  point[1] = -0.339981043584856; point[2] = 0.339981043584856;  point[3] = 0.861135311594053;
      	weight[0] = 0.347854845137454;  weight[1] = 0.652145154862546; weight[2] = 0.652145154862546; weight[3] = 0.347854845137454;
      break;
      default:
      	printf("undefined order is set in the gauss integral\n");
		}
	}
  double point[4],weight[4];
};


int main(){
  FEMF FemF;
  FemF.initialize();
  FemF.main_NavierStokes();

}

void FEMF::main_NavierStokes()
{

  int l = 0;
  int m = 0;
  int n = 0;


  Matrix4_2d dNdr;
  Matrix4_2d x_current;
  Matrix4_2d dNdx;
  Matrix2_2d dxdr;
  Vector4d N;
  MatrixXd u(time,x_size);
  MatrixXd v(time,x_size);
  Vector4d umid;
  Vector4d vmid;
  Vector4d p;
  VectorXd U(x_size*2);
  VectorXd R(x_size*2);
  MatrixXd K(x_size*2, x_size*2);
  VectorXd U2(x_size);
  VectorXd R2(x_size);
  MatrixXd K2(x_size, x_size);
  VectorXd U3(x_size*2);
  VectorXd R3(x_size*2);
  MatrixXd K3(x_size*2, x_size*2);




  for(int t = 0; t < time; t++) {
    cout << "t=" << t << endl;

    for(int i=0; i<4; i++){
      umid(i) = 0e0;
      vmid(i) = 0e0;
      p(i) = 0e0;
    }

    for(int i=0; i<x_size; i++){
      U2(i) = 0e0;
      R2(i) = 0e0;
      u(t,i) = 0e0;
      v(t,i) = 0e0;
    }

    for(int i=0; i<x_size*2; i++){
      U(i) = 0e0;
      R(i) = 0e0;
    }

    for(int i=0; i<x_size; i++){
      u(t,i) = U3(i);
      v(t,i) = U3(i+x_size);
    }

    for(int i=0; i<x_size*2; i++){
      U3(i) = 0e0;
      R3(i) = 0e0;
    }

    for(int i=0; i<x_size*2; i++){
      for(int j=0; j<x_size*2; j++){
        K(i,j) = 0e0;
      }
    }

    for(int i=0; i<x_size; i++){
      for(int j=0; j<x_size; j++){
        K2(i,j) = 0e0;
      }
    }

    for(int i=0; i<x_size*2; i++){
      for(int j=0; j<x_size*2; j++){
        K3(i,j) = 0e0;
      }
    }



    
    for(int i=boundary_element_size; i<boundary_element_size*2; i++){
      u(t,boundaryElement[i].node[0]) = 1e0;
      
        if(i == boundary_element_size*2-1){
          u(t,boundaryElement[i].node[1]) = 1e0;
        }
    }

    /*
    for(int i=boundary_element_size*2; i<boundary_element_size*3; i++){
      u(t,boundaryElement[i].node[0]) = 1e0;
    }

    u(t,boundaryElement[0].node[0]) = 1e0;
    */

    
    for(int ic = 0; ic<numOfElm; ic++){


      int numOfNodeInElm = element[ic].node.size();
    
      for(int i=0;i<numOfNodeInElm;i++){
        for(int j=0;j<2;j++){
          x_current(i,j) = x(element[ic].node[i],j);
        }
      }
      Gauss g1(1),g2(2);
      

      for(int i1=0;i1<3;i1++){
        for(int i2=0;i2<3;i2++){
          C2D4_N(N,g2.point[i1],g2.point[i2]);
          C2D4_dNdr(dNdr,g2.point[i1],g2.point[i2]);
          umidMatrix(u,v,K,R,N,dNdr,x_current,numOfNodeInElm,g2.weight[i1]*g2.weight[i2],ic,t);
        }
      }
    }
  
    bc(K,R,boundary_element_size*2,0e0,0e0,1);//left
    bc(K,R,0,0e0,0e0,1);//bottom
    bc(K,R,boundary_element_size*3,0e0,0e0,1);//right
    bc(K,R,boundary_element_size,1e0,0e0,0);//top

    
    solve_Matrix(K,R,U,x_size*2);
    

    for(int ic = 0; ic<numOfElm; ic++){

      int numOfNodeInElm = element[ic].node.size();

      for(int i=0; i<numOfNodeInElm; i++){
        umid(i) = U(element[ic].node[i]);
        vmid(i) = U(x_size+element[ic].node[i]);
      }
    
      for(int i=0;i<numOfNodeInElm;i++){
        for(int j=0;j<2;j++){
          x_current(i,j) = x(element[ic].node[i],j);
        }
      }
      Gauss g1(1),g2(2);
      

      for(int i1=0;i1<3;i1++){
        for(int i2=0;i2<3;i2++){
          C2D4_N(N,g2.point[i1],g2.point[i2]);
          C2D4_dNdr(dNdr,g2.point[i1],g2.point[i2]);
          poissonMatrix(umid,vmid,K2,R2,N,dNdr,x_current,numOfNodeInElm,g2.weight[i1]*g2.weight[i2],ic);
        }
      }
    }
    
    
    for (int j = 0; j < x_size; j++)
    {
      K2(0, j) = 0e0;
    }
    K2(0, 0) = 1e0;
    R2(0) = 0e0;
    
    
   /*
    for (int i = boundary_element_size*3; i < boundary_element_size*4; i++)
    {
      for (int j = 0; j < x_size ; j++)
      {
        K2(boundaryElement[i].node[0], j) = 0e0;
        if(i==boundary_element_size*4-1){
           K2(boundaryElement[i].node[1], j) = 0e0;
        }
      }
      K2(boundaryElement[i].node[0],boundaryElement[i].node[0]) = 1e0;
      if(i==boundary_element_size*4-1){
        K2(boundaryElement[i].node[1],boundaryElement[i].node[1]) = 1e0;
      }
      R2(boundaryElement[i].node[0]) = 0e0;
      if(i==boundary_element_size*4-1){
        R2(boundaryElement[i].node[1]) = 0e0;
      }
    }*/


    

    solve_Matrix(K2,R2,U2,x_size);


    for(int ic = 0; ic<numOfElm; ic++){

      int numOfNodeInElm = element[ic].node.size();

      for(int i=0; i<numOfNodeInElm; i++){
        umid(i) = U(element[ic].node[i]);
        vmid(i) = U(x_size+element[ic].node[i]);
      }

      for(int i=0; i<numOfNodeInElm; i++){
        p(i) = U2(element[ic].node[i]);
      }
    
      for(int i=0;i<numOfNodeInElm;i++){
        for(int j=0;j<2;j++){
          x_current(i,j) = x(element[ic].node[i],j);
        }
      }
      Gauss g1(1),g2(2);

      for(int i1=0;i1<3;i1++){
        for(int i2=0;i2<3;i2++){
          C2D4_N(N,g2.point[i1],g2.point[i2]);
          C2D4_dNdr(dNdr,g2.point[i1],g2.point[i2]);
          u_nextMatrix(umid,vmid,p,K3,R3,N,dNdr,x_current,numOfNodeInElm,g2.weight[i1]*g2.weight[i2],ic);
        }
      }
    }


    bc(K3,R3,boundary_element_size*2,0e0,0e0,1);//left
    bc(K3,R3,0e0,0e0,0e0,1);//bottom
    bc(K3,R3,boundary_element_size*3,0e0,0e0,1);//right
    bc(K3,R3,boundary_element_size,1e0,0e0,0);//top
    
 
   
    solve_Matrix(K3,R3,U3,x_size*2);
 
    int w=0;
   
    ofstream fv_mid, fv, fp;
    fv_mid.open("Data/vel_mid" + to_string(l) + ".csv");
    fv.open("Data/vel" + to_string(m) + ".csv");
    fp.open("Data/pre" + to_string(n) + ".csv");
    l++;
    m++;
    n++;

    fv_mid << "x" << " " << "y" << " " << "z" << " " << "u" << " " << "v" << " " << "w" << endl;
    for(int i=0; i<x_size; i++){
      fv_mid << x(i,0) << " " << x(i,1) << " " << 0 << " " << U(i) << " " << U(x_size+i) << " " << w << endl;
    }
    
    fv << "x" << " " << "y" << " " << "z" << " " << "u" << " " << "v" << " " << "w" << endl;
    for(int i=0; i<x_size; i++){
      fv << x(i,0) << " " << x(i,1) << " " << 0 << " " << U3(i) << " " << U3(x_size+i) << " " << w << endl;
    }

    fp << "x" << " " << "y" << " " << "z" << " " << "p"  << endl;
    for(int i=0; i<x_size; i++){
      fp << x(i,0) << " " << x(i,1) << " " << 0 << " " << U2(i)  << endl;
    }
    string vtuFile;
   
    vtuFile = "Result/result_"+to_string(t)+".vtu";
    export_vtu(vtuFile,DataType::VELOCITY,U2,U3);
    vtuFile = "Result2/result_pressure_"+to_string(t)+".vtu";
    export_vtu(vtuFile,DataType::PRESSURE,U2,U3);
    

  //exit(1);
  }

  ofstream outputfile("x_center.txt");
  outputfile << "y" << endl;
  for(int i=0; i<ny+1; i++){
    for(int j=0; j<nx+1; j++){
      if(j == nx/2){
        outputfile << x(i*(nx+1)+j,1) << endl;
      }
    }
  }
  outputfile.close();


 ofstream outputfile1("x_center_u.txt");
  outputfile1 << "u" << endl;
  for(int i=0; i<ny+1; i++){
    for(int j=0; j<nx+1; j++){
      if(j == nx/2){
        outputfile1 <<  U(i*(nx+1)+j) << endl;
      }
    }
  }
  outputfile1.close();
}

void C2D4_N(Vector4d &N,const double &g1,const double &g2)
  {
    N(0) = 2.5e-1 * (1e+0-g1) * (1e+0-g2);
    N(1) = 2.5e-1 * (1e+0+g1) * (1e+0-g2);
    N(2) = 2.5e-1 * (1e+0+g1) * (1e+0+g2);
    N(3) = 2.5e-1 * (1e+0-g1) * (1e+0+g2);
  }

void C2D4_dNdr(Matrix4_2d &dNdr,const double &g1,const double &g2)
  {
    dNdr(0,0) = -2.5e-1 * (1e+0-g2);
    dNdr(0,1) = -2.5e-1 * (1e+0-g1);
    dNdr(1,0) =  2.5e-1 * (1e+0-g2);
    dNdr(1,1) = -2.5e-1 * (1e+0+g1);
    dNdr(2,0) =  2.5e-1 * (1e+0+g2);
    dNdr(2,1) =  2.5e-1 * (1e+0+g1);
    dNdr(3,0) = -2.5e-1 * (1e+0+g2);
    dNdr(3,1) =  2.5e-1 * (1e+0-g1);
  }


void calc_dxdr2D(Matrix2_2d &dxdr,Matrix4_2d &dNdr,Matrix4_2d &x1,const int &numOfNodeInElm)
{
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      dxdr(i,j) = 0e0;
      for(int p=0;p<numOfNodeInElm;p++){
        dxdr(i,j) += dNdr(p,j) * x1(p,i);
      }
    }
  }
}

void calc_dNdx2D(Matrix4_2d &dNdx,Matrix4_2d &dNdr,Matrix2_2d &dxdr,const int &numOfNodeInElm)
{
  Matrix2_2d drdx;
  drdx = dxdr.inverse();

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<2;i++){
      dNdx(p,i) = 0e0;
      for(int j=0;j<2;j++) dNdx(p,i) += dNdr(p,j) * drdx(j,i);
    }
  }
}

void FEMF::umidMatrix(MatrixXd &un,MatrixXd &vn,MatrixXd &K,VectorXd &R,Vector4d &N,Matrix4_2d &dNdr,Matrix4_2d &x_current,const int numOfNodeInElm,const double weight,const int ic, const int t)
{
  Matrix4_2d dNdx;
  Matrix2_2d dxdr;
  Matrix8_8d M;
  Matrix8_8d Ms;
  Matrix8_8d Ke;
  Vector8d Re;
  Vector8d A;
  Vector8d D;
  Vector8d As;
  Vector8d Mr;
  Vector8d Mrs;


  calc_dxdr2D(dxdr,dNdr,x_current,numOfNodeInElm);
  double detJ = dxdr.determinant();
  calc_dNdx2D(dNdx,dNdr,dxdr,numOfNodeInElm);


  
  double vel[2]={0e0,0e0};
  for(int p=0;p<numOfNodeInElm;p++){
    vel[0] += N(p)*un(t,element[ic].node[p]);
    vel[1] += N(p)*vn(t,element[ic].node[p]);
  }

  double pre_vel[2]={0e0,0e0};
  if(t==0){
    for(int p=0;p<numOfNodeInElm;p++){
      pre_vel[0] = 0e0;
      pre_vel[1] = 0e0;
    }
  }
  if(t!=0){
    for(int p=0;p<numOfNodeInElm;p++){
      pre_vel[0] += N(p)*un(t-1,element[ic].node[p]);
      pre_vel[1] += N(p)*vn(t-1,element[ic].node[p]);
    }
  }

  double Ad_vel[2]={0e0,0e0};
  Ad_vel[0] = 1.5*vel[0]-5e-1*pre_vel[0];
  Ad_vel[1] = 1.5*vel[1]-5e-1*pre_vel[1];
   
 

  double dudx[2]={0e0,0e0};
  double dvdx[2]={0e0,0e0};

  for(int i=0; i<numOfNodeInElm; i++){
    dudx[0] += dNdx(i,0)*un(element[ic].node[i]);
    dvdx[0] += dNdx(i,0)*vn(element[ic].node[i]);
    dudx[1] += dNdx(i,1)*un(element[ic].node[i]);
    dvdx[1] += dNdx(i,1)*vn(element[ic].node[i]);
  }
  
  /*
  double h = sqrt(dx*dx+dy*dy);
  double uMag = sqrt(Ad_vel[0]*Ad_vel[0]+Ad_vel[1]*Ad_vel[1]);
  //double Renum_e = (uMag*h*Renum)/(2*1e0*1e0);
  double Renum_e = (uMag*h*Renum)/(2e0*1e0*1e0);

  double tau=pow(2/dt,2e0)+pow(2e0*uMag/h,2e0)+pow(4e0/(Renum_e*h*h),2e0);
  //double tau=pow(2/dt,2e0)+pow(2e0*uMag/h,2e0)+9e0*pow(4e0*nu/(h*h),2e0);
  tau = 1e0/sqrt(tau);
  */

 /*
  double tau = 0e0;
  if(t!=0){
    const double norm_v = sqrt(vel[0]*vel[0]+vel[1]*vel[1]);
    double h = sqrt(dx*dx+dy*dy);
    //const double h = sqrt( area / 3.14 )*2;
    const double tau_c = h*0.5/norm_v;
    const double cou_c = norm_v*dt/h;
    const double re_c = 0.5*norm_v*h*Renum;  // 0.5*norm_v*h*rho/myu;
    const double dtmp1 = 1/(cou_c*cou_c)+1+1/(re_c*re_c);
    tau = tau_c / sqrt(dtmp1);
  }
  */
  
  

  /*
  double tau = 0e0;
  
  double term1 = (2e0/dt)*(2e0/dt);
  double term2;
  double term3;

  double uMag = sqrt(Ad_vel[0]*Ad_vel[0]+Ad_vel[1]*Ad_vel[1]);

  if(t == 0){
    term2 = 0e0;
    term3 = 0e0;
  }
  if(t != 0){
    double h;
    h = 0e0;
    for(int i=0; i<3; i++){
      h +=  abs(Ad_vel[0]*dNdx(i,0) + Ad_vel[1]*dNdx(i,1));
    }
    h = 2e0*uMag/h;

    double Renum_e = (uMag*h*Renum)/(2e0*1e0*1e0);
    term2 = pow(2e0*uMag/h,2e0);
    term3 = pow(4/(Renum_e*h*h),2e0);
  }
 
 

  //double Renum_e = (uMag*h*Renum)/(2*1e0*1e0);

  tau = pow(term1+term2+term3,-5e-1);
  */

  double tau = 0e0;
  


  
  Matrix2_2d drdx;
  drdx = dxdr.inverse();

  
  //double tau = 5e0;
  
  double term1 = (2e0/dt)*(2e0/dt);

  double G[2][2];
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      G[i][j] = 0e0;
      for(int k=0;k<2;k++) G[i][j] += drdx(k,i) * drdx(k,j);
    }
  }

  double term2 = 0e0;
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      term2 += Ad_vel[i] * G[i][j] * Ad_vel[j];
    }
  }
  
  double CI;
  CI = 36e0;
  
  double term3 = 0e0;
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      term3 += G[i][j] * G[i][j];
    }
  }
  term3 = CI * term3/Renum;
  
  tau = pow(term1+term2+term3,-5e-1);
  

  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      M(i,j) = 0e0;
      M(i+4,j+4) = 0e0;
      M(i,j) += N(i)*N(j)*detJ*weight;
      M(i+4,j+4) += N(i)*N(j)*detJ*weight;
    }
  }

  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      Ms(i,j) = 0e0;
      Ms(i+4,j+4) = 0e0;
      Ms(i,j) += ((dNdx(i,0)*Ad_vel[0]*N(j)) + (dNdx(i,1)*Ad_vel[1]*N(j)))*tau*detJ*weight;
      Ms(i+4,j+4) += ((dNdx(i,0)*Ad_vel[0]*N(j)) + (dNdx(i,1)*Ad_vel[1]*N(j)))*tau*detJ*weight;
    }
  }

  for(int i=0; i<numOfNodeInElm; i++){
    A(i) = 0e0;
    A(i+4) = 0e0;
    A(i) += ((N(i)*Ad_vel[0]*dudx[0]) + (N(i)*Ad_vel[1]*dudx[1]))*detJ*weight;
    A(i+4) += ((N(i)*Ad_vel[0]*dvdx[0]) + (N(i)*Ad_vel[1]*dvdx[1]))*detJ*weight;
  }

  ARRAY2D<double> Al(numOfNodeInElm*2,numOfNodeInElm*2);
  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      Al(i,j) = 0e0;
      Al(i+4,j+4) = 0e0;
      Al(i,j) += (N(i)*Ad_vel[0]*dNdx(j,0)+ N(i)*Ad_vel[1]*dNdx(j,1))*detJ*weight;
      Al(i+4,j+4) += (N(i)*Ad_vel[0]*dNdx(j,0) + N(i)*Ad_vel[1]*dNdx(j,1))*detJ*weight;
    }
  }


  ARRAY2D<double> D11(numOfNodeInElm,numOfNodeInElm);
  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      D11(i,j) = 0e0;
      D11(i,j) += ((2*dNdx(i,0)*dNdx(j,0)) + (dNdx(i,1)*dNdx(j,1)))*detJ*weight/Renum;
    }
  }

  ARRAY1D<double> D11r(numOfNodeInElm);
  for(int i=0; i<numOfNodeInElm; i++){
    D11r(i) = 0e0;
    D11r(i) += ((2*dNdx(i,0)*dudx[0]) + (dNdx(i,1)*dudx[1]))*detJ*weight/Renum;
  }

  ARRAY2D<double> D12(numOfNodeInElm,numOfNodeInElm);
  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      D12(i,j) = 0e0;
      D12(i,j) += dNdx(i,1)*dNdx(j,0)*detJ*weight/Renum;
    }
  }

  
  ARRAY1D<double> D12r(numOfNodeInElm);
  for(int i=0; i<numOfNodeInElm; i++){
    D12r(i) = 0e0;
    D12r(i) += dNdx(i,1)*dvdx[0]*detJ*weight/Renum;
  }

  ARRAY2D<double> D21(numOfNodeInElm,numOfNodeInElm);
  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      D21(i,j) = 0e0;
      D21(i,j) += dNdx(i,0)*dNdx(j,1)*detJ*weight/Renum;
    }
  }

  ARRAY1D<double> D21r(numOfNodeInElm);
  for(int i=0; i<numOfNodeInElm; i++){
    D21r(i) = 0e0;
    D21r(i) += dNdx(i,0)*dudx[1]*detJ*weight/Renum;
  }

  ARRAY2D<double> D22(numOfNodeInElm,numOfNodeInElm);
  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      D22(i,j) = 0e0;
      D22(i,j) += ((dNdx(i,0)*dNdx(j,0)) + (2*dNdx(i,1)*dNdx(j,1)))*detJ*weight/Renum;
    }
  }

  ARRAY1D<double> D22r(numOfNodeInElm);
  for(int i=0; i<numOfNodeInElm; i++){
    D22r(i) = 0e0;
    D22r(i) += ((dNdx(i,0)*dvdx[0]) + (2*dNdx(i,1)*dvdx[1]))*detJ*weight/Renum;
  }
  
  
  for(int i=0; i<numOfNodeInElm; i++){
    As(i) = 0e0;
    As(i+4) = 0e0;
    As(i) += ((dNdx(i,0)*Ad_vel[0]*Ad_vel[0]*dudx[0]) + (dNdx(i,1)*Ad_vel[1]*Ad_vel[0]*dudx[0]))*tau*detJ*weight;
    As(i) += ((dNdx(i,0)*Ad_vel[0]*Ad_vel[1]*dudx[1]) + (dNdx(i,1)*Ad_vel[1]*Ad_vel[1]*dudx[1]))*tau*detJ*weight;

    As(i+4) += ((dNdx(i,0)*Ad_vel[0]*Ad_vel[0]*dvdx[0]) + (dNdx(i,1)*Ad_vel[1]*Ad_vel[0]*dvdx[0]))*tau*detJ*weight;
    As(i+4) += ((dNdx(i,0)*Ad_vel[0]*Ad_vel[1]*dvdx[1]) + (dNdx(i,1)*Ad_vel[1]*Ad_vel[1]*dvdx[1]))*tau*detJ*weight;
  }

  ARRAY2D<double> Als(numOfNodeInElm*2,numOfNodeInElm*2);
  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      Als(i,j) = 0e0;
      Als(i+4,j+4) = 0e0;
      Als(i,j) += ((dNdx(i,0)*Ad_vel[0]*Ad_vel[0]*dNdx(j,0)) + (dNdx(i,1)*Ad_vel[1]*Ad_vel[0]*dNdx(j,0)))*tau*detJ*weight;
      Als(i,j) += ((dNdx(i,0)*Ad_vel[0]*Ad_vel[1]*dNdx(j,1)) + (dNdx(i,1)*Ad_vel[1]*Ad_vel[1]*dNdx(j,1)))*tau*detJ*weight;

      Als(i+4,j+4) += ((dNdx(i,0)*Ad_vel[0]*Ad_vel[0]*dNdx(j,0)) + (dNdx(i,1)*Ad_vel[1]*Ad_vel[0]*dNdx(j,0)))*tau*detJ*weight;
      Als(i+4,j+4) += ((dNdx(i,0)*Ad_vel[0]*Ad_vel[1]*dNdx(j,1)) + (dNdx(i,1)*Ad_vel[1]*Ad_vel[1]*dNdx(j,1)))*tau*detJ*weight;
    }
  }

  for(int i=0; i<numOfNodeInElm; i++){
    Mr(i) = 0e0;
    Mr(i+4) = 0e0;
    Mr(i) += N(i)*vel[0]*detJ*weight;
    Mr(i+4) += N(i)*vel[1]*detJ*weight;
  }

  for(int i=0; i<numOfNodeInElm; i++){
    Mrs(i) = 0e0;
    Mrs(i+4) = 0e0;
    Mrs(i) += ((dNdx(i,0)*Ad_vel[0]*vel[0]) + (dNdx(i,1)*Ad_vel[1]*vel[0]))*tau*detJ*weight;
    Mrs(i+4) += ((dNdx(i,0)*Ad_vel[0]*vel[1]) + (dNdx(i,1)*Ad_vel[1]*vel[1]))*tau*detJ*weight;
  }


  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
    Ke(i,j) = 0e0;
    Ke(i+4,j+4) = 0e0;
    Ke(i,j+4) = 0e0;
    Ke(i+4,j) = 0e0;
    Ke(i,j) += M(i,j) + 5e-1*D11(i,j)*dt + 5e-1*Al(i,j)*dt + 5e-1*Als(i,j)*dt + Ms(i,j);
    Ke(i,j+4) += 5e-1*D12(i,j)*dt;
    Ke(i+4,j) += 5e-1*D21(i,j)*dt;
    Ke(i+4,j+4) +=  M(i+4,j+4) + 5e-1*D22(i,j)*dt + 5e-1*Al(i+4,j+4)*dt + 5e-1*Als(i+4,j+4)*dt + Ms(i+4,j+4);
    }
  }
  
 
  /*
  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      Ke(i,j) = 0e0;
      Ke(i+4,j+4) = 0e0;
      Ke(i,j) = M(i,j);
      Ke(i+4,j+4) =  M(i+4,j+4);
    }
  }*/
  
  for(int i=0; i<numOfNodeInElm; i++){
    Re(i) = -dt*(5e-1*A(i)+5e-1*D11r(i)+ 5e-1*As(i))+Mr(i)+ Mrs(i) - 5e-1*D12r(i)*dt;
    Re(i+4) = -dt*(5e-1*A(i+4)+5e-1*D22r(i)+5e-1*As(i+4))+Mr(i+4)+ Mrs(i+4) - 5e-1*D21r(i)*dt;
  }
 
   /*
  for(int i=0; i<numOfNodeInElm; i++){
    Re(i) = 0e0;
    Re(i+4) = 0e0;
    Re(i) = -dt*D(i)+Mr(i);
    Re(i+4) = -dt*D(i+4)+Mr(i+4);
  }*/

  if(ic == numOfElm-1){
    /*
    cout << "dt*D" << endl;
    cout << dt*D << endl;
    cout << endl;

    cout << "Mrs" << endl;
    cout << Mrs << endl;
    cout << endl;

    cout << "Mr" << endl;
    cout << Mr << endl;
    cout << endl;
    */  
   /*
    cout << "dNdr" << endl;
    cout << dNdr << endl;
    cout << endl;

    cout << "dxdr" << endl;
    cout << dxdr << endl;
    cout << endl;

    cout << "drdx" << endl;
    cout << dxdr.inverse() << endl;
    cout << endl;

    cout << "x(0)" << " " << "x(1)" << " " <<  "x(2)"  << " " <<  "x(3)" <<endl;
    cout << x_current(0,0) << " " << x_current(1,0) << " " << x_current(2,0) << " " << x_current(3,0)<< endl;
    cout << endl;

    cout << "y(0)" << " " << "y(1)" << " " <<  "y(2)"  << " " <<  "y(3)" <<endl;
    cout << x_current(0,1) << " " << x_current(1,1) << " " << x_current(2,1) << " " << x_current(3,1)<< endl;
    cout << endl;


    cout << "dNdx" << endl;
    cout << dNdx << endl;
    cout << endl;
    */
   /*
    cout << "tau" << endl;
    cout << tau << endl;
    cout << endl;
   
    cout << "dudx" << endl;
    cout << dudx[0] << endl;
    cout << endl;

    cout << "dudy" << endl;
    cout << dudx[1] << endl;
    cout << endl;

    cout << "detJ" << endl;
    cout << detJ << endl;
    cout << endl;

    cout << "weight" << endl;
    cout << weight << endl;
    cout  << endl;

    cout << "vel" << endl;
    cout << vel[0] << endl;
    cout << endl;
    */
  
   /*
    ofstream outputfile_A("A.txt");
    outputfile_A << A;
    outputfile_A.close();

    ofstream outputfile_D("D.txt");
    outputfile_D << D;
    outputfile_D.close();

    ofstream outputfile_As("As.txt");
    outputfile_As << As;
    outputfile_As.close();
    */
  }
  

  /*
  for(int i=0; i<numOfNodeInElm; i++){
    Re(i) = -dt*A(i)+Mr(i);
    Re(i+4) = -dt*A(i+4)+Mr(i+4);
  }*/
  
  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      K(element[ic].node[i], element[ic].node[j]) += Ke(i, j);
      K(element[ic].node[i]+x_size, element[ic].node[j]) += Ke(i+4, j);
      K(element[ic].node[i], element[ic].node[j]+x_size) += Ke(i, j+4);
      K(element[ic].node[i]+x_size, element[ic].node[j]+x_size) += Ke(i+4, j+4);
    }
  }

  for(int i=0; i<numOfNodeInElm; i++){
    R(element[ic].node[i]) += Re(i);
    R(element[ic].node[i]+x_size) += Re(i+4);
  }

}

void FEMF::bc(MatrixXd &K,VectorXd &R,int x,double u_value,double v_value,int check_v){
  for (int i = x; i < x+boundary_element_size; i++)
  {
    for (int j = 0; j < x_size * 2; j++)
    {
      K(boundaryElement[i].node[0], j) = 0e0;
      //if(check_v == 1){
        K(boundaryElement[i].node[0]+x_size,j) = 0e0;
      //}

      if(i == x+boundary_element_size-1){
        K(boundaryElement[i].node[1], j) = 0e0;
        //if(check_v == 1){
        K(boundaryElement[i].node[1]+x_size, j) = 0e0;
        //}
      }
    }
    K(boundaryElement[i].node[0], boundaryElement[i].node[0]) = 1e0;
    //if(check_v == 1){
      K(boundaryElement[i].node[0] + x_size, boundaryElement[i].node[0] + x_size) = 1e0;
    //}
  
    if(i == x+boundary_element_size-1){
      K(boundaryElement[i].node[1], boundaryElement[i].node[1]) = 1e0;
      //if(check_v == 1){
        K(boundaryElement[i].node[1] + x_size, boundaryElement[i].node[1] + x_size) = 1e0;
      //}
      
    }

    R(boundaryElement[i].node[0]) = u_value;
    //if(check_v == 1){
      R(boundaryElement[i].node[0] + x_size) = v_value;
    //}
   
    if(i == x+boundary_element_size-1){
      R(boundaryElement[i].node[1]) = u_value;
       //if(check_v == 1){
        R(boundaryElement[i].node[1] + x_size) = v_value;
       //}
      
    }
  }
}

int FEMF::solve_Matrix(MatrixXd &K,VectorXd &R,VectorXd &a,int s)
{
  SparseMatrix<double> K_sparse(s, s);
  vector<T> tripletList;

  for (int i = 0; i < s; i++)
  {
    for (int j = 0; j < s; j++)
    {
      if (K(i, j) != 0)
      {
        tripletList.push_back(T(i, j, K(i, j)));
      }
    }
  }
  K_sparse.setFromTriplets(tripletList.begin(), tripletList.end());

  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
  solver.compute(K_sparse);
  if (solver.info() != Success)
  {
    // decomposition failed
    cout << "decomposition failed" << endl;
    return -1;
  }

  a = solver.solve(R);
  if (solver.info() != Success)
  {
    // solving failed
    std::cout << "solving failed" << endl;
    return -1;
  }


}

void FEMF::poissonMatrix(Vector4d &umid,Vector4d &vmid,MatrixXd &K2,VectorXd &R2,Vector4d &N,Matrix4_2d &dNdr,Matrix4_2d &x_current,const int numOfNodeInElm,const double weight,const int ic)
{
  Matrix4_2d dNdx;
  Matrix2_2d dxdr;
  Matrix4_4d Ke2;
  Vector4d Re2;

  double umiddNdx[2]={0e0,0e0};
  double vmiddNdx[2]={0e0,0e0};


  calc_dxdr2D(dxdr,dNdr,x_current,numOfNodeInElm);
  double detJ = dxdr.determinant();
  calc_dNdx2D(dNdx,dNdr,dxdr,numOfNodeInElm);

  for(int i=0; i<numOfNodeInElm; i++){
    umiddNdx[0] += dNdx(i,0)*umid(i);
    vmiddNdx[0] += dNdx(i,0)*vmid(i);
    umiddNdx[1] += dNdx(i,1)*umid(i);
    vmiddNdx[1] += dNdx(i,1)*vmid(i);
  }

  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      Ke2(i,j) = 0e0;
      for(int p=0; p<2; p++){
        Ke2(i,j) += dNdx(i,p)*dNdx(j,p)*detJ*weight;
      }
    }
  }

  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      K2(element[ic].node[i],element[ic].node[j]) += Ke2(i,j);
    }
  }

  for(int i=0; i<numOfNodeInElm; i++){
    Re2(i) = 0e0;
    Re2(i) += -(N(i)*umiddNdx[0] + N(i)*vmiddNdx[1])*detJ*weight/dt;
  }

  for(int i=0; i<numOfNodeInElm; i++){
    R2(element[ic].node[i]) += Re2(i);
  }

}

void FEMF::u_nextMatrix(Vector4d &umid,Vector4d &vmid,Vector4d &p,MatrixXd &K3,VectorXd &R3,Vector4d &N,Matrix4_2d &dNdr,Matrix4_2d &x_current,const int numOfNodeInElm,const double weight,const int ic)
{
  Matrix4_2d dNdx;
  Matrix2_2d dxdr;
  Matrix8_8d Ke3;
  Vector8d Re3;


  

  calc_dxdr2D(dxdr,dNdr,x_current,numOfNodeInElm);
  double detJ = dxdr.determinant();
  calc_dNdx2D(dNdx,dNdr,dxdr,numOfNodeInElm);

  double dpdx[2] = {0e0,0e0};

  for(int i=0;i<numOfNodeInElm;i++){
    dpdx[0] += dNdx(i,0)*p(i);
    dpdx[1] += dNdx(i,1)*p(i);
  }/*
  if(ic == numOfElm-1){
    cout << "dNdx" << endl;
    cout << dNdx << endl;
    cout << endl;
    cout << "N" << endl;
    cout << N << endl;
    cout << endl;
  }*/


  double vel[2]={0e0,0e0};
  for(int p=0;p<numOfNodeInElm;p++){
    vel[0] += N(p)*umid(p);
    vel[1] += N(p)*vmid(p);
  }

   for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      Ke3(i,j) = 0e0;
      Ke3(i+4,j+4) = 0e0;
      Ke3(i,j) += N(i)*N(j)*detJ*weight;
      Ke3(i+4,j+4) += N(i)*N(j)*detJ*weight;
    }
  }

  for(int i=0; i<numOfNodeInElm; i++){
    for(int j=0; j<numOfNodeInElm; j++){
      K3(element[ic].node[i],element[ic].node[j]) += Ke3(i,j);
      K3(element[ic].node[i]+x_size, element[ic].node[j]+x_size) += Ke3(i+4, j+4);
    }
  }

  for(int i=0; i<numOfNodeInElm; i++){
    Re3(i) = 0e0;
    Re3(i+4) = 0e0;
    Re3(i) += (N(i)*vel[0] - dt*N(i)*dpdx[0])*detJ*weight;
    Re3(i+4) += (N(i)*vel[1] - dt*N(i)*dpdx[1])*detJ*weight;
  } 

  for(int i=0; i<numOfNodeInElm; i++){
    R3(element[ic].node[i]) += Re3(i);
    R3(element[ic].node[i]+x_size) += Re3(i+4);
  }
   /*
  if(ic == numOfElm-1){
    
    cout << "p" << endl;
    cout << p << endl;
    cout << endl;

    cout << "dNdx" << endl;
    cout << dNdx << endl;
    cout << endl;

    cout << "dpdx" << endl;
    cout << dpdx[0] << endl;
    cout << endl;

    cout << "dpdy" << endl;
    cout << dpdx[1] << endl;
    cout << endl;

  }*/

  
}
void FEMF::export_vtu(const std::string &file,DataType dataType,VectorXd &U2,VectorXd &U3)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", numOfNode, numOfElm);
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ic=0;ic<numOfNode;ic++){
    fprintf(fp,"%e %e 0e0\n",x(ic,0),x(ic,1));
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < numOfElm; i++){
    for (int j = 0; j < element[i].node.size(); j++) fprintf(fp, "%d ", element[i].node[j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < numOfElm; i++)
  {
    num += element[i].node.size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < numOfElm; i++) fprintf(fp, "%d\n", element[i].meshType);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  if(dataType==DataType::VELOCITY){
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"velocity[m/s]\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int ic=0;ic<numOfNode;ic++){
      fprintf(fp,"%e %e 0e0\n",U3(ic),U3(ic+x_size));
    }
    fprintf(fp, "</DataArray>\n");
  }else if(dataType==DataType::PRESSURE){
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(int ic=0;ic<numOfNode;ic++){
      fprintf(fp,"%e\n",U2(ic));
    }
    fprintf(fp, "</DataArray>\n");
  }
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  /*
  if(dataType==DataType::PHI){
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(int ic=0;ic<numOfElm;ic++){
      fprintf(fp,"%e\n",phi(ic));
    }
    fprintf(fp, "</DataArray>\n");
  }*/
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}
