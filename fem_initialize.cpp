#include "fem_fluid.h"

using namespace std;

void FEMF::initialize()
{
  
  inputDomainInfo();
  inputBoundaryInfo();
  
}

void FEMF::inputDomainInfo()
{
  nx=32;
  ny=32;

  Lx=1e0;
  Ly=1e0;

  dx=Lx/nx;
  dy=Ly/ny;

  Renum = 100;
  nu = 0.01;
  rho = 100;
  mu = 0.01;

  dt = 0.1;
  time = 200;

  numOfNode=(nx+1)*(ny+1);
  numOfElm=nx*ny;

  x.allocate(numOfNode,2);
  element.resize(numOfElm);

  for(int ic=0;ic<numOfElm;ic++){
    element[ic].meshType = VTK_QUAD;
    element[ic].node.resize(4);
    element[ic].numOfGaussPoint = 2;
  }

  int tmp2=0;
  
  for(int i=0;i<ny+1;i++){
    for(int j=0;j<nx+1;j++){
      x(tmp2,0)=j*dx;
      x(tmp2,1)=i*dy;
      tmp2++;
    }
  }

  x_size = tmp2;

  tmp2=0;
  for(int i=0;i<ny;i++){
    for(int j=0;j<nx;j++){
      element[tmp2].node[0]=j   +i*(nx+1);
      element[tmp2].node[1]=j+1 +i*(nx+1);
      element[tmp2].node[2]=j+1 +(i+1)*(nx+1);
      element[tmp2].node[3]=j   +(i+1)*(nx+1);
      tmp2++;
    }
  }
}

void FEMF::inputBoundaryInfo()
{
  //boundary information
  boundaryElement.resize(2*nx+2*ny);

    //bottom
    int tmp=0;
    for(int i=0;i<nx;i++){
      boundaryElement[tmp].meshType = VTK_LINE;
      boundaryElement[tmp].node.resize(2);
      boundaryElement[tmp].node[0] = i;
      boundaryElement[tmp].node[1] = i+1;
      tmp++;
    }
    boundary_element_size = tmp;

    //top
    for(int i=0;i<nx;i++){
      boundaryElement[tmp].meshType = VTK_LINE;
      boundaryElement[tmp].node.resize(2);
      boundaryElement[tmp].node[0] = i+(nx+1)*ny;
      boundaryElement[tmp].node[1] = i+1+(nx+1)*ny;
      tmp++;
    }

    //left
    for(int i=0;i<ny;i++){
      boundaryElement[tmp].meshType = VTK_LINE;
      boundaryElement[tmp].node.resize(2);
      boundaryElement[tmp].node[0] = (nx+1)*i;
      boundaryElement[tmp].node[1] = (nx+1)*(i+1);
      tmp++;
    }

    //right
      for(int i=0;i<ny;i++){
      boundaryElement[tmp].meshType = VTK_LINE;
      boundaryElement[tmp].node.resize(3);
      boundaryElement[tmp].node[0] = nx+(nx+1)*i;
      boundaryElement[tmp].node[1] = nx+(nx+1)*(i+1);
      tmp++;
    }
  
}