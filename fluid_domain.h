#ifndef _FLUID_DOMAIN_H_
#define _FLUID_DOMAIN_H_

#include <string>
#include "fem_define.h"
#include "allocation.h"
#include <vector>
using namespace std;

class fluidDomain{
public:
  int nx,ny;
  double dx,dy,Lx,Ly;
  ARRAY2D<double> x;
  int numOfNode,numOfElm,numOfuBd,numOfpBd,numOfBdWall;
  int meshType;
  std::vector<ElementType> element,boundaryElement;

  void set_geometry();
  void set_geometry(const std::string &file1,const std::string &file2,const std::string &file3);
  void set_boundary(const std::string &v_file,const std::string &p_file);

private:
};

#endif