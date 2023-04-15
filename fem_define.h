#ifndef _FEM_DEFINE_H_
#define _FEM_DEFINE_H_


#include <string>
#include <vector>
#include "vtkCellType.h"

using namespace std;

template<typename T>
using VECTOR1D = vector<T>;
template<typename T>
using VECTOR2D = vector<std::vector<T> >;

#define FEM_VERS "1.0" 
//double GRAVITY = 9.80665; // (m/s2)
//double PI = 3.1415926535897932384626;

// general
#define ON          1
#define OFF         0

//data encode
#define INT         0
#define DOUBLE      1
#define ASCII       0
#define BINARY      1

//Constitutive law
#define LinearElasticMaterial 0
#define StVenantMaterial 1
#define NeoHookean 2

// Mesh type
#define C2D4        0
#define C3D4        1
#define C3D8        2

// IO file format
#define ORIGINAL    0
#define HDF5        1

class ElementType{
 public:
  VTKCellType meshType;
  int materialType;
  int numOfGaussPoint;
  VECTOR1D<int> node;
};


#endif // _FB_DEFINE_H_