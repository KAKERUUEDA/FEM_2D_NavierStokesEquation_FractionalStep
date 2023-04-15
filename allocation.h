#ifndef _ALLOCATION_H_
#define _ALLOCATION_H_

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <type_traits>
#include <H5Cpp.h>

template <typename T>
class ARRAY1D
{
public:
    ARRAY1D()
    {
        data = new T[1];
    }
    ARRAY1D(const int X)
    {
        nx = X;
        data = new T[nx];
    }
    ~ARRAY1D()
    {
        delete[] data;
    }
    T &operator()(const int i)
    {
        return data[i];
    }
    void allocate(const int X)
    {
        delete[] data;
        nx = X;
        data = new T[nx];
    }

    void importData(const std::string &file){
      int size = sizeof(data)*(nx);
      std::fstream ofs;
      ofs.open(file.c_str(),std::ios::in | std::ios_base::binary);
      ofs.read((char *)data,size);
      ofs.close();
    }

    void exportData(const std::string &file){
      int size = sizeof(data)*(nx);
      std::fstream ofs;
      ofs.open(file.c_str(),std::ios::out | std::ios_base::binary);
      ofs.write((char *)data,size);
      ofs.close();
    }

    void importHDF5(H5::H5File &file,const std::string &dataName){
      H5::DataSet dataset = file.openDataSet(dataName.c_str());
      H5::DataSpace dataspace = dataset.getSpace();

      hsize_t dims_out[1];
      // int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
      nx = (unsigned long)(dims_out[0]);
      data = new T[nx];
      if(std::is_integral<T>()==true){
        dataset.read(&data[0], H5::PredType::NATIVE_INT);
      }else if(std::is_floating_point<T>()==true){
        dataset.read(&data[0], H5::PredType::NATIVE_DOUBLE);
      }
    }

    void exportHDF5(H5::H5File &file,const std::string &dataName){
      H5std_string DATASET_NAME(dataName.c_str());

      hsize_t dim[1] = {nx}; // dataset dimensions
      H5::DataSpace dataspace(1, dim);

      if(std::is_integral<T>()==true){
        H5::IntType datatype(H5::PredType::NATIVE_INT);
        datatype.setOrder(H5T_ORDER_LE);
        H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
        dataset.write(&data[0], H5::PredType::NATIVE_INT);
      }else if(std::is_floating_point<T>()==true){
          H5::IntType datatype(H5::PredType::NATIVE_DOUBLE);
          datatype.setOrder(H5T_ORDER_LE);
          H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
          dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
      }
    }
private:
    T *data;
    int nx;
};

template <typename T>
class ARRAY2D {
  public:
    ARRAY2D(){
        data = new T[1];
    }
    ARRAY2D(const int X, const int Y,const int GuideCell=0) {
        nx = X;
        ny = Y;
        guideCell = GuideCell;
        data = new T[(nx+2*guideCell) * (ny+2*guideCell)];
    }
    ~ARRAY2D() {
        delete[] data;
    }
    T& operator()(const int i,const int j) {
        return data[(i+guideCell)*(ny+2*guideCell) + (j+guideCell)];
    }
    void allocate(const int X, const int Y,const int GuideCell=0){
        delete[] data;
        nx = X;
        ny = Y;
        guideCell = GuideCell;
        data = new T[(nx+2*guideCell) * (ny+2*guideCell)];
    }
    void importData(const std::string &file){
      int size = sizeof(data)*(nx+2*guideCell)*(ny+2*guideCell);
      std::fstream ofs;
      ofs.open(file.c_str(),std::ios::in | std::ios_base::binary);
      ofs.read((char *)data,size);
      ofs.close();
    }
    void exportData(const std::string &file){
      int size = sizeof(data)*(nx+2*guideCell)*(ny+2*guideCell);
      std::fstream ofs;
      ofs.open(file.c_str(),std::ios::out | std::ios_base::binary);
      ofs.write((char *)data,size);
      ofs.close();
    }

    void importHDF5(H5::H5File &file,const std::string &dataName,const int GuideCell=0){

      guideCell = GuideCell;

      H5::DataSet dataset = file.openDataSet(dataName.c_str());
      H5::DataSpace dataspace = dataset.getSpace();

      hsize_t dims_out[2];
      // int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
      nx = (unsigned long)(dims_out[0])-2*guideCell;
      ny = (unsigned long)(dims_out[1])-2*guideCell;

      data = new T[(nx+2*guideCell) * (ny+2*guideCell)];

      if(std::is_integral<T>()==true){
        dataset.read(&data[0], H5::PredType::NATIVE_INT);
      }else if(std::is_floating_point<T>()==true){
        dataset.read(&data[0], H5::PredType::NATIVE_DOUBLE);
      }
    }

    void exportHDF5(H5::H5File &file,const std::string &dataName){
      H5std_string DATASET_NAME(dataName.c_str());

      hsize_t dim[2] = {nx+2*guideCell,ny+2*guideCell}; // dataset dimensions
      H5::DataSpace dataspace(2, dim);

      if(std::is_integral<T>()==true){
        H5::IntType datatype(H5::PredType::NATIVE_INT);
        datatype.setOrder(H5T_ORDER_LE);
        H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
        dataset.write(&data[0], H5::PredType::NATIVE_INT);
      }else if(std::is_floating_point<T>()==true){
          H5::IntType datatype(H5::PredType::NATIVE_DOUBLE);
          datatype.setOrder(H5T_ORDER_LE);
          H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
          dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
      }
    }

  private:
    T* data;
    int nx;
    int ny;
    int guideCell;
};

template <typename T>
class ARRAY3D {
  public:
    ARRAY3D(){
        data=new T[1];
    }
    ARRAY3D(const int X, const int Y, const int Z,const int GuideCell=0) {
        nz = Z;
        nx = X;
        ny = Y;
        guideCell = GuideCell;
        data = new T[(nz+2*guideCell) * (nx+2*guideCell) * (ny+2*guideCell)];
    }
    ~ARRAY3D() {
        delete[] data;
    }
    T& operator()(const int i,const int j,const int k) {
        return data[(i+guideCell)*(ny+2*guideCell)*(nz+2*guideCell) 
                  + (j+guideCell)*(nz+2*guideCell) 
                  + (k+guideCell)];
    }
    void allocate(const int X, const int Y, const int Z,const int GuideCell=0){
        delete[] data;
        nz = Z;
        nx = X;
        ny = Y;
        guideCell = GuideCell;
        data = new T[(nz+2*guideCell) * (nx+2*guideCell) * (ny+2*guideCell)];
    }
    void importData(const std::string &file){
      int size = sizeof(data)*(nz+2*guideCell)*(nx+2*guideCell)*(ny+2*guideCell);
      std::fstream ofs;
      ofs.open(file.c_str(),std::ios::in | std::ios_base::binary);
      ofs.read((char *)data,size);
      ofs.close();
    }
    void exportData(const std::string &file){
      int size = sizeof(data)*(nz+2*guideCell)*(nx+2*guideCell)*(ny+2*guideCell);
      std::fstream ofs;
      ofs.open(file.c_str(),std::ios::out | std::ios_base::binary);
      ofs.write((char *)data,size);
      ofs.close();
    }

    #ifdef HDF5
    void importHDF5(H5::H5File &file,const std::string &dataName,const int GuideCell=0){

      guideCell = GuideCell;

      H5::DataSet dataset = file.openDataSet(dataName.c_str());
      H5::DataSpace dataspace = dataset.getSpace();

      hsize_t dims_out[3];
      // int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
      nx = (unsigned long)(dims_out[0])-2*guideCell;
      ny = (unsigned long)(dims_out[1])-2*guideCell;
      nz = (unsigned long)(dims_out[2])-2*guideCell;

      data = new T[(nz+2*guideCell) * (nx+2*guideCell) * (ny+2*guideCell)];

      if(std::is_integral<T>()==true){
        dataset.read(&data[0], H5::PredType::NATIVE_INT);
      }else if(std::is_floating_point<T>()==true){
        dataset.read(&data[0], H5::PredType::NATIVE_DOUBLE);
      }
    }

    void exportHDF5(H5::H5File &file,const std::string &dataName){
      H5std_string DATASET_NAME(dataName.c_str());

      hsize_t dim[3] = {nx+2*guideCell,ny+2*guideCell,nz+2*guideCell}; // dataset dimensions
      H5::DataSpace dataspace(3, dim);

      if(std::is_integral<T>()==true){
        H5::IntType datatype(H5::PredType::NATIVE_INT);
        datatype.setOrder(H5T_ORDER_LE);
        H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
        dataset.write(&data[0], H5::PredType::NATIVE_INT);
      }else if(std::is_floating_point<T>()==true){
          H5::IntType datatype(H5::PredType::NATIVE_DOUBLE);
          datatype.setOrder(H5T_ORDER_LE);
          H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
          dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
      }
    }
    #endif

  private:
    T* data;
    int nx;
    int ny;
    int nz;
    int guideCell;
};

template <typename T>
class ARRAY4D {
  public:
    ARRAY4D(){
        data=new T[1];
    }
    ARRAY4D(const int X, const int Y, const int Z,const int H,const int GuideCell=0) {
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        guideCell = GuideCell;
        data = new T[(nh+2*guideCell)* (nz+2*guideCell) * (nx+2*guideCell) * (ny+2*guideCell)];
    }
    ~ARRAY4D() {
        delete[] data;
    }
    T& operator()(const int i,const int j,const int k,const int h) {
        return data[(i+guideCell)*(ny+2*guideCell)*(nz+2*guideCell)*(nh+2*guideCell) 
                  + (j+guideCell)*(nz+2*guideCell)*(nh+2*guideCell) 
                  + (k+guideCell)*(nh+2*guideCell)
                  + (h+guideCell)];
    }
    void allocate(const int X, const int Y, const int Z,const int H,const int GuideCell=0){
        delete[] data;
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        guideCell = GuideCell;
        data = new T[(nh+2*guideCell)* (nz+2*guideCell) * (nx+2*guideCell) * (ny+2*guideCell)];
    }
  private:
    T* data;
    int nx;
    int ny;
    int nz;
    int nh;
    int nl;
    int guideCell;
};

template <typename T>
class ARRAY5D {
  public:
    ARRAY5D(){
        data=new T[1];
    }
    ARRAY5D(const int X, const int Y, const int Z,const int H,const int L,const int GuideCell=0) {
        nl = L;
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        guideCell = GuideCell;
        data = new T[(nl+2*guideCell)*(nh+2*guideCell)* (nz+2*guideCell) * (nx+2*guideCell) * (ny+2*guideCell)];
    }
    ~ARRAY5D() {
        delete[] data;
    }
    T& operator()(const int i,const int j,const int k,const int h,const int l) {
        return data[(i+guideCell)*(ny+2*guideCell)*(nz+2*guideCell)*(nh+2*guideCell)*(nl+2*guideCell) 
                  + (j+guideCell)*(nz+2*guideCell)*(nh+2*guideCell)*(nl+2*guideCell) 
                  + (k+guideCell)*(nh+2*guideCell)+(nl+2*guideCell)
                  + (h+guideCell)*(nl+2*guideCell)
                  + (l+guideCell)];
    }
    void allocate(const int X, const int Y, const int Z,const int H,const int L,const int GuideCell=0){
        delete[] data;
        nl = L;
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        guideCell = GuideCell;
        data = new T[(nl+2*guideCell)*(nh+2*guideCell)* (nz+2*guideCell) * (nx+2*guideCell) * (ny+2*guideCell)];
    }
  private:
    T* data;
    int nx;
    int ny;
    int nz;
    int nh;
    int nl;
    int guideCell;
};

#endif // _ALLOCATION_H_
