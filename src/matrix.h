#ifndef ASC_MATRIX_H
#define ASC_MATRIX_H 1

// simple 2d matrix class
// from Stroustrup
// together with some other stuff

//#include "define.h"

#include <fstream>
#include <assert.h>



using namespace std;

class MatrixInt {
private:
  int d1;
  int d2;
  int * v;
public:
  MatrixInt(int x, int y);
  MatrixInt(const MatrixInt&);
  MatrixInt(ifstream &); // read it from a file stream
  void store(ofstream &) const ; // store it to a file stream
  const MatrixInt& operator=(const MatrixInt&);
  const MatrixInt& operator=(int);
  const MatrixInt& operator*=(int);
  ~MatrixInt();
  
  int size() const { return d1 * d2; }
  int dim1() const { return d1;}
  int dim2() const { return d2;}
  
  int sum() const;
  int max() const;
  
  void dump() const; // dump it to cout
  void dumpNonZero() const; // same but only non-zero elements
  
  int& operator()(int x, int y);
  int operator()(int x, int y) const;
  
  int rowSum(int x) const;
  int columnSum(int y) const;
  
  bool isPermutation() const;
  bool isIdentity() const;
  
  void times(const MatrixInt&, const MatrixInt &) ; // sets MatrixInt to a * b
  void identity(); // sets it to the identity
 private:
  int getV(int x, int y) const{
    return v[(x * d2) + y];
  }
  int& vv(int x, int y) {
    assert(x < d1 && y < d2);
    return v[(x * d2) + y];
  }
  
  void  setV(int x, int y, int v_){
    (*this)(x,y)=v_;
  }
  
};



class Matrix3 {
	double* v;
	int d1,d2,d3;
public:
	Matrix3(int x, int y, int z);
	Matrix3(const Matrix3&);
	Matrix3& operator=(const Matrix3&);
	~Matrix3();

	int size() const { return d1 * d2 *d3; }
	int dim1() const { return d1;}
	int dim2() const { return d2;}
	int dim3() const { return d3;}

	double& operator()(int x, int y, int z);
	double operator()(int x, int y, int z) const;

};

// this returns
// log( e^x + e^y)
// being careful not to mess it up

double logAdd(double x, double y);



class Matrix {
private:
  int d1;
  int d2;
  double * v;
public:
  Matrix(int x, int y);
  Matrix(const Matrix&);
  Matrix(ifstream &); // read it from a file stream
  void store(ofstream &) const ; // store it to a file stream
  const Matrix& operator=(const Matrix&);
  const Matrix& operator=(double);
  const Matrix& operator*=(double);
  ~Matrix();
  
  int size() const { return d1 * d2; }
  int dim1() const { return d1;}
  int dim2() const { return d2;}
  
  double sum() const;
  double max() const;
  
  void dump() const; // dump it to cout
  void dumpNonZero() const; // same but only non-zero elements
  
  double& operator()(int x, int y);
  double operator()(int x, int y) const;
  
  double rowSum(int x) const;
  double columnSum(int y) const;
  // divide each entry by the sum of row and column
  void normaliseRowColumn();
  
  bool isPermutation() const;
  bool isIdentity() const;
  
  Matrix* matrixMinor(int x, int y) const;
  double permanent() const;
  double mi() const;
  double logPermanent() const;
  double determinant() const;
  double diff(const Matrix&) const;
  void inverse(); // returns determinant, calculates inverse. destroys matrix
  
  // returns a new matrix with the optimum weights given that
  // row sums and column sums must be less than 1
  Matrix* simplex() const;
  void times(const Matrix&, const Matrix &) ; // sets Matrix to a * b
  void identity(); // sets it to the identity
 private:
  double getV(int x, int y) const{
    return v[(x * d2) + y];
  }
  double& vv(int x, int y) {
    assert(x < d1 && y < d2);
    return v[(x * d2) + y];
  }
  
  void  setV(int x, int y, double v_){
    (*this)(x,y)=v_;
  }
  
};

#endif

