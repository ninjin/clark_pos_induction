#include "matrix.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <iostream>



const double LOG_LDBL_MAX = log(LDBL_MAX);
const double LOG_LDBL_EPSILON = log(LDBL_EPSILON);


// use
// log(x+y) = log y + log(1 + exp(log x - log y))
// logAdd(log(x),log(y)) = log(x+y)

double logAdd(double logx, double logy)
{
  double diff = logx - logy;
  if (diff  > (-1 * LOG_LDBL_EPSILON))
    return logx;
  else if (diff < LOG_LDBL_EPSILON)
    return logy;
  else
    return logy  + log(1 + exp(diff));
}

//
//
// Stuff for Matrix of doubles
//
//

// copy constructor

Matrix::Matrix(const Matrix & other)
{
  d1 = other.d1;
  d2 = other.d2;
  int size = d1 * d2;
  v = new double[size];
  for (int i = 0; i < size; i++)
    v[i] = other.v[i];
}
  

Matrix::Matrix(int x, int y)
{
  //  cout << "Allocating matrix (" << x << "," << y << ")" << endl;
  assert(x > 0);
  assert(y > 0);
  d1 = x;
  d2 = y;
  int size = x * y;
  v = new double[size];
  for (int i = 0; i < size; i++)
    v[i] = 0.0l;
}

Matrix::~Matrix()
{
  delete[] v;
}



double Matrix::diff(const Matrix & other) const
{
  double d = 0.0l;
  for (int x=0; x < d1; x++){
    for (int y = 0; y < d2; y++){
      d+= fabs(getV(x,y) - other.getV(x,y));
    }
  }
  return d;
}

void Matrix::dump() const
{
  for (int i = 0; i < d1; i++){
    for (int j = 0; j < d2; j++){
      cout << "(" << i << "," << j << ")=" << getV(i,j) << " ";
    }
    cout << endl;
  }
}

void Matrix::dumpNonZero() const
{
  for (int i = 0; i < d1; i++){
    for (int j = 0; j < d2; j++){
      if (getV(i,j)>0)
	cout << "(" << i << "," << j << ")=" << getV(i,j) << " ";
    }
    cout << endl;
  }
}

bool Matrix::isPermutation() const
{
  if (d1 != d2)
    return false;
  for (int i = 0; i < d1; i++){
    for (int j = 0; j < d2; j++){
      double value = getV(i,j);
      if (value > LDBL_EPSILON && value != 1.0l)
	return false;
    }
  }
  return true;
//   for (int i = 0; i < d1; i++){
//     if (rowSum(i) != 1.0l)
//       return false;
//     if (columnSum(i) != 1.0l)
//       return false;
//   }
//   return true;
}

bool Matrix::isIdentity() const
{
  if (d1 != d2)
    return false;
  for (int i = 0; i < d1; i++){
    for (int j = 0; j < d2; j++){
      double value = getV(i,j);
      if (i==j && value!=1.0l)
	return false;
      if (i!=j && value!=0.0l)
	return false;
    }
  }
  return true;
}

double Matrix::rowSum(int x) const
{
  assert(x >= 0 && x < d1);
  double sum = 0.0l;
  for (int y = 0; y < d2; y++){
    sum += getV(x,y);
  }
  return sum;
}	

double Matrix::columnSum(int y) const
{
  assert(y >= 0 && y < d2);
  double sum = 0.0l;
  for (int x = 0; x < d1; x++){
    sum += getV(x,y);
  }
  return sum;
}

const Matrix& Matrix::operator=(const Matrix& other)
{
  assert(d1 == other.d1);
  assert(d2 == other.d2);
  int n = d1*d2;
  for (int i = 0;i< n;i++)
    v[i]= other.v[i];
  return *this;
}

const Matrix& Matrix::operator=(double newValue)
{
  int n = d1*d2;
  for (int i = 0;i< n;i++)
    v[i]= newValue;
  return *this;
}


const Matrix& Matrix::operator*=(double newValue)
{
	int n = d1*d2;
	for (int i = 0;i< n;i++)
	  v[i] *= newValue;
	return *this;
}

double Matrix::operator ()(int x, int y) const
{

	assert(x >= 0 && x < d1);
	assert(y >= 0 && y < d2);
	return v[(x * d2) + y];	
}

double& Matrix::operator ()(int x, int y)
{
	assert(x >= 0 && x < d1);
	assert(y >= 0 && y < d2);
	return v[(x * d2) + y];	
}

// return a crude upper bound of the permanent

double Matrix::permanent() const
{
	assert(d1 == d2);
	dump();
	vector<double> rowSums(d1);
	vector<double> colSums(d1);
	for (int x = 0; x<d1; x++){
	  rowSums[x] = 0.0l;
	  colSums[x] = 0.0l;
	}
	for (int x = 0; x<d1; x++){
	  for (int y = 0; y < d1; y++){
	    double value = getV(x,y);
	    rowSums[x] += value;
	    colSums[y] += value;
	  }
	}
	double rowEst = 1.0l;
	double colEst = 1.0l;
	for (int x = 0; x<d1; x++){
	  rowEst *= rowSums[x];
	  colEst *= colSums[x];
	}
	cout << "Row est = " << rowEst << endl;
	cout << "Col est = " << colEst << endl;

	return (colEst < rowEst) ? colEst : rowEst;
}



double Matrix::logPermanent() const
{
	assert(d1 == d2);
//	dump();
	vector<double> rowSums(d1);
	vector<double> colSums(d1);
	for (int x = 0; x<d1; x++){
	  rowSums[x] = 0.0l;
	  colSums[x] = 0.0l;
	}
	for (int x = 0; x<d1; x++){
	  for (int y = 0; y < d1; y++){
	    double vv = getV(x,y);
	    rowSums[x] += vv;
	    colSums[y] += vv;
	  }
	}
	double rowEst = 0.0l;
	double colEst = 0.0l;
	for (int x = 0; x<d1; x++){
	  rowEst += log(rowSums[x]);
	  colEst += log(colSums[x]);
	}
//	cout << "Row est = " << rowEst << endl;
//	cout << "Col est = " << colEst << endl;
	return (colEst < rowEst) ? colEst : rowEst;
}





Matrix3::Matrix3(int x, int y, int z)
{
	assert(x > 0);
	assert(y > 0);
	assert(z>0);
	d1 = x;
	d2 = y;
	d3 = z;
	int size=x*y*z;
	v = new double[size];
	for (int i = 0; i < size; i++)
	  v[i] = 0.0l;
}

Matrix3::~Matrix3()
{
  delete[] v;
}

double Matrix3::operator ()(int x, int y, int z) const
{
  assert(x >= 0 && x < d1);
  assert(y >= 0 && y < d2);
  assert(z >= 0 && z < d3);
  return v[(x * d2 *d3) + y *d3 + z];	
}

double& Matrix3::operator ()(int x, int y, int z)
{
  assert(x >= 0 && x < d1);
  assert(y >= 0 && y < d2);
  assert(z >= 0 && z < d3);
  return v[(x * d2 *d3) + y *d3 + z];	
}

double Matrix::sum() const
{
  double s = 0.0l;
  int l = size();
  for (int i = 0; i <l; i++)
    s+= v[i];
  return s;
}


double Matrix::max() const
{
  double s = 0.0l;
  int l = size();
  for (int i = 0; i <l; i++){
    double vv = v[i];
    if (vv > s)
      s=vv;
  }
  return s;
}

// destructively sets it to the 
// the identity matrix;

void Matrix::identity()
{
  assert(d1 == d2);
  for (int i = 0; i < d1; i++){
    for (int j = 0; j < d1; j++){
      if (i == j)
	setV(i,j,1.0l);
      else
	setV(i,j,0.0l);
    }
  }
}

void swapDouble(double & a, double & b){
  double v = a;
  a = b;
  b = v;
}

// simple Gauss-Jordan with full pivoting
// see Numerical Recipes in C p. 39

void Matrix::inverse() 
{
  assert(d1 == d2);
  int d = d1;
  std::vector<int> rowPivot(d);
  std::vector<int> colPivot(d);
  std::vector<int> ipivot(d);
  for (int i = 0; i < d; i++){
    ipivot[i]=0;
  }
  for (int i = 0; i < d; i++){
    //		cout << "iter " << i << endl;
    // main loop over columns to be reduced
    // find pivot element
    double max = 0.0;
    int rowIndex = -1;
    int colIndex = -1;
    for (int j = 0; j < d; j++){
      if (ipivot[j]!=1){
	for (int k = 0; k <d; k++){
	  if (ipivot[k]==0){
	    double vv = fabs(getV(j,k));
	    if (vv >= max){
	      max = vv;
	      rowIndex = j;
	      colIndex = k;
	    }
	  }
	  else if (ipivot[k] > 1)
	    {
	      cout << "Error:: Gaussian singular matrix." << endl;
	      exit(-1);
	    }
	}
      }
    }
    // now we have the pivot element
    ++(ipivot[colIndex]);
    //		cout << "pivot point has val " << vv(rowIndex, colIndex) << endl;
    if (rowIndex != colIndex)
      {
	for (int l = 0; l < d; l++){
	  //	cout << vv(rowIndex,l) << endl;
	  swapDouble(vv(rowIndex,l), vv(colIndex,l));
	  //	cout << vv(rowIndex,l) << endl;
	}
      }
    rowPivot[i] = rowIndex;
    colPivot[i] = colIndex;
    //		cout << "pivot point has val " << vv(colIndex, colIndex) << endl;
    if (getV(colIndex,colIndex) == 0.0l){
      cout << "Singular matrix" << endl;
      exit(-1);
    }
    double inverse = 1.0l/getV(colIndex,colIndex);
    vv(colIndex,colIndex) = 1.0l;
    for (int l = 0 ; l<d; l++){
      vv(colIndex,l) *= inverse;
    }
    for (int ll = 0; ll < d; ll++){
      if (ll != colIndex){
	double dum= vv(ll,colIndex);
	vv(ll,colIndex) = 0.0l;
	for (int l = 0; l <d; l++){
	  vv(ll,l) -= vv(colIndex,l) * dum;
	}
      }
    }
  }
  // so we unscramble
  for (int l = d-1; l>=0; l--){
    if (rowPivot[l]!=colPivot[l]){
      for (int k = 0; k< d;k++){
	swapDouble(vv(k,rowPivot[l]), vv(k,colPivot[l]));
      }
    }
  }
}

void Matrix::times(const Matrix& a, const Matrix& b)
{
  assert(a.d2 == b.d1);
  assert(d1 == a.d1);
  assert(d2 == b.d2);
  int mid = a.d2;
  
  for (int x = 0; x < d1; x++){
    for (int y = 0; y < d2; y++){
      // calculate x,y of the result
      double acc = 0.0l;
      for (int z = 0; z < mid; z++){
	acc += a(x,z) * b(z,y);
      }
      setV(x,y,acc);
    }
  }
}

void Matrix::normaliseRowColumn()
{
  // cache the row and column sums so it is n^2 rather than n^3
  vector<double> rowSums(d1);
  vector<double> colSums(d2);
  for (int row = 0; row < d1; row++)
    rowSums[row] = rowSum(row);
  for (int col = 0; col < d2; col++)
    colSums[col] = columnSum(col);
  // do it
  for (int row = 0; row < d1; row++){
    for (int col = 0; col < d2; col++){
      double vv = getV(row,col);
      if (vv>0)
	setV(row,col,vv/(rowSums[row] + colSums[row]));
    }
  }
  // done
}


double Matrix::mi() const
{
  assert(d1 == d2);
  double s = sum();
  if (s == 0)
    return 0.0l;
  vector<double> before(d1,0.0l);
  vector<double> after(d1,0.0l);
  for (int j = 0; j < d1; j++){
    for (int k = 0; k < d1; k++){
      double v = getV(j,k);
      before[j] += v;
      after[k] += v;
    }
  }
  // so now we have p and q
  double answer = 0.0l;
  for (int x = 0; x < d1; x++){
    for (int y = 0; y < d2; y++){
      double p = getV(x,y) / s;
      double q = before[x] * after[y] / (s * s);
      if (p > 0)
	answer += p * log(p/q);
    }
  }
  return answer;

}


//
//
// Stuff for Matrix of integers
//
//

// copy constructor

MatrixInt::MatrixInt(const MatrixInt & other)
{
  d1 = other.d1;
  d2 = other.d2;
  int size = d1 * d2;
  v = new int[size];
  for (int i = 0; i < size; i++)
    v[i] = other.v[i];
}
  

MatrixInt::MatrixInt(int x, int y)
{
  //  cout << "Allocating MatrixInt (" << x << "," << y << ")" << endl;
  assert(x > 0);
  assert(y > 0);
  d1 = x;
  d2 = y;
  int size = x * y;
  v = new int[size];
  for (int i = 0; i < size; i++)
    v[i] = 0;
}

MatrixInt::~MatrixInt()
{
  delete[] v;
}





void MatrixInt::dump() const
{
  for (int i = 0; i < d1; i++){
    for (int j = 0; j < d2; j++){
      cout << "(" << i << "," << j << ")=" << getV(i,j) << " ";
    }
    cout << endl;
  }
}

void MatrixInt::dumpNonZero() const
{
  for (int i = 0; i < d1; i++){
    for (int j = 0; j < d2; j++){
      if (getV(i,j)>0)
	cout << "(" << i << "," << j << ")=" << getV(i,j) << " ";
    }
    cout << endl;
  }
}

bool MatrixInt::isPermutation() const
{
  if (d1 != d2)
    return false;
  for (int i = 0; i < d1; i++){
    for (int j = 0; j < d2; j++){
      int value = getV(i,j);
      if (value > 0 && value != 1)
	return false;
    }
  }
  return true;
}

bool MatrixInt::isIdentity() const
{
  if (d1 != d2)
    return false;
  for (int i = 0; i < d1; i++){
    for (int j = 0; j < d2; j++){
      int value = getV(i,j);
      if (i==j && value!=1)
	return false;
      if (i!=j && value!=0)
	return false;
    }
  }
  return true;
}

int MatrixInt::rowSum(int x) const
{
  assert(x >= 0 && x < d1);
  int sum = 0;
  for (int y = 0; y < d2; y++){
    sum += getV(x,y);
  }
  return sum;
}	

int MatrixInt::columnSum(int y) const
{
  assert(y >= 0 && y < d2);
  int sum = 0;
  for (int x = 0; x < d1; x++){
    sum += getV(x,y);
  }
  return sum;
}

const MatrixInt& MatrixInt::operator=(const MatrixInt& other)
{
  assert(d1 == other.d1);
  assert(d2 == other.d2);
  int n = d1*d2;
  for (int i = 0;i< n;i++)
    v[i]= other.v[i];
  return *this;
}

const MatrixInt& MatrixInt::operator=(int newValue)
{
  int n = d1*d2;
  for (int i = 0;i< n;i++)
    v[i]= newValue;
  return *this;
}


const MatrixInt& MatrixInt::operator*=(int newValue)
{
  int n = d1*d2;
  for (int i = 0;i< n;i++)
    v[i] *= newValue;
  return *this;
}

int MatrixInt::operator ()(int x, int y) const
{

  assert(x >= 0 && x < d1);
  assert(y >= 0 && y < d2);
  return v[(x * d2) + y];	
}

int& MatrixInt::operator ()(int x, int y)
{
  assert(x >= 0 && x < d1);
  assert(y >= 0 && y < d2);
  return v[(x * d2) + y];	
}

// return a crude upper bound of the permanent

