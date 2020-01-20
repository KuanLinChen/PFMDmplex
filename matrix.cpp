#ifndef __QS_MATRIX_CPP
#define __QS_MATRIX_CPP

#include "matrix.hpp"

/*--- Constructor ---*/
template<typename T>
QSMatrix<T>::QSMatrix(unsigned _rows, unsigned _cols, const T& _initial) 
{
	mat.resize(_rows);
	for (unsigned i=0; i<mat.size(); i++) {
		mat[i].resize(_cols, _initial);
	}
	rows = _rows;
	cols = _cols;
}

/* Duplicate the matrix */
template<typename T>
QSMatrix<T>::QSMatrix(const QSMatrix<T>& rhs) 
{
  mat = rhs.mat;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

/*--- (Virtual) Destructor ---*/
template<typename T>
QSMatrix<T>::~QSMatrix() {


}

/*--- Assignment Operator ---*/
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator=(const QSMatrix<T>& rhs) 
{
  if (&rhs == this)
    return *this;

  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();

  mat.resize(new_rows);
  for (unsigned i=0; i<mat.size(); i++) {
    mat[i].resize(new_cols);
  }

  for (unsigned i=0; i<new_rows; i++) {
    for (unsigned j=0; j<new_cols; j++) {
      mat[i][j] = rhs(i, j);
    }
  }
  rows = new_rows;
  cols = new_cols;

  return *this;
}

/*--- Addition of two matrices ---*/
template<typename T>
QSMatrix<T> QSMatrix<T>::operator+(const QSMatrix<T>& rhs) 
{
	unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();
  if ( rows != new_rows or cols != new_cols )  {
  	cout<<"Martix size are not matach. (operator+) "<<endl;
  	exit(1) ;
  }

  QSMatrix result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this->mat[i][j] + rhs(i,j);
    }
  }

  return result;
}

// Cumulative addition of this matrix and another
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator+=(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      this->mat[i][j] += rhs(i,j);
    }
  }

  return *this;
}

/*--- Subtraction of this matrix and another ---*/
template<typename T>
QSMatrix<T> QSMatrix<T>::operator-(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();
  QSMatrix result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this->mat[i][j] - rhs(i,j);
    }
  }

  return result;
}

// Cumulative subtraction of this matrix and another                                                                                                                          
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator-=(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      this->mat[i][j] -= rhs(i,j);
    }
  }

  return *this;
}

// Left multiplication of this matrix and another                                                                                                                              
template<typename T>
QSMatrix<T> QSMatrix<T>::operator*(const QSMatrix<T>& rhs) 
{
  unsigned A_rows = rows ;
  unsigned A_cols = cols ;
  unsigned B_rows = rhs.get_rows();
  unsigned B_cols = rhs.get_cols();
  //unsigned rows = rows ;
  if ( A_cols != B_rows ) {
    cout<<"Error! column of first matrix not equal to row of second. "<<endl;
    exit(1) ;
  }
  QSMatrix result( A_rows, B_cols, 0.0);

  for (unsigned i=0; i<A_rows; i++) {
    for (unsigned j=0; j<B_cols; j++) {
      for (unsigned k=0; k<A_cols; k++) {
        result(i,j) += this->mat[i][k] * rhs(k,j);
      }
    }
  }

  return result;
}

// Cumulative left multiplication of this matrix and another                                                                                                                  
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator*=(const QSMatrix<T>& rhs) {
  QSMatrix result = (*this) * rhs;
  (*this) = result;
  return *this;
}

/*--- Calculate a transpose of this matrix ---*/
template<typename T>
QSMatrix<T> QSMatrix<T>::transpose() 
{

  QSMatrix result(cols, rows, 0.0);
  for (unsigned j=0; j<cols; j++) {
    for (unsigned i=0; i<rows; i++) {
      result(j,i) = this->mat[i][j];
    }
  } 
  return result;
}

/* Get minor matrix for matrix calculation. */
template<typename T>
QSMatrix<T> QSMatrix<T>::minor( unsigned p, unsigned q ) 
{
	QSMatrix results( rows-1, cols-1, 0.0);

	unsigned i = 0, j = 0; 
	unsigned n = rows ; 
	//cout<<"p: "<<p<<"\t"<<"q: "<<q<<"\t"<<"n: "<<n<<endl;
	// Looping for each element of the matrix 
	for (unsigned row = 0; row < n; row++) { 
		for (unsigned col = 0; col < n; col++)  { 
			//cout<<"row: "<<row<<"\t"<<"col: "<<col<<"\t"<<this->mat[row][col]<<endl;
			/* Copying into temporary matrix only those element, which are not in given row and column */
			if (row != p && col != q) { 
				results(i, j) = this->mat[row][col]; 
				//cout<<this->mat[row][col]<<endl;
				j++ ;
				// Row is filled, so increase row index and reset col index 
				if (j == n-1 ) { 
					j = 0; 
					i++ ; 
				}
			}
		}//end column 
	}//end row

	return results;
}
/* */
template<typename T>
T QSMatrix<T>::det()
{ 
	unsigned n = rows ;
	T results=0.0 ; 
	if      ( n==1 ) return this->mat[0][0] ;
	else if ( n==2 ) return this->mat[0][0] * this->mat[1][1] -this->mat[1][0] * this->mat[0][1] ;
	QSMatrix C(n-1, n-1, 0.0 ) ; 

	for ( unsigned f = 0; f < n; f++ ) { 
		C = minor( 0, f ) ; 
		results = results + this->mat[0][f]*pow(-1,f)*C.det() ; 
	} 
	return results ; 
} 



/*--- Adjoint ---*/
template<typename T>
QSMatrix<T> QSMatrix<T>::adjoint() 
{
	unsigned n = rows ;

	//if ( n==1 ) return this->mat ;

	QSMatrix   C(n-1, n-1, 0.0 ) ;
	QSMatrix adj(  n,   n, 0.0 ) ; 

	// temp is used to store cofactors of A[][] 
	int sign = 1; 

	for ( int i=0; i<rows; i++) { 
		for ( int j=0; j<cols; j++) { 

			C = minor( i, j ) ; 
			sign = ((i+j)%2==0)? 1: -1; 

			adj(j,i) = (sign)*(C.det()); 

		}//end cols
	}//end rows 
  return adj;
}



// Matrix/scalar addition                                                                                                                                                     
template<typename T>
QSMatrix<T> QSMatrix<T>::operator+(const T& rhs) {
  QSMatrix result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this->mat[i][j] + rhs;
    }
  }

  return result;
}

// Matrix/scalar subtraction                                                                                                                                                  
template<typename T>
QSMatrix<T> QSMatrix<T>::operator-(const T& rhs) {
  QSMatrix result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this->mat[i][j] - rhs;
    }
  }

  return result;
}

// Matrix/scalar multiplication                                                                                                                                               
template<typename T>
QSMatrix<T> QSMatrix<T>::operator*(const T& rhs) {
  QSMatrix result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this->mat[i][j] * rhs;
    }
  }

  return result;
}

// Matrix/scalar division                                                                                                                                                     
template<typename T>
QSMatrix<T> QSMatrix<T>::operator/(const T& rhs) {
  QSMatrix result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this->mat[i][j] / rhs;
    }
  }

  return result;
}

// Multiply a matrix with a vector                                                                                                                                            
template<typename T>
std::vector<T> QSMatrix<T>::operator*(const std::vector<T>& rhs) {
  std::vector<T> result(rhs.size(), 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result[i] = this->mat[i][j] * rhs[j];
    }
  }

  return result;
}

// Obtain a vector of the diagonal elements                                                                                                                                   
template<typename T>
std::vector<T> QSMatrix<T>::diag_vec() {
  std::vector<T> result(rows, 0.0);

  for (unsigned i=0; i<rows; i++) {
    result[i] = this->mat[i][i];
  }

  return result;
}

// Access the individual elements                                                                                                                                             
template<typename T>
T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) {
  return this->mat[row][col];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) const {
  return this->mat[row][col];
}

// Get the number of rows of the matrix                                                                                                                                       
template<typename T>
unsigned QSMatrix<T>::get_rows() const {
  return this->rows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
unsigned QSMatrix<T>::get_cols() const {
  return this->cols;
}

template<typename T>
void  QSMatrix<T>::print() 
{
	for (unsigned i=0; i<rows; i++){
	 for (unsigned j=0; j<cols; j++){
			printf("% 3.4e   ", this->mat[i][j] ) ;
			//cout<<this->mat[i][j]<<"\t" ;
		} 
		cout<<endl;
	} 
	cout<<endl;
}
#endif
