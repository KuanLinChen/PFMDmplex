#pragma once
using namespace std;
inline double DotProduct( double *A, double *B )
{
	return A[0]*B[0] + A[1]*B[1] + A[2]*B[2] ;
};

inline double DotProduct( double *A, double *B, int ndim )
{
	double value=0.0;
	for(int i = 0 ; i < ndim ; i++ )
	{
		value += A[i]*B[i] ;
	}
	return value ;
};
// double ** MatAlloc( int row, int col )
// {
// 	double ** ptr = new double*[row];
// 	for( int i=0 ; i<row ; i++)
// 	{
// 		ptr[i] = new double[col];
// 	}

// 	for (int i=0 ; i < row ; i++ ){
// 		for (int j=0 ; j < col ; j++ ){
// 			ptr[i][j]=0.0;
// 		}
// 	}
// 	return ptr;
// }
// void MatFree ( double** ptr, int row, int col)
// {
// 	for(int i = 0; i < row; i++)
// 	{
// 		delete [] ptr[i];
// 	}
// 	delete [] ptr;
// }