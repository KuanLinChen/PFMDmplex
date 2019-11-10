#pragma once
#include "common.hpp" 
#include <petsc.h>

using namespace std;
class CSysSolve 
{
  
private:
  
public:
  CSysSolve();
  Mat A ;
  Vec B, solution, solution2 ;
 	KSP ksp;
 	PC	pc ;
  int *d_nnz, *o_nnz, dia_nz, off_nz ;
  int ldata ;
  void Init( DM *, Vec *, PetscInt ) ;
  void Solve() ;
  int Add_Entries( int row, int column, double v ) ;
  int add_source ( int row, double v ) ;
  int MatAZeroEntries() ;
  int VecBZeroEntries() ;
  void ViewMatrix();
  void ViewVector();
};


