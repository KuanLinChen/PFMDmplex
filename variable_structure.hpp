#pragma once
#include "common.hpp" 
#include <petsc.h>
#include "linear_solvers_structure.hpp"
#include "domain_structure.hpp"
using namespace std;
// class Scalar
// {
// private:
  
// public:
// 	PetscScalar *data ;
// 	Vec vector ;

// };
class CVariable 
{
  
private:
  
public:
  CVariable();
  CDomain *m ;
  void Init( CDomain *m ) ;
  //Vec LSQ[3];

  void compute_LSQ_coefficient() ;

  void allocate_variable_vectors() ;
  double *FaceField ;//local
  Vec Field[3] ;
  Vec Potential, ElectricField[3] ;
  Vec Debug ;

};


