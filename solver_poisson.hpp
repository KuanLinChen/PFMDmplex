#pragma once
#include "common.hpp" 
#include <petsc.h>
#include "linear_solvers_structure.hpp"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
using namespace std;
class CPoisson 
{
  
private:
  
public:
  CPoisson();
  CDomain *m ;
  CVariable *var ;
  CSysSolve *s ;
  Vec Gradient[3];
  void Init( CDomain *, CVariable * ) ;
  void Solve();
  //void Solver(CVariable * var);
  void Bulid_A_B();
  void CalculateGraditntLSQ() ;

};


