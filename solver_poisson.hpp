#pragma once
#include "common.hpp" 
#include <petsc.h>
#include "linear_solvers_structure.hpp"
#include "geometry_structure.hpp"
#include "variable_structure.hpp"
using namespace std;
class CPoisson 
{
  
private:
  
public:
  CPoisson();
  CGeometry *m ;
  CVariable *var ;
  CSysSolve *s ;
  
  void Init( CGeometry *, CVariable * ) ;
  void Solve();
  //void Solver(CVariable * var);
  void Bulid_A_B();
  void CalculateGraditntLSQ() ;

};


