#pragma once
#include "common.hpp" 
#include <petsc.h>
#include "linear_solvers_structure.hpp"
#include "geometry_structure.hpp"
#include "variable_structure.hpp"
using namespace std;
class CSpecies 
{
  
private:
  
public:
  CSpecies();
  CGeometry *m ;
  CVariable *var ;
  CSysSolve *s ;




};


