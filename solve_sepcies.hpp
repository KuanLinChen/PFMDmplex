#pragma once
#include "common.hpp" 
#include <petsc.h>
#include "linear_solvers_structure.hpp"
#include "domain_structure.hpp"
#include "variable_structure.hpp"
using namespace std;
class CSpecies 
{
  
private:
  
public:
  CSpecies();

  CDomain *m ;
  CVariable *var ;
  CSysSolve *s ;




};


