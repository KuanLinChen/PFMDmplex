#pragma once
#include <petscdmplex.h>
#include <petscsys.h>
#include <map>
using namespace std;

extern PetscMPIInt mpi_rank, mpi_size ;
extern DM dmMesh, dmCell ;
extern map<string, int> PhysNames, PhysNamesCPU ;