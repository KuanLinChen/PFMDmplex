#pragma once
using namespace std;
/*--- System inforimations ---*/
inline void PrintfMemory() 
{
	PetscLogDouble memory ;
	PetscMemoryGetCurrentUsage( &memory) ;
	PetscPrintf( PETSC_COMM_WORLD, "Mesh dimension: %g [MB]\n", memory/(1.0E6) ) ;
}

