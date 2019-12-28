#include "linear_solvers_structure.hpp"
#include <iostream> 
CSysSolve::CSysSolve() 
{
  
}
void CSysSolve::Init( DM *dmplex, Vec *dupVec, PetscInt local_data_size )
{

	ldata = local_data_size ;

	/*1: Create Vector */
	PetscInt low, high, global_data_size ;
	VecCreate( PETSC_COMM_WORLD, &solution ) ;

	VecSetSizes( solution, local_data_size, PETSC_DECIDE ) ; 
	VecSetFromOptions( solution ) ;
	VecDuplicate ( solution, &B ) ;

	VecGetSize ( solution, &global_data_size) ;
	VecGetOwnershipRange( solution, &low, &high) ;

	dia_nz = 20;
	off_nz = 20;

	d_nnz =	new int  [ local_data_size ] ;
	o_nnz =	new int  [ local_data_size ] ;

	for ( int i = 0 ; i < local_data_size ; i++ ) {
	 	d_nnz[ i ] = 20 ;
	 	o_nnz[ i ] = 20 ;
	}
	/*2: Create Matrix*/

	//cout<<"local_data_size: "<<local_data_size<<endl;
	//cout<<"global_data_size: "<<global_data_size<<endl;
	MatCreateAIJ( PETSC_COMM_WORLD, local_data_size, local_data_size, global_data_size, global_data_size, dia_nz, d_nnz, off_nz, o_nnz, &A )  ;

  MatZeroEntries( A ) ;

  VecDuplicate ( *dupVec, &solution2 ) ;
  //DMCreateGlobalVector( dmMesh, &B ) ;
  //VecZeroEntries( B ) ;
  PetscObjectSetName((PetscObject)B,"source term");

  //VecDuplicate ( *dupVec, &solution ) ;
  //DMCreateGlobalVector( dmMesh, &solution ) ;
  PetscObjectSetName((PetscObject)solution,"solution");

	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetOperators( ksp, A, A ) ;

	KSPGetPC( ksp , &pc );

	PCSetType ( pc , PCASM ) ;

	PCASMSetOverlap( pc, 1 ) ;

	PCASMSetType ( pc, PC_ASM_BASIC ) ;

	KSPSetType( ksp, KSPBCGS );

	KSPSetTolerances( ksp, 1.E-8, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );

	KSPSetFromOptions( ksp ) ;

}
void CSysSolve::Solve()
{

	KSPSolve( ksp, B, solution ) ;

	PetscScalar *sol, *sol2 ;
	VecGetArray( solution, &sol ) ;
	VecGetArray( solution2, &sol2 ) ;
	for( int i = 0 ; i < ldata ; i++ )
	{
		sol2[i]=sol[i];
		PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank[%d] val[%d]: %f \n",mpi_rank, i, sol[i] ) ;
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	VecRestoreArray( solution, &sol ) ;
	VecRestoreArray( solution2, &sol2 ) ;
}
int CSysSolve::Add_Entries( int row, int column, double v )
{

	MatSetValues( A, 1, &row, 1, &column, &v,  ADD_VALUES ) ;
	return 0 ;
}
int CSysSolve::add_source( int row, double v )
{
	VecSetValue ( B, row, v, ADD_VALUES  ) ;
	return 0 ;
}
int CSysSolve::MatAZeroEntries()
{
	MatZeroEntries( A ) ;
	return 0 ;
}
int CSysSolve::VecBZeroEntries()
{
	VecZeroEntries( B ) ;
	return 0 ;
}
void CSysSolve::ViewMatrix()
{
	MatView( A,	PETSC_VIEWER_STDOUT_WORLD	) ;
}
void CSysSolve::ViewVector()
{
	VecView( B,	PETSC_VIEWER_STDOUT_WORLD	) ;
}