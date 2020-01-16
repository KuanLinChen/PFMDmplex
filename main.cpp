//#include "common.hpp"
#include "variable_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_poisson.hpp"
#include "solver_sepcies.hpp"
#include "matrix.hpp"

//#include "tool_function_petsc.hpp"

using namespace std ;

PetscMPIInt mpi_rank,mpi_size ;
map<string, int> PhysNames, PhysNamesCPU ;
DM dmMesh, dmCell ;

int main(int argc, char **argv)
{
	PetscInitialize( &argc, &argv, (char *)0, 0) ;

  QSMatrix<double> m(10, 10, 1.0 ) ;

// PetscEnd() ;

	MPI_Comm_size(PETSC_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(PETSC_COMM_WORLD,&mpi_rank);

  CGeometry *mesh ;
  mesh = new CGeometry ;
  //mesh->ReadMeshFromFile( "./mesh/2d_Structured.msh") ;
  mesh->Init("./mesh/2d_Structured.msh") ;

  
  CVariable *variable ;
  variable = new CVariable() ;
  variable->Init( mesh ) ;
  variable->allocate_variable_vectors();
  //cout<<"AAA2"<<endl;
  
	CPoisson *poisson ;
	poisson = new CPoisson();
	poisson->Init( mesh, variable ) ;
	//cout<<"A"<<endl;
	poisson->Solve() ;
	
  PetscEnd();



	//cout<<"B"<<endl;
 //  CSysSolve *poisson ;
 //  poisson = new CSysSolve() ;
 //  poisson->Init( &dmCell, &mesh->gVec, mesh->local_vector_size ) ;

 //  double Ad_dPN = 0.0 ; 
 //  PetscInt P, N ;
 //  PetscScalar value=0 ;

 //  //for( auto it = mesh->interior_face_loop.cbegin(); it != mesh->interior_face_loop.cend(); ++it) {
 //  for ( unsigned int i = 0 ; i < mesh->interior_face_loop.size() ; i++ ) {
	// 	Ad_dPN= mesh->interior_face_loop[i]->dA / mesh->interior_face_loop[i]->dPN ;
	// 	//C0
	// 	P = mesh->interior_face_loop[i]->gc0 ;
	// 	N = mesh->interior_face_loop[i]->gc1 ;
	// 	poisson->Add_Entries( P, P, -Ad_dPN ) ;
	// 	poisson->Add_Entries( P, N,  Ad_dPN ) ;
	// 	//C1
 //  	P = mesh->interior_face_loop[i]->gc1 ;
 //  	N = mesh->interior_face_loop[i]->gc0 ;
 //  	poisson->Add_Entries( P, P, -Ad_dPN ) ;
	// 	poisson->Add_Entries( P, N,  Ad_dPN ) ;
	// }


	// 	//for(auto it = mesh->ground_face_loop.cbegin(); it != mesh->ground_face_loop.cend(); ++it)  {
	// 	for ( unsigned int i = 0 ; i < mesh->ground_face_loop.size() ; i++ ) {
	// 		Ad_dPN =  mesh->ground_face_loop[i]->dA /  mesh->ground_face_loop[i]->dPN ;
 //   		P =  mesh->ground_face_loop[i]->gc0 ;
 //   		poisson->Add_Entries( P, P, -Ad_dPN ) ;
	// 	}

	// 	double BCValue = 10.0 ;
	// 	//for(auto it = mesh->power_face_loop.cbegin(); it != mesh->power_face_loop.cend(); ++it)  {
	// 	for ( unsigned int i = 0 ; i < mesh->power_face_loop.size() ; i++ ) {
	// 		Ad_dPN=  mesh->power_face_loop[i]->dA /  mesh->power_face_loop[i]->dPN ;
 //   		P =  mesh->power_face_loop[i]->gc0 ;

 //   		poisson->Add_Entries( P, P, -Ad_dPN ) ;
	// 		poisson->add_source ( P, -BCValue*Ad_dPN ) ;
	// 	}
	// 	poisson->Solve();

	// PetscViewer viewer ;
	// PetscViewerCreate(PetscObjectComm((PetscObject)dmCell), &viewer);
	// PetscViewerSetType( viewer, PETSCVIEWERVTK) ;
	// PetscViewerFileSetName( viewer, "flow.vtu");
	// VecView( poisson->B, viewer);
	// VecView( poisson->solution, viewer);
	// PetscViewerDestroy(&viewer);
	// test_global(12345 ) ;
//PetscEnd();
	



	PetscFinalize();
	return 0 ;
}
