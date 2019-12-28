#include "solver_poisson.hpp"
CPoisson::CPoisson() 
{
  
}
void CPoisson::Init( CGeometry *mm, CVariable *vv )
{
  s = new CSysSolve() ;

  m = mm ;
  var = vv ;
  s->Init( &dmCell, &m->gVec, m->cEndInterior) ;

}
void CPoisson::Solve()
{
	Bulid_A_B();

	s->Solve();

	s->ViewMatrix();
	/* update */
	DMGlobalToLocalBegin( dmCell, s->solution2, INSERT_VALUES, var->Potential ) ;
	DMGlobalToLocalEnd  ( dmCell, s->solution2, INSERT_VALUES, var->Potential ) ;


	PetscScalar *sol ;
	VecGetArray( var->Potential, &sol ) ;
	for( int i = m->cStart ; i < m->cEndInterior ; i++ )
	{
		PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank[%d] val[%d]: %f \n",mpi_rank, i, sol[i] ) ;
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	VecRestoreArray( var->Potential, &sol ) ;


	CalculateGraditntLSQ() ;


	PetscViewer viewer ;
	PetscViewerCreate(PetscObjectComm((PetscObject)dmCell), &viewer);
	PetscViewerSetType( viewer, PETSCVIEWERVTK) ;
	PetscViewerFileSetName( viewer, "flow.vtk");

	VecView( var->Potential, viewer);

	VecView( var->ElectricField[0], viewer);
	VecView( var->ElectricField[1], viewer);

	//VecView( var->Debug, viewer);
	PetscViewerDestroy(&viewer);
}
void CPoisson::Bulid_A_B()
{
	s->MatAZeroEntries() ;
	VecZeroEntries( s->B ) ;
	int P, N ; double Ad_dPN ;

  for ( unsigned int i = 0 ; i < m->interior_face_loop.size() ; i++ ) {
		Ad_dPN = m->interior_face_loop[i]->face_nrml_mag / m->interior_face_loop[i]->dPN ;
		//C0
		P = m->interior_face_loop[i]->cgeom[0]->gindex ;
		N = m->interior_face_loop[i]->cgeom[1]->gindex ;
		s->Add_Entries( P, P, -Ad_dPN ) ;
		s->Add_Entries( P, N,  Ad_dPN ) ;
		//C1
  	P = m->interior_face_loop[i]->cgeom[1]->gindex ;
  	N = m->interior_face_loop[i]->cgeom[0]->gindex ;
  	s->Add_Entries( P, P, -Ad_dPN ) ;
		s->Add_Entries( P, N,  Ad_dPN ) ;
	}

	for ( unsigned int i = 0 ; i < m->processors_face_loop.size() ; i++ ) {
		Ad_dPN= m->processors_face_loop[i]->face_nrml_mag / m->processors_face_loop[i]->dPN ;
		//C0
		P = m->processors_face_loop[i]->cgeom[0]->gindex ;
		N = m->processors_face_loop[i]->cgeom[1]->gindex ;
		s->Add_Entries( P, P, -Ad_dPN ) ;
		s->Add_Entries( P, N,  Ad_dPN ) ;
	}

	for ( unsigned int i = 0 ; i < m->ground_face_loop.size() ; i++ ) {
		Ad_dPN =  m->ground_face_loop[i]->face_nrml_mag /  m->ground_face_loop[i]->dPN ;
   	P =  m->ground_face_loop[i]->cgeom[0]->gindex ;
   	s->Add_Entries( P, P, -Ad_dPN ) ;
	}

	double BCValue = 10.0 ;
	for ( unsigned int i = 0 ; i < m->power_face_loop.size() ; i++ ) {
		Ad_dPN=  m->power_face_loop[i]->face_nrml_mag /  m->power_face_loop[i]->dPN ;
  	P =  m->power_face_loop[i]->cgeom[0]->gindex ;
  	var->FaceField[ m->power_face_loop[i]->offsets-m->fStart ] = BCValue ;
  	s->Add_Entries( P, P, -Ad_dPN ) ;
		s->add_source ( P, -BCValue*Ad_dPN ) ;
	}


	MatAssemblyBegin( s->A , MAT_FINAL_ASSEMBLY ) ;
	MatAssemblyEnd  ( s->A , MAT_FINAL_ASSEMBLY ) ;

	VecAssemblyBegin( s->B ) ;
	VecAssemblyEnd  ( s->B ) ;

}
void CPoisson::CalculateGraditntLSQ()
{
	PetscScalar *Gx, *Gy, *Gz, *value ;
	VecZeroEntries( var->Field[0] ) ;
	VecZeroEntries( var->Field[1] ) ;
	//VecZeroEntries( var->Field[2] ) ;

	VecGetArray(var->Field[0], &Gx) ;
	VecGetArray(var->Field[1], &Gy) ;
	//VecGetArray(var->Field[2], &Gz) ;

	VecGetArray(var->Potential, &value) ;

	double dVar=0.0 ;
	int ids=0 ;

	for ( unsigned int i = 0 ; i < m->neumann_face_loop.size() ; i++ ) {
		ids = m->neumann_face_loop[i]->offsets - m->fStart ;
  	var->FaceField[ ids ] = value[ m->neumann_face_loop[i]->cgeom[0]->index ];
	}

	for ( int i = m->cStart ; i < m->cEndInterior ; i++ ) {

		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ ) {

			ids = m->cell_all[i].nghbr_cell[k] ;
			dVar =  value[ ids ]- value[ i ];
			Gx[ i ] = Gx[ i ] +  m->cell_all[i].nghbr_cell_Cx[ k ]*dVar ;
			Gy[ i ] = Gy[ i ] +  m->cell_all[i].nghbr_cell_Cy[ k ]*dVar ;
		}
		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_face.size() ; k++ ) {

			ids = m->cell_all[i].nghbr_face[ k ]-m->fStart ;
			dVar =  var->FaceField[ids] - value[ i ];
			Gx[ i ] = Gx[ i ] +  m->cell_all[i].nghbr_face_Cx[ k ]*dVar ;
			Gy[ i ] = Gy[ i ] +  m->cell_all[i].nghbr_face_Cy[ k ]*dVar ;

		}		
	}//End cell loop

	VecRestoreArray(var->Field[0], &Gx) ;
	VecRestoreArray(var->Field[1], &Gy) ;
	VecRestoreArray(var->Field[2], &Gz) ;
	VecRestoreArray(var->Potential, &value) ;

	/* update */
	DMGlobalToLocalBegin( dmCell, var->Field[0], INSERT_VALUES, var->ElectricField[0] ) ;
	DMGlobalToLocalEnd  ( dmCell, var->Field[0], INSERT_VALUES, var->ElectricField[0] ) ;

	DMGlobalToLocalBegin( dmCell, var->Field[1], INSERT_VALUES, var->ElectricField[1] ) ;
	DMGlobalToLocalEnd  ( dmCell, var->Field[1], INSERT_VALUES, var->ElectricField[1] ) ;
//DMGlobalToLocalBegin( dmCell, var->Field[2], INSERT_VALUES, var->ElectricField[2] ) ;
//DMGlobalToLocalEnd  ( dmCell, var->Field[2], INSERT_VALUES, var->ElectricField[2] ) ;
	
}

