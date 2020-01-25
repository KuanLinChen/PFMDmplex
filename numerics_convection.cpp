#include "numerics_convection.hpp"
CCentJST::CCentJST() 
{
  
}
void CCentJST::Init( CGeometry *_mesh, CVariable *_var )
{
  m   = _mesh ;
  var =  _var ;

  //Taking from SU2
  Param_p = 0.3;
  Param_Kappa_2 = 0.5 ;
  Param_Kappa_4 = 0.02;
	VecDuplicate( m->lVec, &Und_Lapl ) ; PetscObjectSetName((PetscObject)Und_Lapl,"Und_Lapl");
	VecDuplicate( m->gVec, &Residue ) ; PetscObjectSetName((PetscObject)Residue,"Residue");
}
void CCentJST::Undivided_Laplacian( Vec _var )
{
	PetscScalar *gUnd_Lapl, *value ;
	double dVar=0.0 ;
	int iCell=0, jCell=0 ;

	/* Reset the scalar. */
	VecZeroEntries( var->Scalar ) ;

	VecGetArray(var->Scalar, &gUnd_Lapl) ;
	VecGetArray(_var, &value) ;//tmp setup

  for ( unsigned int i = 0 ; i < m->interior_face_loop.size() ; i++ ) {
  	iCell = m->interior_face_loop[i]->cgeom[0]->index ;
  	jCell = m->interior_face_loop[i]->cgeom[1]->index ;

 		dVar =  value[jCell] - value[iCell] ;

 		gUnd_Lapl[iCell] += dVar ;
 		gUnd_Lapl[jCell] -= dVar ;

	}

	for ( unsigned int i = 0 ; i < m->processors_face_loop.size() ; i++ ) {

  	iCell = m->interior_face_loop[i]->cgeom[0]->index ;
  	jCell = m->interior_face_loop[i]->cgeom[1]->index ;

 		dVar =  value[jCell] - value[iCell] ;

 		gUnd_Lapl[iCell] += dVar ;
 		gUnd_Lapl[jCell] -= dVar ;

	}

	VecRestoreArray(var->Scalar, &gUnd_Lapl) ;
	VecRestoreArray(_var, &value) ;

	/* update */
	DMGlobalToLocalBegin( dmCell, var->Scalar, INSERT_VALUES, Und_Lapl ) ;
	DMGlobalToLocalEnd  ( dmCell, var->Scalar, INSERT_VALUES, Und_Lapl ) ;

}
void CCentJST::ComputeResidual()
{

}

/*---------------------------------------------------*/
CUpwTVD::CUpwTVD() 
{
  
}
void CUpwTVD::Init( CGeometry *_mesh, CVariable *_var, string _name )
{
  m   = _mesh ;
  var =  _var ;
  var_name = _name ;
	VecDuplicate( m->lVec, &Gradient[0] ) ; PetscObjectSetName((PetscObject) Gradient[0], (var_name+"_Gx").c_str() ) ;
	VecDuplicate( m->lVec, &Gradient[1] ) ; PetscObjectSetName((PetscObject) Gradient[1], (var_name+"_Gy").c_str() ) ;
	VecDuplicate( m->lVec, &Gradient[2] ) ; PetscObjectSetName((PetscObject) Gradient[2], (var_name+"_Gz").c_str() ) ;
}
void CUpwTVD::CalculateGraditntLSQ()
{
	// PetscScalar *Gx, *Gy, *Gz, *value ;
	// for ( int i = 0 ; i < m->GetDimension() ; i++ ) {
	// 	VecZeroEntries( var->Field[i] ) ;
	// }

	// VecGetArray(var->Field[0], &Gx) ;
	// VecGetArray(var->Field[1], &Gy) ;
	// VecGetArray(var->Field[2], &Gz) ;
	// VecGetArray(var->Potential, &value) ;

	// double dVar=0.0 ;
	// int ids=0 ;

	// for ( int i = m->cStart ; i < m->cEndInterior ; i++ ) {
	// 	/* Cell */
	// 	for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ ) {
	// 		ids = m->cell_all[i].nghbr_cell[k] ;
	// 		dVar =  value[ ids ]- value[ i ];
	// 		Gx[ i ] = Gx[ i ] +  m->cell_all[i].nghbr_cell_Cx[ k ]*dVar ;
	// 		Gy[ i ] = Gy[ i ] +  m->cell_all[i].nghbr_cell_Cy[ k ]*dVar ;
	// 		Gz[ i ] = Gz[ i ] +  m->cell_all[i].nghbr_cell_Cz[ k ]*dVar ;
	// 	}
	// }//End cell loop

	// VecRestoreArray(var->Field[0], &Gx) ;
	// VecRestoreArray(var->Field[1], &Gy) ;
	// VecRestoreArray(var->Field[2], &Gz) ;
	// VecRestoreArray(var->Potential, &value) ;

	// /* update */
	// for (int i = 0 ; i < m->GetDimension() ; i++ ){
	// 	DMGlobalToLocalBegin( dmCell, var->Field[i], INSERT_VALUES, Gradient[i] ) ;
	// 	DMGlobalToLocalEnd  ( dmCell, var->Field[i], INSERT_VALUES, Gradient[i] ) ;
	// }
}

