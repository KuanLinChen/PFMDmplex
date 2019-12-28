#include "variable_structure.hpp"
CVariable::CVariable() 
{
}
void CVariable::Init( CGeometry *mm )
{
	m = mm ;
}
void CVariable::allocate_variable_vectors()
{
	VecDuplicate( m->lVec, &Potential ) ; PetscObjectSetName((PetscObject)Potential,"potential");
	VecDuplicate( m->lVec, &Debug ) ; PetscObjectSetName((PetscObject)Debug,"Debug");

	VecDuplicate( m->lVec, &ElectricField[0] ) ; PetscObjectSetName((PetscObject) ElectricField[0],"Ex");
	VecDuplicate( m->lVec, &ElectricField[1] ) ; PetscObjectSetName((PetscObject) ElectricField[1],"Ey");
	VecDuplicate( m->lVec, &ElectricField[2] ) ; PetscObjectSetName((PetscObject) ElectricField[2],"Ez");
	for ( int i = 0 ; i < 3 ; i++ ){
		VecDuplicate( m->gVec, &Field[i] ) ;
	}
	compute_LSQ_coefficient() ;
	FaceField = new double [ m->fEnd-m->fStart ] ;
}
void CVariable::compute_LSQ_coefficient()
{
	double a11, a12, a21, a22, dx, dy, det, ia11, ia12, ia21, ia22 ;
	int  ids ;
	PetscScalar *val ;
	VecZeroEntries( Debug ) ;

	VecGetArray( Debug, &val) ;
	for ( int i = m->cStart ; i < m->cEnd ; i++ ) {

		/*--- Reset Matrix ---*/
		a11 = 0.0 ; a12 = 0.0 ;
		a21 = 0.0 ; a22 = 0.0 ;

		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ ) {

			ids = m->cell_all[i].nghbr_cell[k] ;

			dx = m->cell_all[ids].centroid[0]  - m->cell_all[i].centroid[0] ;
			dy = m->cell_all[ids].centroid[1]  - m->cell_all[i].centroid[1] ;

			a11 = a11 + dx*dx ; 
			a12 = a12 + dx*dy ;
			a21 = a21 + dx*dy ; 
			a22 = a22 + dy*dy ;
		}

		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_face.size() ; k++ ) {

			ids = m->cell_all[i].nghbr_face[k] - m->fStart ;

			dx =  m->face_all[ids].centroid[0] - m->cell_all[i].centroid[0] ;
			dy =  m->face_all[ids].centroid[1] - m->cell_all[i].centroid[1] ;

			a11 = a11 + dx*dx ; 
			a12 = a12 + dx*dy ;
			a21 = a21 + dx*dy ; 
			a22 = a22 + dy*dy ;
		}
		/*--- Cal. LSQ det. ---*/
		det = a11*a22 - a12*a21 ;

		/*--- Cal. Inverse Matrix ---*/
		ia11 =  a22/det ;
		ia12 = -a21/det ;
		ia21 = -a12/det ;
		ia22 =  a11/det ;

		/*--- Loop over neighbor "cells" ---*/
		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ ) {

			ids = m->cell_all[i].nghbr_cell[k] ;

			dx = m->cell_all[ids].centroid[0]  - m->cell_all[i].centroid[0] ;
			dy = m->cell_all[ids].centroid[1]  - m->cell_all[i].centroid[1] ;

			m->cell_all[i].nghbr_cell_Cx.push_back(ia11*dx + ia12*dy) ;
			m->cell_all[i].nghbr_cell_Cy.push_back(ia21*dx + ia22*dy) ;

			val[i] += 	m->cell_all[i].nghbr_cell_Cx[k] ;
		}
		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_face.size() ; k++ ) {

			ids = m->cell_all[i].nghbr_face[k] - m->fStart ;

			dx =  m->face_all[ids].centroid[0] - m->cell_all[i].centroid[0] ;
			dy =  m->face_all[ids].centroid[1] - m->cell_all[i].centroid[1] ;

			m->cell_all[i].nghbr_face_Cx.push_back(ia11*dx + ia12*dy) ;
			m->cell_all[i].nghbr_face_Cy.push_back(ia21*dx + ia22*dy) ;
			
			val[i] += 	m->cell_all[i].nghbr_face_Cx[k] ;
		}		
	}//end cell loop.
	VecRestoreArray( Debug, &val) ;
}
