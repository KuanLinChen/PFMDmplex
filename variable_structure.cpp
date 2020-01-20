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
	double dx, dy ;
	int  ids ;

	int ndim = m->GetDimension(), nnghbrs=0 ;
	QSMatrix<double> A( ndim, ndim, 0.0 ), A_inv( ndim, ndim, 0.0 ), adj( ndim, ndim, 0.0 ) ;


	for ( int i = m->cStart ; i < m->cEnd ; i++ ) {

		nnghbrs = m->cell_all[i].nghbr_cell.size() + m->cell_all[i].nghbr_face.size() ;

		QSMatrix<double> D( nnghbrs, ndim, 0.0 ), Coeff( nnghbrs, ndim, 0.0 ) ;


		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ ) {

			ids = m->cell_all[i].nghbr_cell[k] ;

			dx = m->cell_all[ids].centroid[0]  - m->cell_all[i].centroid[0] ;
			dy = m->cell_all[ids].centroid[1]  - m->cell_all[i].centroid[1] ;

			D( k, 0 ) = dx ; 
			D( k, 1 ) = dy ; 

		}

		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_face.size() ; k++ ) {

			ids = m->cell_all[i].nghbr_face[k] - m->fStart ;

			dx =  m->face_all[ids].centroid[0] - m->cell_all[i].centroid[0] ;
			dy =  m->face_all[ids].centroid[1] - m->cell_all[i].centroid[1] ;

			D( m->cell_all[i].nghbr_cell.size() + k , 0 ) = dx ; 
			D( m->cell_all[i].nghbr_cell.size() + k , 1 ) = dy ; 
		}

		QSMatrix<double> DTran = D.transpose();

		A = DTran*D ;
		adj = A.adjoint();
		A_inv = adj/A.det();

		Coeff = A_inv*DTran ;

		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ )
		{
			m->cell_all[i].nghbr_cell_Cx.push_back( Coeff(0,k) ) ;
			m->cell_all[i].nghbr_cell_Cy.push_back( Coeff(1,k) ) ;
		}

		for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ )
		{
			m->cell_all[i].nghbr_face_Cx.push_back( Coeff(0,k+m->cell_all[i].nghbr_cell.size())) ;
			m->cell_all[i].nghbr_face_Cy.push_back( Coeff(1,k+m->cell_all[i].nghbr_cell.size())) ;
		}

	}//end cell loop.
}
