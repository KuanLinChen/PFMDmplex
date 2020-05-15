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
	PetscInt ndim = m->GetDimension() ;
	//Potential
	VecDuplicate( m->lVec, &Potential ) ; 
	PetscObjectSetName((PetscObject)Potential,"potential");
	//
	VecDuplicate( m->lVec, &Debug ) ; 
	PetscObjectSetName((PetscObject)Debug,"Debug");
	//Electric field
	VecDuplicateVecs( m->lVec, ndim, &ElectricField);
	PetscObjectSetName((PetscObject) ElectricField[0],"Ex");
	PetscObjectSetName((PetscObject) ElectricField[1],"Ey");
	if ( ndim ==3 )
	PetscObjectSetName((PetscObject) ElectricField[2],"Ez");

	// Field vector
	VecDuplicateVecs( m->gVec, ndim, &Field);

	VecDuplicate( m->gVec, &Scalar ) ;

	compute_LSQ_coefficient() ;
	FaceField = new double [ m->fEnd-m->fStart ] ;
}
void CVariable::compute_LSQ_coefficient()
{
	// PetscInt ndim = m->GetDimension();

	// double *dDist = NULL;
	// dDist = new double[ndim] ;

	// int  ids=0, nnghbrs=0 ;

	// QSMatrix<double> A( ndim, ndim, 0.0 ), A_inv( ndim, ndim, 0.0 ), adj( ndim, ndim, 0.0 ) ;


	// for ( int i = m->cStart ; i < m->cEnd ; i++ ) {

	// 	nnghbrs = m->cell_all[i].nghbr_cell.size() + m->cell_all[i].nghbr_face.size() ;
	// 	cout<<"nnghbrs: "<<nnghbrs<<endl;
	// 	QSMatrix<double> D( nnghbrs, ndim, 0.0 ), Coeff( nnghbrs, ndim, 0.0 ) ;

	// 	//cell
	// 	for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ ) {

	// 		ids = m->cell_all[i].nghbr_cell[k] ;

	// 		for ( int n=0 ; n < ndim ; n++ ){
	// 			dDist[n] = m->cell[ids].coords[n]  - m->cell_all[i].centroid[n] ;
	// 			D( k, n ) = dDist[n] ; 
	// 		}
	// 	}
	// 	//face
	// 	for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_face.size() ; k++ ) {

	// 		ids = m->cell_all[i].nghbr_face[k] - m->fStart ;

	// 		for ( int n = 0 ; n < ndim ; n++ ) {
	// 			dDist[n] = m->face_all[ids].centroid[n]  - m->cell_all[i].centroid[n] ;
	// 			D( m->cell_all[i].nghbr_cell.size()+k, n ) = dDist[n] ; 
	// 		}
	// 	}

	// 	//calculate the lsq-coefficient.
	// 	QSMatrix<double> DTran = D.transpose();
	// 	A = DTran*D ;
	// 	adj = A.adjoint();
	// 	A_inv = adj/A.det();
	// 	Coeff = A_inv*DTran ;

	// 	//cell
	// 	for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ ) {
	// 		for ( int n = 0 ; n < ndim ; n++ ) {
	// 			m->cell_all[i].lsq_cell[n].push_back( Coeff(n,k) ) ;
	// 		}
	// 	}
	// 	//face
	// 	for ( unsigned int k = 0 ; k < m->cell_all[i].nghbr_cell.size() ; k++ ) {
	// 		for ( int n = 0 ; n < ndim ; n++ ) {
	// 			m->cell_all[i].lsq_face[n].push_back( Coeff(n, k+m->cell_all[i].nghbr_cell.size())) ;
	// 		}
	// 	}

	// }//end cell loop.

	// delete [] dDist;
}
