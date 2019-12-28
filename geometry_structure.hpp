#pragma once
#include "common.hpp" 
#include "tool_function_cpp.hpp"
#include "tool_function_petsc.hpp"
#include <cmath>

#include "petscdm.h"          
#include "petscdmlabel.h"     
#include "petscviewer.h" 
#include <petsc.h>
#include <petscsf.h>
#include <iostream>
#include <petscviewerhdf5.h>
#include <algorithm>
#include <fstream>
#include <vector>
#include <../src/sys/classes/viewer/impls/vtk/vtkvimpl.h>
using namespace std;
class CCell
{
	public: 
	CCell()
	{
		index = gindex = -999 ;
		centroid[0] = centroid[1] = centroid[3] =  0.0 ;
		volume = 0.0 ;
	};
	vector<int> nghbr_cell ;
	vector<int> nghbr_face ;
	vector<double> nghbr_cell_Cx;
	vector<double> nghbr_cell_Cy;
	vector<double> nghbr_cell_Cz;

	vector<double> nghbr_face_Cx;
	vector<double> nghbr_face_Cy;
	vector<double> nghbr_face_Cz;

	int  index,  /*!< \brief  local cell index. */
	    gindex ; /*!< \brief global cell index. */

	double centroid[3], auxiliary[3] ;
	double volume ;
	int nnghbrs_cell()
	{
		return nghbr_cell.size();
	}
	int nnghbrs_face()
	{
		return nghbr_face.size();
	}
};
/*
!          v2=Left(2)
!        o---o---------o       face(j,:) = [i,k,v2,v1]
!       .    .          .
!      .     .           .
!     .      .normal      .
!    .  Left .--->  Right  .
!   .   c1   .       c2     .
!  .         .               .
! o----------o----------------o
!          v1=Right(1)
*/
class CFace
{
	public:
		CFace()
		{
			face_nrml_mag = 0.0 ;
			face_nrml[0]=0.0;
			face_nrml[1]=0.0;
			face_nrml[2]=0.0;
			 centroid[0]=0.0;
			 centroid[1]=0.0;
			 centroid[2]=0.0;
		} ;

		int offsets ;

		CCell *cgeom[2] ;

		int nnghbrs_cell ; 		/*!< \brief number of neighbors */ 


		double dL, dR ;/*!< \brief left and right cell centroid to face centroid index. */
		//Vector 
		double  PN[3], /*!< \brief Vector of PN. */
						Pf[3], /*!< \brief Vector of Pf. */
						Nf[3], /*!< \brief Vector of Nf. */
					 PPP[3], /*!< \brief Vector of PP'. */
					 NNP[3], /*!< \brief Vector of NN'. */
					 PPf[3], /*!< \brief Vector of PP'. */
					 NPf[3], /*!< \brief Vector of NN'. */
						  dPN, /*!< \brief Distance between P' and f. */
						 dPPf, /*!< \brief Distance between P' and f. */
				 		 dNPf; /*!< \brief Distance between N' and f. */

		//face_nrml_mag(nfaces) = sqrt( ( x1-x2 )**2 + ( -(y1-y2) )**2 )
		double face_nrml_mag ; //, dA, nA[3] ;
		double face_nrml[3] ;
		double  centroid[3] ;
		// void calculate_normal_mag( PetscInt iface )
		// {
		// 	Vec            coordinates;
		// 	PetscScalar   *coords = NULL;
		// 	PetscInt       coordSize ;
		// 	PetscSection   coordSection;

		// 	DMGetCoordinatesLocal ( dmMesh, &coordinates  ) ;
		// 	DMGetCoordinateSection( dmMesh, &coordSection ) ;

		// 	DMPlexVecGetClosure( dmMesh, coordSection, coordinates, iface, &coordSize, &coords);
		// 	PetscScalar x1 = coords[0] ;
		// 	PetscScalar y1 = coords[1] ;
		// 	PetscScalar x2 = coords[2] ;
		// 	PetscScalar y2 = coords[3] ;
		// 	face_nrml_mag = sqrt( pow( x1-x2, 2.0 ) + pow( -(y1-y2), 2.0 ) );
		// 	DMPlexVecRestoreClosure( dmMesh, coordSection, coordinates, iface, &coordSize, &coords);
		// }
		void calculate_normal_mag()
		{
			face_nrml_mag = sqrt( pow( face_nrml[0], 2.0 ) + pow( face_nrml[1], 2.0 )+ pow( face_nrml[2], 2.0 ) );
			face_nrml[0] = face_nrml[0]/face_nrml_mag ;
			face_nrml[1] = face_nrml[1]/face_nrml_mag ;
			face_nrml[2] = face_nrml[2]/face_nrml_mag ;
		}
		void exchange()
		{	
			CCell *tmp0, *tmp1 ;
			tmp0 = cgeom[0] ;
			tmp1 = cgeom[1] ;
			cgeom[0] = tmp1 ;
			cgeom[1] = tmp0 ;
		}
};

class CGeometry 
{
	public:
		CGeometry();

		void Init(string);
		DM dm, var_dm ;

		void ReadMeshFromFile( string ) ;
		void ViewDMLabelsIndex();
		//CCell_Face *cell_face ;


		Vec CellGeometryVec, FaceGeometryVec ;/*!< \brief Petsc vector for store the cell and face geometry. (like cell center, face normal...)*/
		


		Vec gVec, lVec ;
		PetscInt global_vector_size, local_vector_size ;
		/* coordinate informations */
			int fStart, fEnd ; /*!< \brief 'face' start & end index on each processor. */ 
			int eStart, eEnd ; /*!< \brief 'element' start & end index on each processor. (For 3D only.)*/ 
			int vStart, vEnd ; /*!< \brief 'vertices' start & end index on each processor. */ 
			int cStart,        /*!< \brief 'cell' start index on each processor, usually is zero . */
					cEnd,          /*!< \brief 'cell'   end index on each processor including the processor boundary ghost cells. */
					cEndInterior ; /*!< \brief 'cell'   end index on each processor excluding the processor boundary ghost cells. */

		/* Cell & Face data */
			CCell *cell_all ;
			CFace *face_all ;

			vector<CFace*> power_face_loop, ground_face_loop, neumann_face_loop, processors_face_loop, 
			interior_face_loop ;//The face data contain only the interior faces, which are shared by two cells.



			void BulidFaceCellLoopMap() ;
			void ExtractCellGeomInformations() ;
			void CreateCellNeighborVector();
			void ExtractFaceGeomInformations() ;

			vector<int> interior_face_ids, interior_bface_ids ;

			void construct_cell_face() ;

		/*--- Cell data ---*/
			/* Create the cell dm for cell vector alloc. */
			void CreateCellDM() ;


			/* Mapping local_id to global_id. */
			map<int, int> mapLocal2GlobalIds ;
			int *l2g ;

			inline int Local2GlobalIds( int key )
			{
				int iter=0 ;
				if (mapLocal2GlobalIds.find(key) == mapLocal2GlobalIds.end() ) {
					cout<<"ERROR!!!"<<endl; PetscEnd() ;
				} else {
				 iter = mapLocal2GlobalIds[key];	
				 //cout<<iter<<endl;
				}
				 return iter ;
			}
			void CreateLocalToGlobalCellIdMapping() ;

			void ComputeFaceCellInformations() ;
			/* cell-cell adjacency */
			PetscInt *offsets, *adjacency, nnghbrs, ncells ; 
			/* number of cell neighbors */
			inline PetscInt cell_cell_nnghbrs( PetscInt i ) {
				return offsets[i+1]-offsets[ i ] ;
			}
			inline PetscInt cell_face_nnghbrs( PetscInt i ) {
				PetscInt cell_nnghbrs_face ;
				DMPlexGetConeSize(dmMesh, i, &cell_nnghbrs_face) ;
				return cell_nnghbrs_face ;
			}

			/* list of cell neighbors */
			inline PetscInt cell_nghbr( PetscInt i, PetscInt j ) {
				return adjacency[ j+offsets[i] ] ;
			}

		/* To get the global cell index */
			PetscSection global_ids_section ;
			inline int GetGlobalCellId( int i ) { 
				PetscInt global_cell_ids ;
				PetscSectionGetOffset( global_ids_section, i, &global_cell_ids ) ;
				return global_cell_ids ;
			}


		void ReadBoundaryCellMarkersFromFile(string filename) ;

		//DMLabel CellLabels, FaceLabels ;

		DMLabel cell_label, face_label ;
		PetscInt cell_label_size, face_label_size ;
		IS cell_IS, face_IS ;

		DMLabel face_marker_labels, cell_marker_labels ;
		PetscInt nface_markers, ncell_markers ;
		IS cell_marker_IS, face_marker_IS ;

		void ViewDMLabels();

		map<string,int> BC_Marker, Cell_Marker ;

		inline int GetDimension(){ 
			int nDim ; 
			DMGetDimension( dmMesh, &nDim ) ;
			return nDim ;
		}

		inline PetscInt face_nnghbrs_cell( PetscInt iface ){ 
			PetscInt nnghbrs ; 
			DMPlexGetSupportSize(dmMesh, iface, &nnghbrs );
			return nnghbrs ;
		}

		string& remove_chars(string& s, const string& chars) {
		    s.erase(remove_if(s.begin(), s.end(), [&chars](const char& c) {
		        return chars.find(c) != string::npos;
		    }), s.end());
		    return s;
		}
};

