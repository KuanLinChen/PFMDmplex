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
	CCell(){};
	vector<int> nghbr_cell ;
	vector<int> nghbr_face ;
	vector<double> nghbr_cell_Cx;
	vector<double> nghbr_cell_Cy;
	vector<double> nghbr_cell_Cz;

	vector<double> nghbr_face_Cx;
	vector<double> nghbr_face_Cy;
	vector<double> nghbr_face_Cz;
	int nnghbrs_cell()
	{
		return nghbr_cell.size();
	}
	int nnghbrs_face()
	{
		return nghbr_cell.size();
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
			c0 = -9999 ;
			c1 = -9999 ;
			for (int i = 0 ; i < 3 ; i++ )      centroid[i] = 0.0 ;
			for (int i = 0 ; i < 3 ; i++ )     face_nrml[i] = 0.0 ;
			face_nrml_mag = 0.0 ;
		};
		int offsets ;
		//PetscInt face_type ; //1: Interior face, 2: boundary face, 3: processors boundary face.
		int nnghbrs ; 		/*!< \brief number of neighbors */ 

		int  c0,  c1 ;		/*!< \brief left and right cell index. (local ) */
		int gc0, gc1 ;		/*!< \brief left and right cell index. (global) */

		double c0_centroid[3], c1_centroid[3] ; //cell centroid point.
		double   c0_volume, c1_volume   ;
		double dL, dR ;/*!< \brief left and right cell centroid to face centroid index. */

		double c0_auxiliary[3], c1_auxiliary[3] ; 
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
		double face_nrml[3], face_nrml_mag, dA, centroid[3], nA[3] ;
		void calculate_normal_mag( PetscInt iface )
		{
			Vec            coordinates;
			PetscScalar   *coords = NULL;
			PetscInt       coordSize ;
			PetscSection   coordSection;

			DMGetCoordinatesLocal ( dmMesh, &coordinates  ) ;
			DMGetCoordinateSection( dmMesh, &coordSection ) ;

			DMPlexVecGetClosure( dmMesh, coordSection, coordinates, iface, &coordSize, &coords);
			PetscScalar x1 = coords[0] ;
			PetscScalar y1 = coords[1] ;
			PetscScalar x2 = coords[2] ;
			PetscScalar y2 = coords[3] ;
			face_nrml_mag = sqrt( pow( x1-x2, 2.0 ) + pow( -(y1-y2), 2.0 ) );
			DMPlexVecRestoreClosure( dmMesh, coordSection, coordinates, iface, &coordSize, &coords);
		}
		void calculate_normal_mag2()
		{
			face_nrml_mag = sqrt( pow( face_nrml[0], 2.0 ) + pow( face_nrml[1], 2.0 )+ pow( face_nrml[2], 2.0 ) );
		}
		void calculate_face_area()
		{
			dA = sqrt( pow( face_nrml[0], 2.0 ) + pow( face_nrml[1], 2.0 )+ pow( face_nrml[2], 2.0 ) );
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


		Vec CellGeometryVec, FaceGeometryVec ;
		Vec gVec, lVec ;
		PetscInt global_vector_size, local_vector_size ;
		/* coordinate informations */
			int fStart, fEnd ;
			int eStart, eEnd ;
			int vStart, vEnd ;
			int cStart, //start cell index on each processor, usually is zero 
					cEnd, 	//end cell indel on each processor including the processor boundary ghost cells.
					cEndInterior ;//end cell indel on each processor excluding the processor boundary ghost cells.

		/* Face data */
			map<int,CFace*> face_map_all ;
			map<int,CCell*> cell_map_all ;

			//map<int,CFace*> power_face_loop, ground_face_loop, neumann_face_loop, 
			//interior_face_loop ;//The face data contain only the interior faces, which are shared by two cells.
			vector<CFace*> power_face_loop, ground_face_loop, neumann_face_loop, processors_face_loop, 
			interior_face_loop ;//The face data contain only the interior faces, which are shared by two cells.


			//double *gid ;
			//vector<CFace*> power_face_loop2 ;

			void BulidFaceCellLoopMap() ;

			void ExtractCellGeomInformations() ;
			double *volume, *centroid ;
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

