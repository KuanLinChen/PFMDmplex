#include "geometry_structure.hpp"
using namespace std ;
CGeometry::CGeometry()
{
}
void CGeometry::Init( string filename )
{
	/*--- Read mesh file ---*/
	ReadMeshFromFile( filename ) ;

	/*--- Create cell DM ---*/
	CreateCellDM() ;
	CreateLocalToGlobalCellIdMapping() ; //Also calculate the ghost cell index.
	PetscPrintf( PETSC_COMM_SELF, "rank: %d, cStart: %d, cEnd: %d, cEndInterior: %d\n", mpi_rank, cStart, cEnd, cEndInterior ) ;

	/*--- Computer the cell & face geometry using PETSc default function ---*/
	DMPlexComputeGeometryFVM( dmMesh, &CellGeometryVec, &FaceGeometryVec ) ;

	ExtractCellGeomInformations() ;
	ExtractFaceGeomInformations() ;

	BulidFaceCellLoopMap();

	/* Output mesh */
	Vec partition ;
	VecDuplicate ( gVec, &partition );

	PetscViewer viewer ;
	PetscViewerCreate(PetscObjectComm((PetscObject)dmCell), &viewer);
	PetscViewerSetType( viewer, PETSCVIEWERVTK) ;
	PetscViewerFileSetName( viewer, "Mesh.vtu");
}
void CGeometry::ReadMeshFromFile( string filename )
{
	PetscPrintf(PETSC_COMM_WORLD,"\n---------------------- Read Grid File Information -----------------------\n") ;
	PetscPrintf( PETSC_COMM_WORLD, "The mesh file: %s\n", filename.c_str() ) ;

	/*--- Read mesh file & Create PETSc DMPLEX object.  ---*/
		DMPlexCreateFromFile( PETSC_COMM_WORLD, filename.c_str(), PETSC_TRUE, &dmMesh ) ;

	/* Notes :
    	FEM:   Two points p and q are adjacent if q \in closure(star(p)),   useCone = PETSC_FALSE, useClosure = PETSC_TRUE
    	FVM:   Two points p and q are adjacent if q \in support(p+cone(p)), useCone = PETSC_TRUE,  useClosure = PETSC_FALSE
    	FVM++: Two points p and q are adjacent if q \in star(closure(p)),   useCone = PETSC_TRUE,  useClosure = PETSC_TRUE
	*/
		DMSetBasicAdjacency( dmMesh, PETSC_TRUE, PETSC_TRUE ) ;


	/*--- Distribute mesh over processes ---*/
		DM dmDist ; 
		PetscInt overlap = 1 ;
		DMPlexDistribute( dmMesh, overlap, NULL, &dmDist ) ;
		if ( dmDist ) {//If run one processes, this part will be ingnore.
			PetscPrintf( PETSC_COMM_WORLD, "Distribute mesh over %d processes with %d overlapping !\n", mpi_size, overlap ) ;
			DMDestroy( &dmMesh );
			dmMesh   = dmDist;
		}

		PetscPrintf( PETSC_COMM_WORLD, "%d dimensional problem.\n", GetDimension() ) ;

		DM gdm ; 
		PetscInt ghoscell ;
		DMLabel vtkLabel ;

	/*--- Using 'DMPlexConstructGhostCells' to construct the 'vtkLabel' for output. ---*/
		DMCreateLabel(dmMesh, "dummy") ;
		DMPlexConstructGhostCells(dmMesh, "dummy", &ghoscell, &gdm); 
		DMGetLabel( gdm, "vtk", &vtkLabel ) ;
		DMAddLabel(dmMesh, vtkLabel);
		DMDestroy( &gdm );

		#define monitor_dmMesh false
		#if ( monitor_dmMesh == true )
		PetscPrintf( PETSC_COMM_WORLD, "\n") ;
		DMView( dmMesh, PETSC_VIEWER_STDOUT_WORLD) ;
		PetscPrintf( PETSC_COMM_WORLD, "\n") ;
		#endif


	/*--- Get the range of cell indices ( including cpu boundary overlapping cells ) ---*/
		cStart = cEnd = vStart = vEnd = fStart = fEnd = eStart = eEnd = 0 ;
		DMPlexGetHeightStratum( dmMesh, 0, &cStart, &cEnd ) ; /* cell */
		DMPlexGetDepthStratum ( dmMesh, 0, &vStart, &vEnd ) ; /* verties */
		DMPlexGetHeightStratum( dmMesh, 1, &fStart, &fEnd ) ; /* face */
		if( GetDimension() == 3 ) DMPlexGetDepthStratum ( dmMesh, 1, &eStart, &eEnd ) ; /* edge */
		 

		PetscPrintf( PETSC_COMM_WORLD, "rank: %d, cStart: %d, cEnd: %d\n", mpi_rank, cStart, cEnd ) ;
		PetscPrintf( PETSC_COMM_WORLD, "rank: %d, vStart: %d, vEnd: %d\n", mpi_rank, vStart, vEnd ) ;
		PetscPrintf( PETSC_COMM_WORLD, "rank: %d, fStart: %d, fEnd: %d\n", mpi_rank, fStart, fEnd ) ;
		PetscPrintf( PETSC_COMM_WORLD, "rank: %d, eStart: %d, eEnd: %d\n", mpi_rank, eStart, eEnd ) ;

		ReadBoundaryCellMarkersFromFile( filename.c_str() ) ;
}
void CGeometry::CreateCellDM()
{
  PetscSF        sfPoint;
  PetscSection   coordSection;
  Vec            coordinates, gcoordinates;
  PetscSection   sectionCell;
  PetscInt       Start, End;

	DMClone( dmMesh, &dmCell );

	DMGetCoordinateSection( dmMesh, &coordSection);
	DMSetCoordinateSection( dmCell, PETSC_DETERMINE, coordSection);

	DMGetCoordinates      ( dmMesh,  &gcoordinates) ;
	DMSetCoordinates     ( dmCell,   gcoordinates ) ;

	DMGetCoordinatesLocal ( dmMesh,  &coordinates) ;
	DMSetCoordinatesLocal( dmCell,   coordinates ) ;

	DMGetPointSF( dmMesh, &sfPoint);
	DMSetPointSF( dmCell, sfPoint);

	PetscSectionCreate(PetscObjectComm((PetscObject) dmMesh), &sectionCell);

	PetscSectionSetNumFields(sectionCell, 1) ; PetscSectionSetFieldName(sectionCell, 0, " ");

	DMPlexGetHeightStratum( dmCell, 0, &Start, &End) ;

	PetscSectionSetChart( sectionCell, Start, End) ;

	for ( PetscInt c = Start; c < End; ++c ) 
	{
	  PetscSectionSetDof( sectionCell, c, 1 ) ;
	  PetscSectionSetFieldDof( sectionCell, c, 0, 1 ) ;
	}
	PetscSectionSetUp(sectionCell);

	DMSetLocalSection( dmCell, sectionCell);

	PetscPrintf(PETSC_COMM_WORLD, "create cell dm ... done\n") ;

	DMCreateGlobalVector( dmCell, &gVec ) ;
	DMCreateLocalVector ( dmCell, &lVec ) ;
	VecGetSize( gVec, &global_vector_size);
 	VecGetSize( lVec, &local_vector_size);

	PetscSynchronizedPrintf( PETSC_COMM_WORLD, "rank[%d]-> gVec size->%d, lVec size->%d\n", mpi_rank, global_vector_size, local_vector_size ) ;
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

}
void CGeometry::CreateLocalToGlobalCellIdMapping()
{
		
	/* Note: This gets a borrowed reference, so the user should not destroy this PetscSection. */
		DMGetGlobalSection( dmCell, &global_ids_section ) ;

 		PetscInt GhostCellCount=0, gid=0;	

		for ( PetscInt i = cStart ; i < cEnd ; i++ ) {
			//If the negative number is n, the corresponding global number is -(n+1).
			gid = GetGlobalCellId( i ); 
			if ( gid < 0 ) {
				GhostCellCount++ ;
				mapLocal2GlobalIds[ i ] = -(gid+1) ;
			} else {
				VecSetValue( gVec, gid, gid, INSERT_VALUES ) ;
				mapLocal2GlobalIds[i] = gid ;
			}
		}
		VecAssemblyBegin(gVec); VecAssemblyEnd  (gVec);

		l2g = new int [ cEnd-cStart ] ;
		for ( PetscInt i = cStart ; i < cEnd ; i++ ) 
		{
			l2g[ i ] = Local2GlobalIds(i) ; 
		}

		#define print_Local2GlobalMapping false
		#if( print_Local2GlobalMapping == true)
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank[%d]\n",mpi_rank) ;
			for ( int i = cStart ; i < cEnd ; i++ ) {
				gid = Local2GlobalIds(i) ; 
				PetscSynchronizedPrintf( PETSC_COMM_WORLD,"lid-> %d \t gid-> %d, offsets->%d\n",i, gid, gid-i) ;
			}
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
			PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
		#endif

		cEndInterior = cEnd - GhostCellCount ;

		PetscSynchronizedPrintf( PETSC_COMM_WORLD, "rank[%d]-> number of ghost cell %d\n", mpi_rank, GhostCellCount ) ;
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);


		DMGlobalToLocalBegin( dmCell, gVec, INSERT_VALUES, lVec ) ;
		DMGlobalToLocalEnd( dmCell, gVec, INSERT_VALUES, lVec ) ;

	/* calculate the interior face index */
		PetscInt cell_nnghbrs_face ;
		const PetscInt *cell_nghbr_face ;

		for ( PetscInt i = cStart ; i < cEndInterior ; i++ ) {
			DMPlexGetConeSize( dmMesh, i, &cell_nnghbrs_face) ;
			DMPlexGetCone    ( dmMesh, i, &cell_nghbr_face  ) ;
			for (PetscInt j = 0 ; j < cell_nnghbrs_face ; j++ ){
				interior_face_ids.push_back(cell_nghbr_face[j]);
			}
		}

		#define print_interior_face_ids false
		#if( print_interior_face_ids == true)
			PetscSynchronizedPrintf( PETSC_COMM_WORLD, "mpi_rank-> %d, total interior face (before): %d\n", mpi_rank, interior_face_ids.size() ) ;
		#endif

			sort( interior_face_ids.begin(), interior_face_ids.end() ) ;

			vector<int>::iterator it;
			it = std::unique( interior_face_ids.begin(), interior_face_ids.end() ) ;
			interior_face_ids.resize( distance(interior_face_ids.begin(),it) ); 

		#if( print_interior_face_ids == true)
				PetscSynchronizedPrintf( PETSC_COMM_WORLD, "total interior face (after): %d\n",  interior_face_ids.size() ) ;
				for ( unsigned int i = 0 ; i < interior_face_ids.size() ; i++ ) {
					PetscSynchronizedPrintf( PETSC_COMM_WORLD, "iface -> %d\n", interior_face_ids[i] ) ;
				}
			PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
		#endif
}

void CGeometry::ExtractCellGeomInformations()
{

	/*--- Create Cell-Cell relaction ---*/
		DMPlexCreateNeighborCSR(dmMesh, 0, &ncells, &offsets, &adjacency) ;/* Notes: The neighbor cell indexs are the local indexs. */
	
		#define print_neighborCSR false
		#if ( print_neighborCSR == true )
				for ( PetscInt i = cStart ; i < cEnd ; i++ ) {
					PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, Cell[%d]-> nnghbrs: %d\n", mpi_rank, i, cell_cell_nnghbrs(i) ) ;
					for ( PetscInt j = 0 ; j < cell_cell_nnghbrs(i) ; j++ ) {
						PetscSynchronizedPrintf( PETSC_COMM_WORLD,"nghbr[%d]-> %d \n",j, cell_nghbr(i,j) ) ;
					}
					PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
				}
				PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
		#endif 

	/*--- Extract the cell volume and centroid ---*/
		const PetscScalar *cgeom ;
		DM dmCell ;
		VecGetDM( CellGeometryVec, &dmCell);
		VecGetArrayRead(CellGeometryVec, &cgeom);

		cell_all = new CCell [ cEnd-cStart ] ; 

		for( PetscInt i=cStart ; i <cEnd ; i++ ) {
			PetscFVCellGeom *cg ;
			//typedef struct { PetscReal centroid[3]; PetscReal volume; } PetscFVCellGeom;
			DMPlexPointLocalRead(dmCell, i, cgeom, &cg ) ;

			cell_all[ i ].index 			= i ;
			cell_all[ i ].centroid[0] = cg->centroid[0] ;
			cell_all[ i ].centroid[1] = cg->centroid[1] ;
			cell_all[ i ].centroid[2] = cg->centroid[2] ;
			cell_all[ i ].volume      = cg->volume ;
		}
		VecRestoreArrayRead(CellGeometryVec, &cgeom);

	#define print_volume_centroid false
	#if ( print_volume_centroid == true )
		for ( PetscInt i = cStart ; i < cEnd ; i++ ) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, volume[%d]-> %e\n", mpi_rank, i, cell_all[i].volume ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"cx[%d]->%e, cy[%d]->%e, cz[%d]->%e \n", i, cell_all[i].centroid[0], i, cell_all[i].centroid[1], i, cell_all[i].centroid[2] ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif

	/*---  create cell-to-node lists for convenience. ---*/
		cell_node = new vector<int> [ cEnd-cStart ] ;
		int ids ;
		for ( PetscInt i=cStart ; i <cEnd ; i++ ) {

		  PetscInt *points = NULL, numPoints = 0;
			DMPlexGetTransitiveClosure(dmMesh, i, PETSC_TRUE, &numPoints, &points);
			//cout<<"iCell: "<< i <<"\t"<<"numPoints: "<<numPoints<<endl;
			for (int k=0 ; k < numPoints*2 ; k++ ){
				ids = points[k];
				//cout<<points[k]<<endl;
				if( ids >= vStart and ids < vEnd ) {
					cell_node[i].push_back( ids ) ;
				}
			}
			//cout<<endl;
			DMPlexRestoreTransitiveClosure(dmMesh, i, PETSC_FALSE, &numPoints, &points);
		}
		#define monitor_cell_node_list false
		#if ( monitor_cell_node_list == true )
	 	Vec            coordinates;
	 	PetscScalar   *coords = NULL;
	 	PetscInt       coordSize ;
	 	PetscSection   coordSection;

	 	DMGetCoordinatesLocal ( dmMesh, &coordinates  ) ;
	 	DMGetCoordinateSection( dmMesh, &coordSection ) ;

		for (int i = cStart ; i < cEnd ; i++ ) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, i: %d\n", mpi_rank, i) ;
			for ( unsigned int k = 0 ; k < cell_node[i].size() ; k++ ) {
				DMPlexVecGetClosure( dmMesh, coordSection, coordinates, cell_node[i][k], &coordSize, &coords);
				PetscScalar x1 = coords[0] ;
				PetscScalar y1 = coords[1] ;
				PetscScalar z1 = coords[2] ;
				PetscSynchronizedPrintf( PETSC_COMM_WORLD,"node: %d, x->%e, y->%e, z->%e \n", cell_node[i][k], x1, y1, z1 ) ;
				DMPlexVecRestoreClosure( dmMesh, coordSection, coordinates, cell_node[i][k], &coordSize, &coords);
			}
				PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
		#endif


	/*---  create node-to-cell lists for convenience. ---*/
		node_cell = new vector<int> [ vEnd-vStart ] ;
		for ( int i = cStart ; i < cEnd ; i++ ) {
			for ( unsigned int k = 0 ; k < cell_node[i].size() ; k++ ) {
				node_cell[ cell_node[i][k] - vStart ].push_back(i) ;
			}
		}

		#define monitor_node_cell_list true
		#if ( monitor_node_cell_list == true )
		for (int i = vStart ; i < vEnd ; i++ ) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, i: %d\n", mpi_rank, i) ;
			for ( unsigned int k = 0 ; k < node_cell[i-vStart].size() ; k++ ) {
				PetscSynchronizedPrintf( PETSC_COMM_WORLD,"%d\n", node_cell[i-vStart][k]) ;
			}
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
		#endif
//-------------------------------------------------------------




	PetscEnd();

	CreateCellNeighborVector();

}
void CGeometry::CreateCellNeighborVector()
{
	string Type="cell_face" ;
	//string Type="cell_cell" ;
	int cell_nnghbrs_cell, face_nnghbrs_cell ;

	int             cell_nghbr_cell_index ;
	const PetscInt *cell_nghbr_face_index ;
	PetscInt cell_nnghbrs_face ;//, face_nnghbrs_node ;
	//const PetscInt *cell_nghbr_node_list ;

	if ( Type == "cell_face" ) {

		for( PetscInt i=cStart ; i <cEnd ; i++ ) {

			/* cell-cell */
			cell_nnghbrs_cell  = cell_cell_nnghbrs(i) ;

			for ( int k = 0 ; k < cell_nnghbrs_cell ; k++ ) {

				cell_nghbr_cell_index = cell_nghbr(i,k);
				cell_all[ i ].nghbr_cell.push_back(cell_nghbr_cell_index) ;

			}

			/* cell-face */
			DMPlexGetConeSize( dmMesh, i, &cell_nnghbrs_face) ;
			DMPlexGetCone    ( dmMesh, i, &cell_nghbr_face_index ) ;

			for ( int k = 0 ; k < cell_nnghbrs_face ; k++ ) {

				DMPlexGetSupportSize( dmMesh, cell_nghbr_face_index[k], &face_nnghbrs_cell) ;

				if ( face_nnghbrs_cell == 1 )  {
					cell_all[ i ].nghbr_face.push_back(cell_nghbr_face_index[k]) ;
				}

			}
		}//end cell loop

	} else if ( Type == "cell_cell" ) {

		// for( PetscInt i=cStart ; i <cEndInterior ; i++ ) {

		// 	DMPlexGetConeSize( dmMesh, i, &cell_nnghbrs_face) ;
		// 	DMPlexGetCone    ( dmMesh, i, &cell_nghbr_face_index ) ;	

		// 	for ( int k=0 ; k < cell_nghbr_face ; k++ ) {

		// 		DMPlexGetConeSize( dmMesh, cell_nghbr_face_index[k], &face_nnghbrs_node ) ;
		// 		cout<<"face nnghbrs node: "<<face_nnghbrs_node<<endl;
		// 		DMPlexGetCone    ( dmMesh, i, &cell_nghbr_node_list ) ;
		// 		for (int l=0 ; l < face_nnghbrs_node ; l++ ) {
		// 			cout<<"list["<<k<<"]: "<<cell_nghbr_node_list[k]<<endl;
		// 		}
		// 	}
		// }

		// for ( PetscInt i=cStart ; i <cEnd ; i++ ) {
	 //    PetscInt *points = NULL, numPoints, p, dof, cldof = 0;

		// 	DMPlexGetTransitiveClosure(dmMesh, i, PETSC_FALSE, &numPoints, &points);
		// 	cout<<"i: "<< i <<endl;
		// 	for (int k=0 ; k < numPoints ; k++ ){
		// 		cout<<points[k]<<endl;
		// 	}cout<<endl;
		// 	DMPlexRestoreTransitiveClosure(dmMesh, i, PETSC_FALSE, &numPoints, &points);
		// }



		// for ( PetscInt i=fStart ; i < fEnd ; i++ ) {
		// 	DMPlexGetConeSize(dmMesh, i, &face_nnghbrs_node ) ;
		// 	cout<<"face nnghbrs node: "<<face_nnghbrs_node<<endl;
		// 	DMPlexGetCone    ( dmMesh, i, &cell_nghbr_node_list ) ;
		// 	for (int k=0 ; k < face_nnghbrs_node ; k++ ) {
		// 		cout<<"list["<<k<<"]: "<<cell_nghbr_node_list[k]<<endl;
		// 	}
		// }


	}//End cell_cell
	#define print_cell_information false
	#if ( print_cell_information == true )
		for ( PetscInt i = cStart ; i < cEnd ; i++ ) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, nnghbrs cell: %d, nnghbrs face: %d\n", mpi_rank, cell_all[i].nnghbrs_cell(), cell_all[i].nnghbrs_face()) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
		//PetscEnd();
	#endif 
}
void CGeometry::ExtractFaceGeomInformations()
{

	face_all = new CFace[ fEnd-fStart ] ;

	const PetscInt    *adjacent_cell ;
	const PetscScalar *fgeom ;

	DM dmFace ;
	VecGetDM( FaceGeometryVec, &dmFace);
	VecGetArrayRead(FaceGeometryVec, &fgeom);

	for ( PetscInt i = fStart ; i < fEnd ; i++ ) {

		face_all[i-fStart].offsets = i ;
		/* Number of adjacent cells at the face. */
		face_all[i-fStart].nnghbrs_cell = face_nnghbrs_cell(i) ;

		DMPlexGetSupport( dmMesh, i, &adjacent_cell) ;

		PetscFVFaceGeom *fg ;
		DMPlexPointLocalRead(dmFace, i, fgeom, &fg ) ;

		face_all[i-fStart].face_nrml[0] = fg->  normal[0] ;
		face_all[i-fStart].face_nrml[1] = fg->  normal[1] ;
		face_all[i-fStart].face_nrml[2] = fg->  normal[2] ;

		face_all[i-fStart]. centroid[0] = fg->centroid[0] ;
		face_all[i-fStart]. centroid[1] = fg->centroid[1] ;
		face_all[i-fStart]. centroid[2] = fg->centroid[2] ;

		face_all[i-fStart].calculate_normal_mag();


		if ( face_all[i-fStart].nnghbrs_cell == 1 ) {//means that the face is boundary face or processor boundary face.

			face_all[i-fStart].cgeom[0] = &cell_all[ adjacent_cell[0] ] ;
			face_all[i-fStart].cgeom[1] = &cell_all[ adjacent_cell[0] ] ;

			face_all[i-fStart].cgeom[0]->gindex = Local2GlobalIds(  adjacent_cell[0] )  ;
			face_all[i-fStart].cgeom[0]->gindex = Local2GlobalIds(  adjacent_cell[0] )  ;

		} else {


			face_all[i-fStart].cgeom[0] = &cell_all[ adjacent_cell[0] ] ;
			face_all[i-fStart].cgeom[1] = &cell_all[ adjacent_cell[1] ] ;

			face_all[i-fStart].cgeom[0]->gindex = Local2GlobalIds(  adjacent_cell[0] )  ;
			face_all[i-fStart].cgeom[1]->gindex = Local2GlobalIds(  adjacent_cell[1] )  ;

		}
	}//End face loop

	#define print_faceNormal  false
	#if (print_faceNormal == true)
		for ( PetscInt i =fStart ; i < fEnd ; i++ ) {

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, face[%d]-> nnghbrs: %d\n", mpi_rank, i, face_all[i-fStart].nnghbrs_cell ) ;

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"nghbr[0]-> %d \n", face_all[i-fStart].cgeom[0]->index ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"nghbr[1]-> %d \n", face_all[i-fStart].cgeom[1]->index ) ;

			for (PetscInt d = 0 ; d < GetDimension() ; d++ )
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"face_nrml[%d]-> %e \n",d, face_all[i-fStart].face_nrml[d] ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"face_nrml_mag-> %e \n",face_all[i-fStart].face_nrml_mag) ;

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif

	ComputeFaceCellInformations();
}
void CGeometry::ComputeFaceCellInformations()
{
	// //int C0, C1 ;

	 for ( PetscInt i = fStart ; i < fEnd ; i++ ) {

	 	if( face_all[i-fStart].nnghbrs_cell == 2 ) {

			for ( int k = 0 ; k < 3 ; k++ ) {
				face_all[i-fStart].PN[ k ] = face_all[ i-fStart ].cgeom[1]->centroid[k] - face_all[i-fStart].cgeom[0]->centroid[k] ;
			}

			if ( DotProduct( face_all[i-fStart].PN, face_all[i-fStart].face_nrml ) < 0.0 )
			{
				PetscPrintf( PETSC_COMM_SELF, "flip face\n" ) ;
				face_all[i-fStart].exchange() ;
				PetscEnd() ;
			}

			for ( int k = 0 ; k < 3 ; k++ ) {

				/*--- PN Vector ---*/
				face_all[i-fStart].PN[ k ] = face_all[i-fStart].cgeom[1]->centroid[k] - face_all[i-fStart].cgeom[0]->centroid[k] ;

				/*--- Pf Vector ---*/
				face_all[i-fStart].Pf[ k ] = face_all[i-fStart].centroid[k] - face_all[i-fStart].cgeom[0]->centroid[k] ;

				/*--- Nf Vector ---*/
				face_all[i-fStart].Pf[ k ] = face_all[i-fStart].centroid[k] - face_all[i-fStart].cgeom[1]->centroid[k] ;

				/*--- P'f Vector ---*/
				face_all[i-fStart].PPf[ k ] = DotProduct( face_all[i-fStart].Pf, face_all[i-fStart].face_nrml )*face_all[i-fStart].face_nrml[k] ;

				/*--- N'f Vector ---*/
				face_all[i-fStart].NPf[ k ] = DotProduct( face_all[i-fStart].Nf, face_all[i-fStart].face_nrml )*face_all[i-fStart].face_nrml[k] ;

				/*--- PP'  Vector ( PP' = Pf - P'f ) ---*/
				face_all[i-fStart].PPP[ k ] = face_all[i-fStart].Pf[ k ] - face_all[i-fStart].PPf[ k ] ;

				/*--- NN'  Normal Vector ( NN' = Nf - N'f ) ---*/
				face_all[i-fStart].NNP[ k ] = face_all[i-fStart].Nf[ k ] - face_all[i-fStart].NPf[ k ] ;

			}//end x,y,z

			face_all[i-fStart].dPN     = sqrt( DotProduct( face_all[i-fStart].PN, face_all[i-fStart].PN   ) ) ; 
			face_all[i-fStart].dPPf    = sqrt( DotProduct( face_all[i-fStart].PPf, face_all[i-fStart].PPf ) ) ; 
			face_all[i-fStart].dNPf    = sqrt( DotProduct( face_all[i-fStart].NPf, face_all[i-fStart].NPf ) ) ; 

		
		} else if( face_all[i-fStart].nnghbrs_cell == 1) {/* For boundary */

			for ( int k = 0 ; k < 3 ; k++ ) {

				face_all[i-fStart].PN[ k ] = face_all[i-fStart].centroid[k] - face_all[i-fStart].cgeom[0]->centroid[k] ;

				/*--- Pf Vector ---*/
				face_all[i-fStart].Pf[ k ] = face_all[i-fStart].PN[ k ] ;

				/*--- Nf Vector ---*/
				face_all[i-fStart].Nf[ k ] = 0.0 ;

				/*--- P'f Vector ---*/
				face_all[i-fStart].PPf[ k ] = DotProduct( face_all[i-fStart].Pf, face_all[i-fStart].face_nrml )*face_all[i-fStart].face_nrml[k] ;

				/*--- N'f Vector ---*/
				face_all[i-fStart].NPf[ k ] = 0.0 ;

				/*--- PP'  Vector ( PP' = Pf - P'f ) ---*/
				face_all[i-fStart].PPP[ k ] = face_all[i-fStart].Pf[ k ] - face_all[i-fStart].PPf[ k ] ;

				/*--- NN'  Normal Vector ( NN' = Nf - N'f ) ---*/
				face_all[i-fStart].NNP[ k ] = 0.0 ;

			}//end x,y,z
			face_all[i-fStart].dPPf    = sqrt( DotProduct( face_all[i-fStart].PPf, face_all[i-fStart].PPf ) ) ; 
			face_all[i-fStart].dPN     = face_all[i-fStart].dPPf ;
			face_all[i-fStart].dNPf    = 0.0 ; 

		} else {
			PetscPrintf( PETSC_COMM_SELF, "ERROR @ ComputeFaceCellInformations\n" ) ;
			PetscEnd();
		}
	}
	PetscBarrier(NULL);
}
void CGeometry::construct_cell_face()
{
		/*--- extract cell & face properties (Label) form DMPLEX ---*/
		DMGetLabel( dmMesh, "Cell Sets", &cell_label ) ;
		DMGetLabel( dmMesh, "Face Sets", &face_label ) ;

		/*Get number of cell types in the dmplex */
		DMGetLabelSize( dmMesh, "Cell Sets", &cell_label_size ) ;

		/*Get number of face types in the dmplex */
		DMGetLabelSize( dmMesh, "Face Sets", &face_label_size ) ;

		/*Get cell &face label index sets (IS) */
		DMLabelGetValueIS( cell_label, &cell_IS ) ;
		DMLabelGetValueIS( face_label, &face_IS ) ;

		#define print_markers false
		#if (print_markers == true)
			PetscPrintf(PETSC_COMM_WORLD,"\n***************************\n");
			ISView( cell_IS, PETSC_VIEWER_STDOUT_WORLD );
			ISView( face_IS, PETSC_VIEWER_STDOUT_WORLD );
		#endif
}
void CGeometry::ReadBoundaryCellMarkersFromFile(string filename)
{
	/*--- Get face & cell label informations. ---*/
		DMGetLabel    ( dmMesh, "Cell Sets", &cell_label      ) ;
		DMGetLabelSize( dmMesh, "Cell Sets", &cell_label_size ) ;

		DMGetLabel    ( dmMesh, "Face Sets", &face_label      ) ;
		DMGetLabelSize( dmMesh, "Face Sets", &face_label_size ) ;

	/*--- Get face & cell index sets. ---*/
		DMLabelGetValueIS( cell_label, &cell_IS ) ;
		DMLabelGetValueIS( face_label, &face_IS ) ;

	#define monitor_markers false
	#if (monitor_markers == true)
		PetscPrintf(PETSC_COMM_WORLD,"\n***************************\n");
		ISView( cell_IS, PETSC_VIEWER_STDOUT_WORLD );
		ISView( face_IS, PETSC_VIEWER_STDOUT_WORLD );
	#endif
	PetscPrintf(PETSC_COMM_WORLD,"\n---------------------- Mesh boundary & cell markers ---------------------\n") ;
						/* example
						$PhysicalNames
						5
						1 2 "AAA"
						1 3 "AAA"
						1 4 "CCC"
						1 5 "DDD"
						2 1 "plasma"
						$EndPhysicalNames
						 */
	/*--- Read mesh markers. ---*/
		fstream file ;
		string  line ;

		file.open( filename, ios::in ) ; 
		if (!file) {
			printf(" Fail to open mesh file to read the markers" ) ;
			PetscEnd();
		}
		int PhysicalNames = 0, type, index ; string name, name2;

	  while( getline( file, line ) ){
	     if(line == "$PhysicalNames")
	     	break ;
	  }
	  file >> PhysicalNames ;

	  for (int i = 0 ; i < PhysicalNames ; i++ )
	  {
	  	file >> type >> index >> name ;
	  	name2= remove_chars( name, "\"" ) ;

	  	PhysNames[name] = index ;
	  }
		file.close() ; 	file.clear() ;

		PetscPrintf( PETSC_COMM_WORLD,"%d surface & cell markers.\n", PhysNames.size()) ;
		PetscPrintf( PETSC_COMM_WORLD,"+---------------------+\n") ;
		PetscPrintf( PETSC_COMM_WORLD,"|ids|   Physical names|\n") ;
		PetscPrintf( PETSC_COMM_WORLD,"+---------------------+\n") ;
		for(auto it = PhysNames.cbegin(); it != PhysNames.cend(); ++it) {
			PetscPrintf( PETSC_COMM_WORLD,"|%3d|%17s|\n", it->second, it->first.c_str() ) ;
		}
		PetscPrintf( PETSC_COMM_WORLD,"+---------------------+\n") ;
		//PetscPrintf( PETSC_COMM_WORLD,"\n" ) ;

		IS all_face_IS ;
		DMLabelGetValueIS( face_label, &all_face_IS ) ;
		PetscInt location ;

		for(auto it = PhysNames.cbegin(); it != PhysNames.cend(); ++it) {
			ISLocate( all_face_IS, PhysNames[ it->first ], &location) ;
			if ( location < 0 ) {
				//NOT FOUND.
			} else {
				PhysNamesCPU[it->first] =  it->second ;
			}
		}//Loop of physical names within the mesh file.

		#define monitor_markers_map_CPU false
		#if ( monitor_markers_map_CPU == true )
				PetscPrintf( PETSC_COMM_WORLD,"PhysNamesCPU\n") ;
			//PetscPrintf( PETSC_COMM_WORLD,"mpi_rank: %d\n", mpi_rank ) ;
			for(auto it = PhysNamesCPU.cbegin(); it != PhysNamesCPU.cend(); ++it) {
				PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, name: %s, ids: %d\n", mpi_rank, it->first.c_str(), it->second ) ;
			}
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n" ) ;
			PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT ) ;
		#endif

}
 void CGeometry::BulidFaceCellLoopMap()
{
	IS index_stratum_IS ;
	const PetscInt *ids ;
	PetscInt index_stratum_size ;
	int face_id=0 ;

	/*--- power electroded boundary face id ---*/
		if( PhysNamesCPU.find("POWER") == PhysNamesCPU.end() ) {
			//PetscPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, POWER NOT FOUND\n",mpi_rank) ;
		} else {
			DMLabelGetStratumIS  ( face_label, PhysNamesCPU["POWER"], &index_stratum_IS   ) ;
			DMLabelGetStratumSize( face_label, PhysNamesCPU["POWER"], &index_stratum_size ) ;
			ISGetIndices( index_stratum_IS, &ids) ;
			
			for ( int i = 0 ; i < index_stratum_size ; i++ ) {
				face_id = ids[ i ] - fStart ;
					if ( face_all[face_id].cgeom[0]->index < cEndInterior ) {
						power_face_loop.push_back( &face_all[ face_id ] ) ;
					}
			}
	
			ISRestoreIndices(index_stratum_IS, &ids) ;

			#define print_power_markers_map false
			#if (print_power_markers_map == true)
				// PetscPrintf( PETSC_COMM_SELF,"mpi_rank: %d, POWER FOUND\n",mpi_rank) ;
				for ( unsigned int i = 0 ; i < power_face_loop.size() ; i++ ) {
		  	  cout << power_face_loop[i]->offsets << "\n" ;
				}
			#endif
		}//end if power.

	/*--- Ground electroded boundary face id ---*/
		if ( PhysNamesCPU.find("GROUND") == PhysNamesCPU.end() ) {
			//PetscPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, GROUND NOT FOUND\n",mpi_rank) ;
		} else {

			/* power electroded boundary face id */
			DMLabelGetStratumIS  ( face_label, PhysNamesCPU["GROUND"], &index_stratum_IS   ) ;
			DMLabelGetStratumSize( face_label, PhysNamesCPU["GROUND"], &index_stratum_size ) ;
			ISGetIndices( index_stratum_IS, &ids) ;

			face_id=0 ;
			for ( int i = 0 ; i < index_stratum_size ; i++ ) {
				face_id = ids[ i ] - fStart ;
				if ( face_all[face_id].cgeom[0]->index < cEndInterior ) {
					ground_face_loop.push_back( &face_all[ face_id ] ) ;
				}
			}
			ISRestoreIndices(index_stratum_IS, &ids) ;

			#define print_ground_markers_map false
			#if (print_ground_markers_map == true)
				PetscPrintf( PETSC_COMM_SELF,"mpi_rank: %d, GROUND FOUND\n",mpi_rank) ;
				for ( unsigned int i = 0 ; i < ground_face_loop.size() ; i++ ) {
		  	  cout << ground_face_loop[i]->offsets << "\n" ;
				}
			#endif
		}//end if power.

	/*--- Neumann boundary face id ---*/
		if ( PhysNamesCPU.find("NEUMANN") == PhysNamesCPU.end()  ) {
		//	PetscPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, NEUMANN NOT FOUND\n",mpi_rank) ;
		} else {

			/* power electroded boundary face id */
			DMLabelGetStratumIS  ( face_label, PhysNamesCPU["NEUMANN"], &index_stratum_IS   ) ;
			DMLabelGetStratumSize( face_label, PhysNamesCPU["NEUMANN"], &index_stratum_size ) ;
			ISGetIndices( index_stratum_IS, &ids) ;

			face_id=0 ;
			for ( int i = 0 ; i < index_stratum_size ; i++ ) {
				face_id = ids[ i ] - fStart ;
				if ( face_all[face_id].cgeom[0]->index < cEndInterior ) {
					neumann_face_loop.push_back( &face_all[ face_id ] ) ;
				}
			}
			ISRestoreIndices(index_stratum_IS, &ids) ;

			#define print_neumann_markers_map false
			#if (print_neumann_markers_map == true)
				PetscPrintf( PETSC_COMM_SELF,"mpi_rank: %d, NEUMANN FOUND\n",mpi_rank) ;
				for ( unsigned int i = 0 ; i < neumann_face_loop.size() ; i++ ) {
		  	  cout << neumann_face_loop[i]->offsets << "\n" ;
				}
			#endif
		}//end if power.


	/*--- interior face ---*/
		int m = 0 ;
		for ( unsigned int i = 0 ; i < interior_face_ids.size() ; i++ ){
			m = interior_face_ids[i] - fStart  ;
			if ( face_all[m].nnghbrs_cell == 2 )
			{
				if ( face_all[m].cgeom[1]->index >= cEndInterior ) {
					processors_face_loop.push_back( &face_all[ m ] ) ;
				}else	if ( face_all[m].cgeom[0]->index >= cEndInterior ) {
					processors_face_loop.push_back( &face_all[ m ] ) ;
				} else {
					interior_face_loop.push_back( &face_all[ m ] ) ;
				}
			}
		}

		#define print_interior_markers_map false
		#if (print_interior_markers_map == true)
			PetscPrintf( PETSC_COMM_SELF,"mpi_rank: %d, INTERIOR FOUND\n",mpi_rank) ;
			for ( unsigned int i = 0 ; i < interior_face_loop.size() ; i++ ) {
		    cout << interior_face_loop[i]->offsets << "\n" ;
			}
		#endif

		#define print_processors_markers_map false
		#if (print_processors_markers_map == true)
			PetscPrintf( PETSC_COMM_SELF,"mpi_rank: %d, processors FOUND\n",mpi_rank) ;
			for ( unsigned int i = 0 ; i < processors_face_loop.size() ; i++ ) {
		    cout << processors_face_loop[i]->offsets << "\n" ;
			}
		#endif

		int check_face = interior_face_ids.size() - power_face_loop.size() - ground_face_loop.size()-neumann_face_loop.size()-interior_face_loop.size() - processors_face_loop.size();
		if ( check_face != 0 ) {
			PetscPrintf( PETSC_COMM_SELF,"%d,Loop map error!!!!\n", check_face ) ;
			PetscEnd() ;
		}
		PetscPrintf( PETSC_COMM_WORLD,"BulidFaceCellLoopMap done ...\n",mpi_rank) ;

}
void CGeometry::ViewDMLabels()
{
	  DMLabel        label;
	  const char    *labelName;
	  PetscInt       numLabels ;//, l;

	  /* query the number and name of labels*/
	  DMGetNumLabels( dm, &numLabels);
	  PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD, "Number of labels: %d\n", numLabels);

	  for ( int i = 0 ; i < numLabels ; i++ ) {

	    DMGetLabelName( dm, i, &labelName ) ;
	    PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD, "Label %d: name: %s\n", i, labelName);
	    PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD, "\n");

	    DMGetLabel(dm, labelName, &label);
			DMLabelView( label, PETSC_VIEWER_STDOUT_WORLD ) ;
	  }

}
void CGeometry::ViewDMLabelsIndex()
{
	const PetscInt *ids, *ISid ;
	PetscInt LabelSize ;
	DMLabel vtkLabel ;

	DMGetLabel    ( dmMesh, "vtk",  &vtkLabel ) ;
	DMGetLabelSize( dmMesh, "vtk", &LabelSize ) ;
	PetscSynchronizedPrintf( PETSC_COMM_WORLD,"vtk label size=%d\n",LabelSize ) ;

	IS vtkIS;
	DMLabelGetValueIS( vtkLabel, &vtkIS ) ;

	IS index_stratum_IS ;
	PetscInt index_stratum_size ;

	ISGetIndices( vtkIS, &ISid) ;
	for (int i = 0 ; i < LabelSize ; i++ ) 
	{
		PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank[%d]->IS id: %d \n",mpi_rank, ISid[i]) ;

		DMLabelGetStratumIS  ( vtkLabel, ISid[i], &index_stratum_IS   ) ;
		DMLabelGetStratumSize( vtkLabel, ISid[i], &index_stratum_size ) ;
		ISGetIndices( index_stratum_IS, &ids) ;		

		for ( int k = 0 ; k < index_stratum_size ; k++ ) 
		{
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"id-> %d \n",ids[k]) ;
		}
		ISRestoreIndices(index_stratum_IS, &ids) ;

	}
	ISRestoreIndices( vtkIS, &ISid) ;
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}