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

	CreateCellConnectivity() ;
	CreateFaceConnectivity() ;
	CreateNodeConnectivity() ;

	/*--- Computer the cell & face geometry using PETSc default function ---*/
	ExtractCellFaceInformations() ;
	ExtractNodeInformations() ;
	ComputeFaceCellInformations();




	BulidFaceCellLoopMap();

	/* Output mesh */
	Vec partition ;
	VecDuplicate ( gVec, &partition );	
	cout<<"A"<<endl;
	PetscViewer viewer ;
	PetscViewerCreate(PetscObjectComm((PetscObject)dmCell), &viewer);
	PetscViewerSetType( viewer, PETSCVIEWERVTK) ;
	PetscViewerFileSetName( viewer, "Mesh2.vtk");
	VecView( partition, viewer);
	PetscViewerDestroy(&viewer);
	PetscEnd() ;
}
void CGeometry::ReadMeshFromFile( string filename )
{
	PetscPrintf(PETSC_COMM_WORLD,"\n---------------------- Read Grid File Information -----------------------\n") ;
	PetscPrintf( PETSC_COMM_WORLD, "The mesh file: %s\n", filename.c_str() ) ;

	/*--- Read mesh file & Create PETSc DMPLEX object.  ---*/
		DMPlexCreateFromFile( PETSC_COMM_WORLD, filename.c_str(), PETSC_TRUE, &dmMesh ) ;
		//DMPlexCreateFluentFromFile( PETSC_COMM_WORLD, filename.c_str(), PETSC_TRUE, &dmMesh ) ;
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
		int dimension ;
		DMGetDimension( dmMesh, &dimension ) ;
		cStart = cEnd = vStart = vEnd = fStart = fEnd = eStart = eEnd = 0 ;

		DMPlexGetHeightStratum( dmMesh, 0, &cStart, &cEnd ) ; /* cell */
		DMPlexGetDepthStratum ( dmMesh, 0, &vStart, &vEnd ) ; /* verties */
		DMPlexGetHeightStratum( dmMesh, 1, &fStart, &fEnd ) ; /* face */
		if( dimension == 3 ) DMPlexGetDepthStratum ( dmMesh, 1, &eStart, &eEnd ) ; /* edge */
		 

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
void CGeometry::CreateCellConnectivity()
{
	/*--- create mesh connectivity  ---*/
		int dimension ;
		DMGetDimension( dmMesh, &dimension ) ;

		PetscInt *off, *adj, ncells ; 
		DMPlexCreateNeighborCSR(dmMesh, 0, &ncells, &off, &adj) ;

		cell = new CCell[ ncells ] ;
		for (int i=0; i < ncells ; i++ ){
			cell[i].init( dimension ) ;
		}

	/* create cell-cell list */
		int local_id=0 ;
		for ( int i=0 ; i < ncells ; i++ ) {
			cell[i].adj_cell_list_size = off[i+1] - off[ i ] ;
			for ( int j=0 ; j < cell[i].adj_cell_list_size ; j++ ) {

				local_id = adj[j+off[i]] ;
				cell[i].adj_cell_list.push_back( local_id ) ;

			}//j-loop
		}//i-loop

	/* create cell-node list */
		for ( PetscInt i=0 ; i < ncells ; i++ ) {

		  PetscInt *points = NULL, numPoints = 0;
			DMPlexGetTransitiveClosure( dmMesh, i, PETSC_TRUE, &numPoints, &points);
			//cout<<"iCell: "<< i <<"\t"<<"numPoints: "<<numPoints<<endl;

			for ( int k=0 ; k < numPoints*2 ; k+=2 ) {

				local_id = points[k];
				//Node
				if( local_id >= vStart and local_id < vEnd ) {
					cell[i].adj_node_list.push_back( local_id ) ;
				}
				//Face
				if ( local_id >= fStart and local_id < fEnd) {
					cell[i].adj_face_list.push_back( local_id ) ;
				}

			}//k-loop

			cell[i].adj_node_list_size = cell[i].adj_node_list.size();
			cell[i].adj_face_list_size = cell[i].adj_face_list.size();
			DMPlexRestoreTransitiveClosure(dmMesh, i, PETSC_FALSE, &numPoints, &points);
		}//i-loop


	#if ( false )
		for ( PetscInt i = 0 ; i < ncells ; i++ ) {

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, Cell[%d]\n", mpi_rank, i ) ;
			//Cell-Cell
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"adj_cell_list_size = %d\n", cell[i].adj_cell_list_size ) ;
			for ( PetscInt j = 0 ; j < cell[i].adj_cell_list_size ; j++ ) {
				PetscSynchronizedPrintf( PETSC_COMM_WORLD," nghbr[%d]-> %d \n",j, cell[i].adj_cell_list[j] ) ;
			}
			//Cell-Node
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"adj_node_list_size = %d\n", cell[i].adj_node_list_size ) ;
			for ( PetscInt j = 0 ; j < cell[i].adj_node_list_size ; j++ ) {
				PetscSynchronizedPrintf( PETSC_COMM_WORLD," nghbr[%d]-> %d \n",j, cell[i].adj_node_list[j] ) ;
			} 
			//Cell-Face
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"adj_face_list_size = %d\n", cell[i].adj_face_list_size ) ;
			for ( PetscInt j = 0 ; j < cell[i].adj_face_list_size ; j++ ) {
				PetscSynchronizedPrintf( PETSC_COMM_WORLD," nghbr[%d]-> %d \n",j, cell[i].adj_face_list[j] ) ;
			} 
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif 

}
void CGeometry::CreateFaceConnectivity()
{
	int dimension ;
	DMGetDimension( dmMesh, &dimension ) ;

	int nfaces = fEnd-fStart, local_id=0 ;
	face = new CFace[ nfaces ] ;
	for (int i=0; i < nfaces ; i++ ){
		face[i].init( dimension ) ;
	}

	PetscInt nnghbrs ;
	const PetscInt *adjacent_cell ;

	for ( PetscInt i = fStart ; i < fEnd ; i++ ) {

		face[i-fStart].offsets = i ;
		
		/* Number of adjacent cells at the face. */
		DMPlexGetSupportSize( dmMesh, i, &nnghbrs );
		face[i-fStart].adj_cell_list_size = nnghbrs ;
		DMPlexGetSupport( dmMesh, i, &adjacent_cell) ;

		//Cell
		for (int j = 0 ; j < face[i-fStart].adj_cell_list_size ; j++ ) {
			local_id = adjacent_cell[ j ] ;
			face[i-fStart].adj_cell_list.push_back( local_id ) ;
		}

		//Node
		PetscInt *points = NULL, numPoints = 0;
		DMPlexGetTransitiveClosure( dmMesh, i, PETSC_TRUE, &numPoints, &points);
		for ( int k=0 ; k < numPoints*2 ; k+=2 ) {
			local_id = points[k];
			if( local_id >= vStart and local_id < vEnd ) {
				//cout<<points[k]<<endl;
			 	face[i-fStart].adj_node_list.push_back( local_id ) ;
			}
		}//k-loop
		face[i-fStart].adj_node_list_size = face[i-fStart].adj_node_list.size();
		DMPlexRestoreTransitiveClosure(dmMesh, i, PETSC_FALSE, &numPoints, &points);
	}//End face loop

	#if ( false ) 
		for ( PetscInt i = 0 ; i < nfaces ; i++ ) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, face[%d]\n", mpi_rank, i ) ;
			//Cell
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"adj_cell_list_size = %d\n", face[i].adj_cell_list_size ) ;
			for ( int j = 0 ; j < face[i].adj_cell_list_size ; j++ ) {
				PetscSynchronizedPrintf( PETSC_COMM_WORLD," nghbr cell[%d]-> %d \n", j, face[i].adj_cell_list[j] ) ;
			}
			//Node
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"adj_node_list_size = %d\n", face[i].adj_node_list_size ) ;
			for ( int j = 0 ; j < face[i].adj_node_list_size ; j++ ) {
				PetscSynchronizedPrintf( PETSC_COMM_WORLD," nghbr node[%d]-> %d \n", j, face[i].adj_node_list[j] ) ;
			}
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif

}
void CGeometry::CreateNodeConnectivity()
{
	int dimension ;
	DMGetDimension( dmMesh, &dimension ) ;

	int nnode = vEnd-vStart, local_id=0 ;
	node = new CNode[ nnode ] ;
	for( int i = 0 ; i < nnode ; i++ ) {
		node[i].init(dimension);
	}

	for ( PetscInt i = vStart ; i < vEnd ; i++ ) {

		node[i-vStart].offsets = i ;
		
		PetscInt *points = NULL, numPoints = 0;
		DMPlexGetTransitiveClosure( dmMesh, i, PETSC_FALSE, &numPoints, &points);
		for ( int k=0 ; k < numPoints*2 ; k+=2 ) {
			local_id = points[k];
			//cout<<points[k]<<endl;
			if( local_id >= fStart and local_id < fEnd ) {
			 	node[i-vStart].adj_face_list.push_back( local_id ) ;
				//cout<<"face: "<<local_id<<endl;
			}

			if( local_id >= cStart and local_id < cEnd ) {
			 	node[i-vStart].adj_cell_list.push_back( local_id ) ;
				//cout<<"Cell: "<<local_id<<endl;
			}

		}//k-loop
		node[i-vStart].adj_face_list_size = node[i-vStart].adj_face_list.size();
		//cout<<"nFace: "<< node[i-vStart].adj_face_list_size<<endl;
		node[i-vStart].adj_cell_list_size = node[i-vStart].adj_cell_list.size();
		//cout<<"nCell: "<< node[i-vStart].adj_cell_list_size<<endl<<endl;
		DMPlexRestoreTransitiveClosure(dmMesh, i, PETSC_FALSE, &numPoints, &points);
	}//End node loop

	#if ( false )
		for ( PetscInt i = 0 ; i < nnode ; i++ ) {

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, Cell[%d]\n", mpi_rank, i ) ;
			//Node-Cell
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"adj_cell_list_size = %d\n", node[i].adj_cell_list_size ) ;
			for ( PetscInt j = 0 ; j < node[i].adj_cell_list_size ; j++ ) {
				PetscSynchronizedPrintf( PETSC_COMM_WORLD," nghbr[%d]-> %d \n",j, node[i].adj_cell_list[j] ) ;
			}
			//Node-Face
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"adj_face_list_size = %d\n", node[i].adj_face_list_size ) ;
			for ( PetscInt j = 0 ; j < node[i].adj_face_list_size ; j++ ) {
				PetscSynchronizedPrintf( PETSC_COMM_WORLD," nghbr[%d]-> %d \n",j, node[i].adj_face_list[j] ) ;
			} 
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif 
}
void CGeometry::ExtractCellFaceInformations()
{
	Vec cVec, fVec ;
	DM  fdm, cdm;
	DMPlexComputeGeometryFVM( dmMesh, &cVec, &fVec ) ;
	double area=0.0;
	int dimension ;
	DMGetDimension( dmMesh, &dimension ) ;

	const PetscScalar *cgeom, *fgeom ;
	VecGetDM( cVec, &cdm) ;
	VecGetDM( fVec, &fdm) ;

	/* Cell */
	VecGetArrayRead(cVec, &cgeom);
	for ( PetscInt i=cStart ; i <cEnd ; i++ ) {
		PetscFVCellGeom *cg ;
		DMPlexPointLocalRead(cdm, i, cgeom, &cg ) ;

		cell[ i ].index 	= i ;
		cell[ i ].gindex 	= Local2GlobalIds(i) ;

		for ( int j = 0 ; j < dimension ; j++ ) {
			cell[ i ].coords[j] = cg->centroid[j] ;
		}
		cell[ i ].volume      = cg->volume ;
	}
	VecRestoreArrayRead(cVec, &cgeom);

	/* Face */
	VecGetArrayRead(fVec, &fgeom) ;
	for( PetscInt i=fStart ; i <fEnd ; i++ ) {
		area = 0.0 ;
		PetscFVFaceGeom *fg ;
		DMPlexPointLocalRead(fdm, i, fgeom, &fg ) ;

		for ( int j = 0 ; j < dimension ; j++ ){
			face[ i-fStart ].coords[j] = fg->centroid[j] ;
			face[ i-fStart ].normal[j] = fg->normal[j] ;
			area += face[ i-fStart ].normal[j]*face[ i-fStart ].normal[j] ;
		}
		area = sqrt(area) ;
		face[ i-fStart ].area      = area ;
		//normalized
		for( int j = 0 ; j < dimension ; j++ ) {
			face[ i-fStart ].normal[j] =  face[ i-fStart ].normal[j]/area ;
		}
	}
	VecRestoreArrayRead(fVec, &fgeom);

	VecDestroy(&cVec) ;
	VecDestroy(&fVec) ;
	DMDestroy(&cdm) ;
	DMDestroy(&fdm) ;


	//Cell
	#if ( false )
		for ( PetscInt i = cStart ; i < cEnd ; i++ ) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, volume[%d]-> %e\n", mpi_rank, i, cell[i].volume ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD," cx[%d]->%e\n", i, cell[i].coords[0] ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD," cy[%d]->%e\n", i, cell[i].coords[1] ) ;
			if (dimension==3)
			PetscSynchronizedPrintf( PETSC_COMM_WORLD," cz[%d]->%e\n", i, cell[i].coords[2] ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif

	//Face	
	#if ( false )
		for ( PetscInt i = 0 ; i < fEnd-fStart ; i++ ) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, area[%d]-> %e\n", mpi_rank, i, face[i].area ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD," fx[%d]->%e\n", i, face[i].coords[0] ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD," fy[%d]->%e\n", i, face[i].coords[1] ) ;
			if (dimension==3)
			PetscSynchronizedPrintf( PETSC_COMM_WORLD," fz[%d]->%e\n", i, face[i].coords[2] ) ;

			PetscSynchronizedPrintf( PETSC_COMM_WORLD," nx[%d]->%e\n", i, face[i].normal[0] ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD," ny[%d]->%e\n", i, face[i].normal[1] ) ;
			if (dimension==3)
			PetscSynchronizedPrintf( PETSC_COMM_WORLD," nz[%d]->%e\n", i, face[i].normal[2] ) ;

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif
}
void CGeometry::ExtractNodeInformations()
{
	Vec            coordinates;
	PetscScalar   *coords = NULL;
	PetscInt       coordSize ;
	PetscSection   coordSection;
	int dimension ;
	DMGetDimension( dmMesh, &dimension ) ;

	DMGetCoordinatesLocal ( dmMesh, &coordinates  ) ;// This is a borrowed reference, so the user should NOT destroy this vector
	DMGetCoordinateSection( dmMesh, &coordSection ) ;

	for (int i = vStart ; i < vEnd ; i++ ) {
		DMPlexVecGetClosure( dmMesh, coordSection, coordinates, i, &coordSize, &coords);
		for ( int j = 0 ; j < dimension ; j++ ){
			node[i-vStart].coords[j] = coords[j] ; 
		}
		DMPlexVecRestoreClosure( dmMesh, coordSection, coordinates, i, &coordSize, &coords);
	}
}
void CGeometry::CreateCellNeighborVector()
{
	// //string Type="cell_face" ;
	// string Type="cell_cell" ;
	// int cell_nnghbrs_cell, face_nnghbrs_cell ;

	// int             cell_nghbr_cell_index ;
	// const PetscInt *cell_nghbr_face_index ;
	// PetscInt cell_nnghbrs_face ;//, face_nnghbrs_node ;
	// //const PetscInt *cell_nghbr_node_list ;
	// unsigned int iCell, iNode, nCell, nNode ;
	// if ( Type == "cell_face" ) {

	// 	for( PetscInt i=cStart ; i <cEnd ; i++ ) {

	// 		/* cell-cell */
	// 		cell_nnghbrs_cell  = cell_cell_nnghbrs(i) ;

	// 		for ( int k = 0 ; k < cell_nnghbrs_cell ; k++ ) {

	// 			cell_nghbr_cell_index = cell_nghbr(i,k);
	// 			cell_all[ i ].nghbr_cell.push_back(cell_nghbr_cell_index) ;

	// 		}

	// 		/* cell-face */
	// 		DMPlexGetConeSize( dmMesh, i, &cell_nnghbrs_face) ;
	// 		DMPlexGetCone    ( dmMesh, i, &cell_nghbr_face_index ) ;

	// 		for ( int k = 0 ; k < cell_nnghbrs_face ; k++ ) {

	// 			DMPlexGetSupportSize( dmMesh, cell_nghbr_face_index[k], &face_nnghbrs_cell) ;

	// 			if ( face_nnghbrs_cell == 1 )  {
	// 				cell_all[ i ].nghbr_face.push_back(cell_nghbr_face_index[k]) ;
	// 			}

	// 		}
	// 	}//end cell loop

	// } else if ( Type == "cell_cell" ) {

	// 	for( PetscInt i=cStart ; i <cEndInterior ; i++ ) {

	// 		nNode = cell_node[i].size() ;

	// 		for ( unsigned int n=0 ; n < nNode; n++ ) {
	// 			iNode = cell_node[i][n]-vStart ;
	// 			nCell = node_cell[iNode].size() ;
	// 			for (unsigned int c=0 ; c < nCell ; c++ ) {
	// 				iCell = node_cell[iNode][c] ;
	// 				cell_all[ i ].nghbr_cell.push_back(iCell) ;
	// 			}//end node_cell list
	// 		}//end cell_node list

	// 		sort( cell_all[ i ].nghbr_cell.begin(), cell_all[ i ].nghbr_cell.end() ) ;

	// 		vector<int>::iterator it;
	// 		it = std::unique( cell_all[ i ].nghbr_cell.begin(), cell_all[ i ].nghbr_cell.end() ) ;
	// 		cell_all[ i ].nghbr_cell.resize( distance(cell_all[ i ].nghbr_cell.begin(),it) ); 
	// 	}//End cell-loop.


	// }//End cell_cell
	// #define print_cell_information false
	// #if ( print_cell_information == true )
	// 	for ( PetscInt i = cStart ; i < cEnd ; i++ ) {
	// 		PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, nnghbrs cell: %d, nnghbrs face: %d\n", mpi_rank, cell_all[i].nnghbrs_cell(), cell_all[i].nnghbrs_face()) ;
	// 		PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
	// 	}
	// 	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	// 	PetscEnd();
	// #endif 
}

void CGeometry::ComputeFaceCellInformations()
{
	int dimension ;
	DMGetDimension( dmMesh, &dimension ) ;
	int cL=0, cR=0 ;

	for ( PetscInt i = fStart ; i < fEnd ; i++ ) {

	 	switch( face[i-fStart].adj_cell_list_size ) {
	 		case 1:

	 			cL = face[i-fStart].adj_cell_list[0] ;
	 			cR = face[i-fStart].adj_cell_list[0] ;

				for ( int k = 0 ; k < dimension ; k++ ) {

					face[i-fStart].PN[ k ] = face[i-fStart].coords[k] - cell[cL].coords[k] ;
					/*--- Pf Vector ---*/
					face[i-fStart].Pf[ k ] = face[i-fStart].PN[ k ] ;

					/*--- Nf Vector ---*/
					face[i-fStart].Nf[ k ] = 0.0 ;

					/*--- P'f Vector ---*/
					face[i-fStart].PPf[ k ] = DotProduct( face[i-fStart].Pf, face[i-fStart].normal, dimension )*face[i-fStart].normal[k] ;

					/*--- N'f Vector ---*/
					face[i-fStart].NPf[ k ] = 0.0 ;

					/*--- PP'  Vector ( PP' = Pf - P'f ) ---*/
					face[i-fStart].PPP[ k ] = face[i-fStart].Pf[ k ] - face[i-fStart].PPf[ k ] ;

					/*--- NN'  Normal Vector ( NN' = Nf - N'f ) ---*/
					face[i-fStart].NNP[ k ] = 0.0 ;

				}//end x,y,z
				face[i-fStart].dPPf    = sqrt( DotProduct( face[i-fStart].PPf, face[i-fStart].PPf, dimension ) ) ; 
				face[i-fStart].dPN     = face[i-fStart].dPPf ;
				face[i-fStart].dNPf    = 0.0 ; 

	 		break ;

    	case 2:

    		cL = face[i-fStart].adj_cell_list[0] ;
	 			cR = face[i-fStart].adj_cell_list[1] ;

				for ( int k = 0 ; k < dimension ; k++ ) {

					/*--- PN Vector ---*/
					face[i-fStart].PN[ k ] = cell[cR].coords[k] - cell[cL].coords[k] ;

					/*--- Pf Vector ---*/
					face[i-fStart].Pf[ k ] = face[i-fStart].coords[k] - cell[cL].coords[k] ;

					/*--- Nf Vector ---*/
					face[i-fStart].Pf[ k ] = face[i-fStart].coords[k] - cell[cR].coords[k] ;

					/*--- P'f Vector ---*/
					face[i-fStart].PPf[ k ] = DotProduct( face[i-fStart].Pf, face[i-fStart].normal, dimension )*face[i-fStart].normal[k] ;

					/*--- N'f Vector ---*/
					face[i-fStart].NPf[ k ] = DotProduct( face[i-fStart].Nf, face[i-fStart].normal, dimension )*face[i-fStart].normal[k] ;

					/*--- PP'  Vector ( PP' = Pf - P'f ) ---*/
					face[i-fStart].PPP[ k ] = face[i-fStart].Pf[ k ] - face[i-fStart].PPf[ k ] ;

					/*--- NN'  Normal Vector ( NN' = Nf - N'f ) ---*/
					face[i-fStart].NNP[ k ] = face[i-fStart].Nf[ k ] - face[i-fStart].NPf[ k ] ;
				}

				if ( DotProduct( face[i-fStart].PN, face[i-fStart].normal, dimension ) < 0.0 )
				{
					PetscPrintf( PETSC_COMM_SELF, "flip face\n" ) ;
					//face[i-fStart].exchange() ;
					PetscEnd() ;
				}

				face[i-fStart].dPN     = sqrt( DotProduct( face[i-fStart].PN , face[i-fStart].PN , dimension  ) ) ; 
				face[i-fStart].dPPf    = sqrt( DotProduct( face[i-fStart].PPf, face[i-fStart].PPf, dimension  ) ) ; 
				face[i-fStart].dNPf    = sqrt( DotProduct( face[i-fStart].NPf, face[i-fStart].NPf, dimension  ) ) ; 
      break;
			default:
				PetscPrintf( PETSC_COMM_SELF, "ERROR @ ComputeFaceCellInformations\n" ) ;
				PetscEnd();
      break;
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

		#if ( false )
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
					if ( face[face_id].adj_cell_list[0] < cEndInterior ) {
						power_face_loop.push_back( &face[ face_id ] ) ;
					}
			}
	
			ISRestoreIndices(index_stratum_IS, &ids) ;

			#if ( false )
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

			DMLabelGetStratumIS  ( face_label, PhysNamesCPU["GROUND"], &index_stratum_IS   ) ;
			DMLabelGetStratumSize( face_label, PhysNamesCPU["GROUND"], &index_stratum_size ) ;
			ISGetIndices( index_stratum_IS, &ids) ;

			face_id=0 ;
			for ( int i = 0 ; i < index_stratum_size ; i++ ) {
				face_id = ids[ i ] - fStart ;
				if ( face[face_id].adj_cell_list[0] < cEndInterior ) {
					ground_face_loop.push_back( &face[ face_id ] ) ;
				}
			}
			ISRestoreIndices(index_stratum_IS, &ids) ;

			#if ( false )
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
				if ( face[face_id].adj_cell_list[0] < cEndInterior ) {
					neumann_face_loop.push_back( &face[ face_id ] ) ;
				}
			}
			ISRestoreIndices(index_stratum_IS, &ids) ;

			#if ( false )
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
			if ( face[m].adj_cell_list_size == 2 )
			{
				if ( face[m].adj_cell_list[1] >= cEndInterior ) {

					processors_face_loop.push_back( &face[ m ] ) ;

				}else	if ( face[m].adj_cell_list[0] >= cEndInterior ) {

					processors_face_loop.push_back( &face[ m ] ) ;

				} else {

					interior_face_loop.push_back( &face[ m ] ) ;

				}
			}
		}

		#if ( false )
			PetscPrintf( PETSC_COMM_SELF,"mpi_rank: %d, INTERIOR FOUND\n",mpi_rank) ;
			for ( unsigned int i = 0 ; i < interior_face_loop.size() ; i++ ) {
		    cout << interior_face_loop[i]->offsets << "\n" ;
			}
		#endif

		#if ( false )
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