#include "domain_structure.hpp"
using namespace std ;
CDomain::CDomain()
{
}
void CDomain::ReadMeshFromFile( string filename )
{
	PetscPrintf( PETSC_COMM_WORLD, "Read the mesh file, the file path: %s\n", filename.c_str() ) ;

	/*--- Create DMPLEX for gmsh file ---*/
		DMPlexCreateFromFile( PETSC_COMM_WORLD, filename.c_str(), PETSC_TRUE, &dmMesh ) ;

	/* Notes :
    	FEM:   Two points p and q are adjacent if q \in closure(star(p)),   useCone = PETSC_FALSE, useClosure = PETSC_TRUE
    	FVM:   Two points p and q are adjacent if q \in support(p+cone(p)), useCone = PETSC_TRUE,  useClosure = PETSC_FALSE
    	FVM++: Two points p and q are adjacent if q \in star(closure(p)),   useCone = PETSC_TRUE,  useClosure = PETSC_TRUE
	*/
		DMSetBasicAdjacency( dmMesh, PETSC_TRUE, PETSC_FALSE ) ;


	/*--- Distribute mesh over processes ---*/
		DM dmDist ; 
		PetscInt overlap = 1 ;
		DMPlexDistribute( dmMesh, overlap, NULL, &dmDist ) ;
		if ( dmDist ) {//If run one processes, this part will be ingnore.
			PetscPrintf( PETSC_COMM_WORLD, "Distribute mesh over %d processes with %d overlapping !\n", mpi_size, overlap ) ;
			DMDestroy( &dmMesh );
			dmMesh   = dmDist;
		}
		PetscPrintf( PETSC_COMM_WORLD, "Mesh dimension: %d\n", GetDimension() ) ;

		DM gdm ; 
		PetscInt ghoscell ;
		DMLabel vtkLabel ;

		DMCreateLabel(dmMesh, "dummy") ;
		DMPlexConstructGhostCells(dmMesh, "dummy", &ghoscell, &gdm); 
		DMGetLabel( gdm, "vtk", &vtkLabel ) ;
		DMAddLabel(dmMesh, vtkLabel);
		DMDestroy( &gdm );
		//ViewDMLabelsIndex();

		#define print_dmMesh true
		#if ( print_dmMesh == true )
		PetscPrintf( PETSC_COMM_WORLD, "\n") ;
		DMView( dmMesh, PETSC_VIEWER_STDOUT_WORLD) ;
		PetscPrintf( PETSC_COMM_WORLD, "\n") ;
		#endif

		ReadBoundaryCellMarkersFromFile( filename.c_str() ) ;

	/*--- Create cell DM ---*/
		CreateCellDM() ;

	/*--- Get the range of cell indices ( including cpu boundary overlapping cells ) ---*/
		DMPlexGetHeightStratum( dmMesh, 0, &cStart, &cEnd ) ; 
		DMPlexGetDepthStratum ( dmMesh, 0, &vStart, &vEnd ) ;
		DMPlexGetHeightStratum( dmMesh, 1, &fStart, &fEnd ) ; 
		if( GetDimension() == 3 ) DMPlexGetDepthStratum ( dmMesh, 1, &eStart, &eEnd ) ;
		PetscPrintf( PETSC_COMM_SELF, "rank: %d, cStart: %d, cEnd: %d, cEndInterior: %d\n", mpi_rank, cStart, cEnd, cEndInterior ) ;
		PetscPrintf( PETSC_COMM_SELF, "rank: %d, vStart: %d, vEnd: %d,\n", mpi_rank, vStart, vEnd ) ;
		PetscPrintf( PETSC_COMM_SELF, "rank: %d, fStart: %d, fEnd: %d,\n", mpi_rank, fStart, fEnd ) ;
		if( GetDimension() == 3 ) PetscPrintf( PETSC_COMM_SELF, "rank: %d, eStart: %d, eEnd: %d,\n", mpi_rank, eStart, eEnd ) ;

		CreateLocalToGlobalCellIdMapping() ; //Also calculate the ghost cell index.


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
void CDomain::ViewDMLabels()
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
void CDomain::ViewDMLabelsIndex()
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
void CDomain::CreateLocalToGlobalCellIdMapping()
{
		
		/* Note: This gets a borrowed reference, so the user should not destroy this PetscSection. */
		DMGetGlobalSection( dmCell, &global_ids_section ) ;

 		//Vec gVec, lVec ;
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
		//cout<<cEnd<<endl;
		//cout<<GhostCellCount<<endl;
		//PetscEnd() ;
		cEndInterior = cEnd - GhostCellCount ;
		//DMPlexSetHybridBounds(dmCell, cEndInterior, PETSC_DETERMINE, PETSC_DETERMINE, PETSC_DETERMINE);

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
	//PetscEnd();
}
void CDomain::CreateCellDM()
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
void CDomain::ExtractCellGeomInformations()
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
		//typedef struct { PetscReal centroid[3]; PetscReal volume; } PetscFVCellGeom;

		DM dmCell ;
		VecGetDM( CellGeometryVec, &dmCell);
		VecGetArrayRead(CellGeometryVec, &cgeom);
		volume   = new double [ cEnd-cStart ] ;
		centroid = new double [ (cEnd-cStart)*3 ] ;
		//for( PetscInt i=cStart ; i <cEnd ; i++ ) centroid[ i ] = new double [3] ;

		for( PetscInt i=cStart ; i <cEnd ; i++ )
		{
			PetscFVCellGeom *cg ;
			DMPlexPointLocalRead(dmCell, i, cgeom, &cg);
			volume[i] = cg->volume ;
			centroid[i*3+0] = cg->centroid[0] ;
			centroid[i*3+1] = cg->centroid[1] ;
			centroid[i*3+2] = cg->centroid[2] ;
		}
		VecRestoreArrayRead(CellGeometryVec, &cgeom);

	#define print_volume_centroid false
	#if ( print_volume_centroid == true )
		for ( PetscInt i = cStart ; i < cEnd ; i++ ) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, volume[%d]-> nnghbrs: %d\n", mpi_rank, i, volume[i] ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"centroid[0]->%e, centroid[1]->%e, centroid[2]->%e \n",centroid[i*3+0],centroid[i*3+1],centroid[i*3+2]  ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif 

	CreateCellNeighborVector();

}
void CDomain::CreateCellNeighborVector()
{
	string Type="cell_face" ;
	//string Type="cell_cell" ;
	int cell_nnghbrs_cell, face_nnghbrs_cell ;

	int             cell_nghbr_cell_index ;
	const PetscInt *cell_nghbr_face_index ;
	PetscInt cell_nnghbrs_face ;

	if ( Type == "cell_face" ) {


		for( PetscInt i=cStart ; i <cEnd ; i++ ) {

			cell_map_all.insert( make_pair( i, new CCell() ) ); 

			cell_nnghbrs_cell  = cell_cell_nnghbrs(i) ;

			DMPlexGetConeSize( dmMesh, i, &cell_nnghbrs_face) ;

			DMPlexGetCone    ( dmMesh, i, &cell_nghbr_face_index ) ;

			/* cell-cell */
			for ( int k = 0 ; k < cell_nnghbrs_cell ; k++ ) {
				cell_nghbr_cell_index = cell_nghbr(i,k);
				cell_map_all[ i ]->nghbr_cell.push_back(cell_nghbr_cell_index) ;
			}

			/* cell-face */
			for ( int k = 0 ; k < cell_nnghbrs_face ; k++ ) {
				DMPlexGetSupportSize( dmMesh, cell_nghbr_face_index[k], &face_nnghbrs_cell) ;
				//cout<<face_nnghbrs_cell<<endl;
				if ( face_nnghbrs_cell == 1 ) {
					cell_map_all[i]->nghbr_face.push_back(cell_nghbr_face_index[k]) ;
				}
			}
		}//end cell loop
		//cout<<"A"<<endl;PetscEnd();
	}else if(Type == "cell_cell"){

	}
	#define print_cell_information false
	#if ( print_cell_information == true )
		for ( PetscInt i = cStart ; i < cEnd ; i++ ) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, nnghbrs cell: %d, nnghbrs face: %d\n", mpi_rank, cell_map_all[i]->nnghbrs_cell(), cell_map_all[i]->nnghbrs_face()) ;
			//PetscSynchronizedPrintf( PETSC_COMM_WORLD,"centroid[0]->%e, centroid[1]->%e, centroid[2]->%e \n",centroid[i*3+0],centroid[i*3+1],centroid[i*3+2]  ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif 
}
void CDomain::ExtractFaceGeomInformations()
{
	//face = new CFace[ fEnd-fStart ] ;
	PetscInt check_face=0.0 ;

	const PetscInt       *cc ;
	const PetscScalar *fgeom ;

	DM dmFace ;
	VecGetDM( FaceGeometryVec, &dmFace);
	VecGetArrayRead(FaceGeometryVec, &fgeom);

	PetscInt cell_nnghbrs_face ;
	const PetscInt *cell_nghbr_face, *orientations ;

	for ( PetscInt i = fStart ; i < fEnd ; i++ ) {

		face_map_all.insert( make_pair( i, new CFace() ) ); 
		face_map_all[ i ]->nnghbrs = face_nnghbrs_cell(i) ;
		face_map_all[ i ]->offsets = i ;

		DMPlexGetSupport( dmMesh, i, &cc) ;
		if( face_map_all[ i ]->nnghbrs==2 and cc[0] > cc[1] ) cout<<"flip"<<endl;
		PetscFVFaceGeom *fg ;
		DMPlexPointLocalRead(dmFace, i, fgeom, &fg);

		/*Copy face normal & centroid */
		for (PetscInt d = 0 ; d < GetDimension() ; d++ ) {
			face_map_all[ i ]->face_nrml[d] = fg->  normal[d] ;
			face_map_all[ i ]->nA[d] 				= fg->  normal[d] ;
			face_map_all[ i ]-> centroid[d] = fg->centroid[d] ;
		}

		//face[i-fStart].calculate_normal_mag2();
		face_map_all[ i ]->calculate_normal_mag2();
		face_map_all[ i ]->calculate_face_area();

		/* normal vector*/
		for (PetscInt d = 0 ; d < GetDimension() ; d++ ) {
			//face[i-fStart].face_nrml[d] /= face[i-fStart].face_nrml_mag ;
			face_map_all[ i ]->face_nrml[d] /= face_map_all[ i ]->face_nrml_mag ;
		}

		if ( face_map_all[ i ]->nnghbrs == 1 ) {//means that the face is boundary face or processor boundary face.

				face_map_all[ i ]->c0 = cc[0] ;
				face_map_all[ i ]->c1 = cc[0] ;
			// }
			//PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, bcface[%d]->c0-> %d \n", mpi_rank, i, face_map_all[ i ]->c0 ) ;

			face_map_all[ i ]->gc0 = Local2GlobalIds( face_map_all[ i ]->c0 )  ;
			face_map_all[ i ]->gc1 = Local2GlobalIds( face_map_all[ i ]->c1 )  ;

			face_map_all[ i ]->c0_volume = volume[face_map_all[ i ]->c0] ;
			face_map_all[ i ]->c1_volume = volume[face_map_all[ i ]->c1] ;

			face_map_all[ i ]->c0_centroid[0] = centroid[ face_map_all[ i ]->c0*3+0 ] ;
			face_map_all[ i ]->c0_centroid[1] = centroid[ face_map_all[ i ]->c0*3+1 ] ;
			face_map_all[ i ]->c0_centroid[2] = centroid[ face_map_all[ i ]->c0*3+2 ] ;

			face_map_all[ i ]->c1_centroid[0] = centroid[ face_map_all[ i ]->c1*3+0 ] ;
			face_map_all[ i ]->c1_centroid[1] = centroid[ face_map_all[ i ]->c1*3+1 ] ;
			face_map_all[ i ]->c1_centroid[2] = centroid[ face_map_all[ i ]->c1*3+2 ] ;

		} else {
			if( cc[0] > cc[1] ) {
				cout<<"cc[0] > cc[1]"<<endl;
				PetscEnd();
			}
			//Get how many faces arround the cell cc[0].
			DMPlexGetConeSize(dmMesh, cc[0], &cell_nnghbrs_face) ;
			check_face = cell_nnghbrs_face ;

			//Get the face index arrays & orientation array at cc[0] cell.
			DMPlexGetCone    ( dmMesh, cc[0], &cell_nghbr_face ) ;
			DMPlexGetConeOrientation( dmMesh, cc[0], &orientations );

			for (PetscInt j=0 ; j < cell_nnghbrs_face ; j++ ) {
				//cout<<"face : "<<i<<", cell_nghbr_face: "<<cell_nghbr_face[j]<<endl;
				if( cell_nghbr_face[j] == i ) {

					if ( orientations[j] >= 0.0 ) {

						face_map_all[ i ]->c0 = cc[0] ;
						face_map_all[ i ]->c1 = cc[1] ;


						face_map_all[ i ]->gc0 = Local2GlobalIds( face_map_all[ i ]->c0 )  ;
						face_map_all[ i ]->gc1 = Local2GlobalIds( face_map_all[ i ]->c1 )  ;

						face_map_all[ i ]->c0_volume = volume[face_map_all[ i ]->c0] ;
						face_map_all[ i ]->c1_volume = volume[face_map_all[ i ]->c1] ;
						
						face_map_all[ i ]->c0_centroid[0] = centroid[ face_map_all[ i ]->c0*3+0 ] ;
						face_map_all[ i ]->c0_centroid[1] = centroid[ face_map_all[ i ]->c0*3+1 ] ;
						face_map_all[ i ]->c0_centroid[2] = centroid[ face_map_all[ i ]->c0*3+2 ] ;

						face_map_all[ i ]->c1_centroid[0] = centroid[ face_map_all[ i ]->c1*3+0 ] ;
						face_map_all[ i ]->c1_centroid[1] = centroid[ face_map_all[ i ]->c1*3+1 ] ;
						face_map_all[ i ]->c1_centroid[2] = centroid[ face_map_all[ i ]->c1*3+2 ] ;

						//cout<<"Case A"<<endl;
					} else {
						cout<<"AAAAAA"<<endl;	
						face_map_all[ i ]->c0 = cc[0] ;
						face_map_all[ i ]->c1 = cc[1] ;	

						face_map_all[ i ]->gc0 = Local2GlobalIds( face_map_all[ i ]->c0 )  ;
						face_map_all[ i ]->gc1 = Local2GlobalIds( face_map_all[ i ]->c1 )  ;

						face_map_all[ i ]->c0_volume = volume[face_map_all[ i ]->c0] ;
						face_map_all[ i ]->c1_volume = volume[face_map_all[ i ]->c1] ;
						
						face_map_all[ i ]->c0_centroid[0] = centroid[ face_map_all[ i ]->c0*3+0 ] ;
						face_map_all[ i ]->c0_centroid[1] = centroid[ face_map_all[ i ]->c0*3+1 ] ;
						face_map_all[ i ]->c0_centroid[2] = centroid[ face_map_all[ i ]->c0*3+2 ] ;

						face_map_all[ i ]->c1_centroid[0] = centroid[ face_map_all[ i ]->c1*3+0 ] ;
						face_map_all[ i ]->c1_centroid[1] = centroid[ face_map_all[ i ]->c1*3+1 ] ;
						face_map_all[ i ]->c1_centroid[2] = centroid[ face_map_all[ i ]->c1*3+2 ] ;

						for (PetscInt d = 0 ; d < GetDimension() ; d++ ) {
							face_map_all[ i ]->face_nrml[d] = (-1.0)*fg->  normal[d] ;
						}
						//cout<<"Case B"<<endl;
					}
					check_face=check_face-1;
				}//Find face
			}//End cell nnghbrs face
			if ( check_face == cell_nnghbrs_face ) {
				cout<<"check face error"<<endl;
			}
		}
	}//End face loop
	//PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	//VecRestoreArrayRead(FaceGeometryVec, &fgeom);
//PetscEnd();
	#define print_faceNormal  false
	#if (print_faceNormal == true)
		for ( PetscInt i =fStart ; i < fEnd ; i++ ) {

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, face[%d]-> nnghbrs: %d\n", mpi_rank, i, face_map_all[i]->nnghbrs ) ;

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"nghbr[0]-> %d \n", face_map_all[i]->c0 ) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"nghbr[1]-> %d \n", face_map_all[i]->c1 ) ;

			for (PetscInt d = 0 ; d < GetDimension() ; d++ )
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"face_nrml[%d]-> %e \n",d, face_map_all[i]->face_nrml[d] ) ;

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"face_nrml_mag-> %e \n",face_map_all[i]->face_nrml_mag) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"volume[0]-> %e \n",face_map_all[i]->c0_volume) ;
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"volume[1]-> %e \n",face_map_all[i]->c1_volume) ;

			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n") ;
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	#endif
	ComputeFaceCellInformations();
	//PetscEnd();
}
void CDomain::ComputeFaceCellInformations()
{
	for ( PetscInt i = fStart ; i < fEnd ; i++ ) {

		if( face_map_all[ i ]->nnghbrs == 2 ) {

			for ( int k = 0 ; k < 3 ; k++ ) {
				face_map_all[ i ]->PN[ k ] = face_map_all[ i ]->c1_centroid[k] - face_map_all[ i ]->c0_centroid[k] ;

				/*--- Pf Vector ---*/
				face_map_all[ i ]->Pf[ k ] = face_map_all[ i ]->centroid[k] - face_map_all[ i ]->c0_centroid[k] ;
				/*--- Nf Vector ---*/
				face_map_all[ i ]->Nf[ k ] = face_map_all[ i ]->centroid[k] - face_map_all[ i ]->c1_centroid[k] ;

				/*--- P'f Vector ---*/
				face_map_all[ i ]->PPf[ k ] = DotProduct( face_map_all[ i ]->Pf, face_map_all[ i ]->face_nrml )*face_map_all[ i ]->face_nrml[k] ;

				/*--- N'f Vector ---*/
				face_map_all[ i ]->NPf[ k ] = DotProduct( face_map_all[ i ]->Nf, face_map_all[ i ]->face_nrml )*face_map_all[ i ]->face_nrml[k] ;

				/*--- PP'  Vector ( PP' = Pf - P'f ) ---*/
				face_map_all[ i ]->PPP[ k ] = face_map_all[ i ]->Pf[ k ] - face_map_all[ i ]->PPf[ k ] ;

				/*--- NN'  Normal Vector ( NN' = Nf - N'f ) ---*/
				face_map_all[ i ]->NNP[ k ] = face_map_all[ i ]->Nf[ k ] - face_map_all[ i ]->NPf[ k ] ;

			}//end x,y,z

			face_map_all[ i ]->dPN     = sqrt( DotProduct(  face_map_all[ i ]->PN, face_map_all[ i ]->PN ) ) ; 
			face_map_all[ i ]->dPPf    = sqrt( DotProduct( face_map_all[ i ]->PPf, face_map_all[ i ]->PPf ) ) ; 
			face_map_all[ i ]->dNPf    = sqrt( DotProduct( face_map_all[ i ]->NPf, face_map_all[ i ]->NPf ) ) ; 

		/* For boundary */
		} else {
			for ( int k = 0 ; k < 3 ; k++ ) {
				//face_map_all[ i ]->PN[ k ] = face_map_all[ i ]->c1_centroid[k] - face_map_all[ i ]->c0_centroid[k] ;

				/*--- Pf Vector ---*/
				face_map_all[ i ]->Pf[ k ] = face_map_all[ i ]->centroid[k] - face_map_all[ i ]->c0_centroid[k] ;

				/*--- Nf Vector ---*/
				//face_map_all[ i ]->Nf[ k ] = face_map_all[ i ]->centroid[k] - face_map_all[ i ]->c1_centroid[k] ;

				/*--- P'f Vector ---*/
				face_map_all[ i ]->PPf[ k ] = DotProduct( face_map_all[ i ]->Pf, face_map_all[ i ]->face_nrml )*face_map_all[ i ]->face_nrml[k] ;

				/*--- N'f Vector ---*/
				//face_map_all[ i ]->NPf[ k ] = DotProduct( face_map_all[ i ]->Nf, face_map_all[ i ]->face_nrml )*face_map_all[ i ]->face_nrml[k] ;

				/*--- PP'  Vector ( PP' = Pf - P'f ) ---*/
				face_map_all[ i ]->PPP[ k ] = face_map_all[ i ]->Pf[ k ] - face_map_all[ i ]->PPf[ k ] ;

				/*--- NN'  Normal Vector ( NN' = Nf - N'f ) ---*/
				//face_map_all[ i ]->NNP[ k ] = face_map_all[ i ]->Nf[ k ] - face_map_all[ i ]->NPf[ k ] ;
			}//end x,y,z
			face_map_all[ i ]->dPPf    = sqrt( DotProduct( face_map_all[ i ]->PPf, face_map_all[ i ]->PPf ) ) ; 
			face_map_all[ i ]->dPN     = face_map_all[ i ]->dPPf ;
			//face_map_all[ i ]->dNPf    = sqrt( DotProduct( face_map_all[ i ]->NPf, face_map_all[ i ]->NPf ) ) ; 

		}//end face
	}
	PetscBarrier(NULL);
}

void CDomain::construct_cell_face()
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
void CDomain::ReadBoundaryCellMarkersFromFile(string filename)
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
	PetscPrintf(PETSC_COMM_WORLD,"ReadBoundaryCellMarkersFromFile\n") ;
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
	file.close() ;
	file.clear() ;

	#define print_markers_map true
	#if ( print_markers_map == true )
			PetscPrintf( PETSC_COMM_WORLD,"PhysNames\n") ;
		for(auto it = PhysNames.cbegin(); it != PhysNames.cend(); ++it) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, name: %s, ids: %d\n", mpi_rank, it->first.c_str(), it->second ) ;
		}
		PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n" ) ;
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT ) ;
	#endif

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

	#define print_markers_map_CPU true
	#if ( print_markers_map_CPU == true )
			PetscPrintf( PETSC_COMM_WORLD,"PhysNamesCPU\n") ;
		//PetscPrintf( PETSC_COMM_WORLD,"mpi_rank: %d\n", mpi_rank ) ;
		for(auto it = PhysNamesCPU.cbegin(); it != PhysNamesCPU.cend(); ++it) {
			PetscSynchronizedPrintf( PETSC_COMM_WORLD,"mpi_rank: %d, name: %s, ids: %d\n", mpi_rank, it->first.c_str(), it->second ) ;
		}
		PetscSynchronizedPrintf( PETSC_COMM_WORLD,"\n" ) ;
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT ) ;
	#endif

}
 void CDomain::BulidFaceCellLoopMap()
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
			
//PrintfMemory(); 
			for ( int i = 0 ; i < index_stratum_size ; i++ ) {
				face_id = ids[ i ] ;
					if ( face_map_all[face_id]->c0 < cEndInterior ) {
						//power_face_loop.insert( make_pair( face_id, face_map_all[ face_id ] ) ); 
						//cout<<face_map_all[ face_id ]<<endl;
						power_face_loop.push_back( face_map_all[ face_id ] ) ;
						//cout<<power_face_loop2[i]<<endl<<endl;
					}
			}
//PrintfMemory() ;
	
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
				face_id = ids[ i ] ;
				//PetscPrintf(PETSC_COMM_WORLD, "index: %d\n", face_id) ;
				if ( face_map_all[face_id]->c0 < cEndInterior ) {
					//ground_face_loop.insert( make_pair( face_id, face_map_all[ face_id ] ) ); 
					ground_face_loop.push_back( face_map_all[ face_id ] ) ;
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
				face_id = ids[ i ] ;
				if ( face_map_all[face_id]->c0 < cEndInterior ) {
					neumann_face_loop.push_back( face_map_all[ face_id ] ) ;
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
			m = interior_face_ids[i] ;
			if ( face_map_all[m]->nnghbrs == 2 )
			{
				if ( face_map_all[m]->c1 >= cEndInterior ) {
					processors_face_loop.push_back( face_map_all[ m ] ) ;
				}else	if ( face_map_all[m]->c0 >= cEndInterior ) {
					//cout<<"face: "<<m<<"\t"<<"c0: "<<face_map_all[m]->c0<<"\t"<<"c1: "<<face_map_all[m]->c1<<", int: "<<cEndInterior<<endl;
					//cout<<"error"<<endl;
					processors_face_loop.push_back( face_map_all[ m ] ) ;
				} else {
					interior_face_loop.push_back( face_map_all[ m ] ) ;
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
			//interior_face_loop
		//BulidFaceCellLoopMap
		PetscPrintf( PETSC_COMM_WORLD,"BulidFaceCellLoopMap done ...\n",mpi_rank) ;

}
