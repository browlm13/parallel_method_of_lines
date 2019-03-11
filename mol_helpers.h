#ifndef MOL_HELPERS
#define MOL_HELPERS

// macros
#define TRUE 1
#define FALSE 0

////
//   Setting Mesh / Lattice Dimensions

int 	set_C_dims( 	int n, & int [N_DIMS] C_dims, int exact );

void 	set_local_dims( int [N_DIMS] global_dims, int [N_DIMS] C_dims,
						int [N_DIMS] C_coords, & int [N_DIMS] local_dims );

void 	set_local_memory_dims( int [N_DIMS][N_FACES] neighbor_ranks, int [MAX_FEILD][N_FACES] boundry_padding,
						int [MAX_FEILD][MAX_AXIS] feild_communications,
						int [N_DIMS] local_dims, & int [N_FEILDS][N_DIMS] local_memory_dims  );

////
//   Setting Process Communication Lattice Variables :  MPI Communicator, process id, process coordinates, neighboring process ides

void 	set_comm_lattice_vars( & MPI_Comm C_comm, & int rank, & int [N_DIMS] C_coords, & int [N_DIMS][N_FACES] neighbor_ranks );


#endif
