#ifndef MOL_HELPERS
#define MOL_HELPERS

// macros
#define TRUE 1
#define FALSE 0

////
//   Setting Mesh / Lattice Dimensions

int 	set_comm_dims( 	int n, & int [3] comm_dims, int exact );

void 	set_local_dims( int [3] global_dims, int [3] comm_dims,
						int [3] comm_coords, & int [3] local_dims );

void 	set_L_mem_dims( int [3][2] neighbor_pids, int [MAX_FEILD][2] boundry_padding,
						int [MAX_FEILD][MAX_AXIS] feild_communications,
						int [3] local_dims, & int [4][3] L_mem_dims  );

////
//   Setting Process Communication Lattice Variables :  MPI Communicator, process id, process coordinates, neighboring process ides

void 	set_comm_lattice_vars( & MPI_Comm C_comm, & int pid, & int [3] comm_coords, & int [3][2] neighbor_pids );


#endif
