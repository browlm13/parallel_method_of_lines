#include <stdlib.h>
#include <math.h>

////
// macros

#define TRUE 1
#define FALSE 0

////
// Enums

enum axis_indices {
  X = 0,
  Y = 1,
  Z = 2,
  MIN_AXIS = 0,
  MAX_AXIS = 2
};

enum feild_indices {
  P_INDEX = 0,
  U_INDEX = 1,
  V_INDEX = 2,
  W_INDEX = 3,
  MIN_FEILD = 0,
  MAX_FEILD = 3
};

enum faces {
  FACE_1 = 0,
  FACE_2 = 1,
  MIN_FACE = 0,
  MAX_FACE = 1
};

int set_C_dims( int n, & int [N_DIMS] C_dims, int exact ){
  /*

      Chooses dimensions of process comunication lattice for input parameter to MPI_Cart_create,

                                int [3] dims = [ nx , ny, nz ]


      Select nx, ny, nz by minimizin the cost function, J(n, nx, ny, nz), given a fixed n,

                                    J(n,nx,ny,nz) =

          (nx-ny)**( 2 * weight_1 ) + (ny-nz)**( 2 * weight_1 ) + (product-n)**( 2 * weight_2 )

                                         min   ( J )
                                       nx,ny,nz

      where the two weights in the exponets, weight_1 and weight_2, correspond to the two prefrences:

                                1.) nx ~= ny ~= nz,
                                2.) nx*ny*nz ~= n


      int n                   ~ total number of processes available
      int [3] C_dims            ~ [nx, ny, nz]
      int exact               ~ boolean whether or not to allow unused processes
      returns int remainder   ~ return number of unused processes
  */

  // declarations
  int old_cost; int new_cost;
  int product;
  int remainder;

  // cost function weights
  int weight_1 = 1;                                          // preference  : nx ~= ny ~= nz
  int weight_2 = 1;                                          // preference  : nx*ny*nz ~= n

  // initilize ranks per dimension and their cost
  old_cost = 2*(n-1)**2;
  C_dims[0] = n; C_dims[1] = 1; C_dims[2] = 1;      // init values : [nx,ny,nz] = [n,1,1]

  int s = ceil((double)n/2.0);
  for nz in range(1, s+1){
  for ny in range(1, s+1){
  for nx in range(1, s+1){

      // compute product of nx, ny and nz
      product = nx*ny*nz;
      if (product <= n){

        // compute cost function J(n,nx,ny,nz)
        new_cost = (nx-ny)**(2*weight_1) + (ny-nz)**(2*weight_1) + (product-n)**(2*weight_2);

        // if use of all processes is required (obliterated by primes)
        if (exact){
          if (product - n !=0){ new_cost = old_cost; }
        }

        // update minimum
        if(new_cost < old_cost){
          C_dims[0] = i; C_dims[1] = j; C_dims[2] = k;
          old_cost = new_cost;
        }

      }

  } // nx
  } // ny
  } // nz

  // compute number of processes not used
  remainder = n - C_dims[0] * C_dims[1] * C_dims[2];

  // return number of left over processes
  return remainder;
}


void set_comm_lattice_vars( & MPI_Comm C_comm, & int rank, & int [N_DIMS] C_coords, & int [N_DIMS][N_FACES] neighbor_ranks ){
  // Note: local_dims are process specific
  // Init:    rank,  C_coords,  neighbor_ranks

  // set rank ~ process id / rank
  MPI_Comm_rank( C_comm, & rank );

  // set C_coords ~ process coordinates in communication lattice
  MPI_Cart_coords( C_comm, rank, N_DIMS, & C_coords );

  // set neighbor_ranks [axis] [face] ~ neighbor process ids / ranks
  int disp = 1;
  MPI_Cart_shift( C_comm, X, disp, & neighbor_ranks[X][FACE_1], & neighbor_ranks[X][FACE_2] );
  MPI_Cart_shift( C_comm, Y, disp, & neighbor_ranks[Y][FACE_1], & neighbor_ranks[Y][FACE_2] );
  MPI_Cart_shift( C_comm, Z, disp, & neighbor_ranks[Z][FACE_1], & neighbor_ranks[Z][FACE_2] );

}

void set_local_dims( int [N_DIMS] global_dims, int [N_DIMS] C_dims, int [N_DIMS] C_coords, & int [N_DIMS] local_dims ){
  // Note: local_dims are process specific

  // declarations outside if statments for readability
  int [N_DIMS] remainders;

  for ( int axis  = MIN_AXIS;  axis  <= MAX_AXIS;  axis++  ){

      // compute local lattice dimensions ( ignoring remainders )
      local_dims[axis] = floor( (double)global_dims[axis] / (double)C_dims[axis] );

      // remainders outside if statments for readability
      remainders[axis] = global_dims[axis] - local_dims[axis] * C_dims[axis];

      // add remainders to exception faces
      if ( C_coords[axis] == C_dims[axis] - 1 ): local_dims[axis] += remainders[axis];
  }
}

void set_local_memory_dims(  int [N_DIMS][N_FACES] neighbor_ranks, int [N_DIMS][N_FACES] boundry_padding, int [N_FEILDS][N_DIMS] feild_communications, int [N_DIMS] local_dims, & int [N_FEILDS][N_DIMS] local_memory_dims  ){
  /*
      Sets up memory indexing for process.

        - Can add storage for boundry values for edge processes.
        - Can add storage for neighbor communications.
        - Can add storage for multiple feilds.

      Note: local_memory_dims are process specific.
  */

  // communication padding
  int difference_method_order = 4;

  for ( int feild = MIN_FEILD; feild <= MAX_FEILD; feild++ ){
  for ( int axis  = MIN_AXIS;  axis  <= MAX_AXIS;  axis++  ){

    // copy local dims
    local_memory_dims[feild][axis] = local_dims[axis];

    // add communication padding if feild communicates along the axis AND process is not an edge
    if ( feild_communications[feild][axis] == TRUE ){
      if ( neighbor_ranks[axis][FACE_1] != MPI_PROC_NULL && neighbor_ranks[axis][FACE_2] != MPI_PROC_NULL ): local_memory_dims[feild][axis] += difference_method_order;
    }

    // add boundry padding for edge processes
    if ( neighbor_ranks[axis][FACE_1] == MPI_PROC_NULL ): local_memory_dims[feild][axis] += boundry_padding[feild][FACE_1] + difference_method_order/2;
    if ( neighbor_ranks[axis][FACE_2] == MPI_PROC_NULL ): local_memory_dims[feild][axis] += boundry_padding[feild][FACE_2] + difference_method_order/2;
  }
  }
}

/*

      Maps / Coordinate Conversions:

        - memory_to_local
        - local_to_memory
        - local_to_global
        - global_to_local

*/

int [3] memory_to_local( int feild_index, & int [N_DIMS] process_coors, & int [N_DIMS] memory_coors, & int [N_FEILDS][N_DIMS] local_memory_dims, & int [N_DIMS] local_dims ){
  /*
      Converts local memory coordinates for specified feild into local coordinates.

      Fm[xm][ym][zm] -> Fl[xl][yl][zl]

      local_memory_dims[X] - local_dims[X] = x_pad;
      local_memory_dims[Y] - local_dims[Y] = y_pad;

      xm
  */

  int [3] local_coors = {
    // [TODO]: remove local_memory_dims to local_memory_dims discrepency and you have local coordinates

    //local_memory_dims[feild_index]


  };

  return local_coors;
}

int [3] local_to_global( & int [N_DIMS] process_coors, & int [N_DIMS] local_coors, & int [N_DIMS] local_dims ){

  int [3] global_coors = {
      (double)local_dims[X] * (double)process_coors[X] + (double)local_coors[X],    // xg = lnx*xc + xl
      (double)local_dims[Y] * (double)process_coors[Y] + (double)local_coors[Y],    // yg = lny*yc + yl
      (double)local_dims[Z] * (double)process_coors[Z] + (double)local_coors[Z]     // zg = lnz*zc + zl
  };

  return global_coors;
}


//int                     is_ghost  ( int [N_DIMS] local_coors, int [3] process_coors, int [3] local_dims );
//int [N_FEILDS][N_FACES] comm_faces( int [N_DIMS] local_coors, int [3] process_coors, int [3] local_dims),
//                                    int [N_FEILDS][N_DIMS] feild_communications );



int [3] global_to_world( & int [N_DIMS] global_coors, double h ){

  int [3] world_coors = {
      h * (double)global_coors[0],    // xw = h*gx
      h * (double)global_coors[1],    // yw = h*gy
      h * (double)global_coors[2]     // zw = h*gz
  };

  return world_coors;
}
