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

int set_comm_dims( int n, & int [3] comm_dims, int exact ){
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
      int [3] comm_dims            ~ [nx, ny, nz]
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
  comm_dims[0] = n; comm_dims[1] = 1; comm_dims[2] = 1;      // init values : [nx,ny,nz] = [n,1,1]

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
          comm_dims[0] = i; comm_dims[1] = j; comm_dims[2] = k;
          old_cost = new_cost;
        }

      }

  } // nx
  } // ny
  } // nz

  // compute number of processes not used
  remainder = n - comm_dims[0] * comm_dims[1] * comm_dims[2];

  // return number of left over processes
  return remainder;
}

void set_comm_lattice_vars( & MPI_Comm C_comm, & int pid, & int [3] comm_coords, & int [3][2] neighbor_pids ){
  // Note: local_dims are process specific
  // Init:    pid,  comm_coords,  neighbor_pids

  // set pid ~ process id / rank
  MPI_Comm_rank( C_comm, & pid );

  // set comm_coords ~ process coordinates in communication lattice
  MPI_Cart_coords( C_comm, pid, n_dims, & comm_coords );

  // set neighbor_pids [axis] [face] ~ neighbor process ids / ranks
  int disp = 1;
  MPI_Cart_shift( C_comm, X, disp, & neighbor_pids[X][FACE_1], & neighbor_pids[X][FACE_2] );
  MPI_Cart_shift( C_comm, Y, disp, & neighbor_pids[Y][FACE_1], & neighbor_pids[Y][FACE_2] );
  MPI_Cart_shift( C_comm, Z, disp, & neighbor_pids[Z][FACE_1], & neighbor_pids[Z][FACE_2] );

}

void set_local_dims( int [n_dims] global_dims, int [n_dims] comm_dims, int [n_dims] comm_coords, & int [n_dims] local_dims ){
  // Note: local_dims are process specific

  // declarations outside if statments for readability
  int [n_dims] remainders;

  for ( int axis  = MIN_AXIS;  axis  <= MAX_AXIS;  axis++  ){

      // compute local lattice dimensions ( ignoring remainders )
      local_dims[axis] = floor( (double)global_dims[axis] / (double)comm_dims[axis] );

      // remainders outside if statments for readability
      remainders[axis] = global_dims[axis] - local_dims[axis] * comm_dims[axis];

      // add remainders to exception faces
      if ( comm_coords[axis] == comm_dims[axis] - 1 ): local_dims[axis] += remainders[axis];
  }
}

void set_L_mem_dims(  int [n_dims][n_faces] neighbor_pids, int [n_dims][n_faces] boundry_padding, int [n_feilds][n_dims] feild_communications, int [n_dims] local_dims, & int [n_feilds][n_dims] L_mem_dims  ){
  /*
      Sets up memory indexing for process.

        - Can add storage for boundry values for edge processes.
        - Can add storage for neighbor communications.
        - Can add storage for multiple feilds.

      Note: L_mem_dims are process specific.
  */

  // communication padding
  int difference_method_order = 4;

  for ( int feild = MIN_FEILD; feild <= MAX_FEILD; feild++ ){
  for ( int axis  = MIN_AXIS;  axis  <= MAX_AXIS;  axis++  ){

    // copy local dims
    L_mem_dims[feild][axis] = local_dims[axis];

    // add communication padding if feild communicates along the axis AND process is not an edge
    if ( feild_communications[feild][axis] == TRUE ){
      if ( neighbor_pids[axis][FACE_1] != MPI_PROC_NULL && neighbor_pids[axis][FACE_2] != MPI_PROC_NULL ): L_mem_dims[feild][axis] += difference_method_order;
    }

    // add boundry padding for edge processes
    if ( neighbor_pids[axis][FACE_1] == MPI_PROC_NULL ): L_mem_dims[feild][axis] += boundry_padding[feild][FACE_1] + difference_method_order/2;
    if ( neighbor_pids[axis][FACE_2] == MPI_PROC_NULL ): L_mem_dims[feild][axis] += boundry_padding[feild][FACE_2] + difference_method_order/2;
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

int [3] memory_to_local( int feild_index, & int [n_dims] process_coors, & int [n_dims] memory_coors, & int [n_feilds][n_dims] L_mem_dims, & int [n_dims] local_dims ){
  /*
      Converts local memory coordinates for specified feild into local coordinates.

      Fm[xm][ym][zm] -> Fl[xl][yl][zl]

      L_mem_dims[X] - local_dims[X] = x_pad;
      L_mem_dims[Y] - local_dims[Y] = y_pad;

      xm
  */

  int [3] local_coors = {
    // [TODO]: remove L_mem_dims to L_mem_dims discrepency and you have local coordinates

    //L_mem_dims[feild_index]


  };

  return local_coors;
}

int [3] local_to_global( & int [3] process_coors, & int [3] local_coors, & int [3] local_dims ){

  int [3] global_coors = {
      (double)local_dims[0] * (double)process_coors[0] + (double)local_coors[0],    // xg = lnx*xc + xl
      (double)local_dims[1] * (double)process_coors[1] + (double)local_coors[1],    // yg = lny*yc + yl
      (double)local_dims[2] * (double)process_coors[2] + (double)local_coors[2]     // zg = lnz*zc + zl
  };

  return global_coors;
}


//int                     is_ghost  ( int [n_dims] local_coors, int [3] process_coors, int [3] local_dims );
//int [n_feilds][n_faces] comm_faces( int [n_dims] local_coors, int [3] process_coors, int [3] local_dims),
//                                    int [n_feilds][n_dims] feild_communications );



int [3] global_to_world( & int [3] global_coors, double h ){

  int [3] world_coors = {
      h * (double)global_coors[0],    // xw = h*gx
      h * (double)global_coors[1],    // yw = h*gy
      h * (double)global_coors[2]     // zw = h*gz
  };

  return world_coors;
}
