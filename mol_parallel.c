// internal
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// external
#include <mpi.h>

// mylib
#include "mol_helpers.h"

/*


          *  *  *  *
          * *  *  ** *
          * * *  * **  *
          * * *  * **  *
          * *** ** **  *
            * ** ** ** *
              *  *  *  *

  C ~ process communication lattice
  L ~ Local lattice mesh

*/


////
// Mesh Settings

// Mesh Size -- m
const int m = 30;    // [TODO]: feild dimensions as input parameter

////
// Settings

const int N_DIMS = 3;                                           // number of dimensions in topology
const int N_FEILDS = 4;
const int N_FACES = 2;

int main( int argc, char** argv ) {

  ////
  // Global Mesh Dimensions (Does Not Correspond To Memory Allocated)

  const int [N_DIMS] global_dims = {m-1, m-1, m-1};               // global mesh interior dimesions ( no walls, no ghosts ) -- "Not the actuall memory allocated"

  ////
  // Boundry Conditions (Periodicity Of Feilds On Given Axes)                       [TODO]: Make flexible
  // int [N_FEILDS][N_DIMS] feild_periodicty = {FALSE};           // Not periodic

  ////
  // Extra Memory Allocation For Boundry Conditions

  int [N_FEILDS][N_FACES] boundry_padding = {0};                  // velocity_feild_dims = {m-1, m-1, m-1};
  boundry_padding[P_INDEX][FACE_1] = 1;                           // pressure_feild_dims = {m+1, m+1, m+1};
  boundry_padding[P_INDEX][FACE_2] = 1;

  ////
  // Feild Axes of Communications With Neighboring Faces

  int [N_FEILDS][N_DIMS] feild_communications = {FALSE};          // feild axes of communication:
  feild_communications[P_INDEX][X] = TRUE;                        //         P on (x,
  feild_communications[P_INDEX][Y] = TRUE;                        //               y,
  feild_communications[P_INDEX][Z] = TRUE;                        //               z)
  feild_communications[U_INDEX][X] = TRUE;                        //         U on (x),
  feild_communications[V_INDEX][Y] = TRUE;                        //         V in (y),
  feild_communications[W_INDEX][Z] = TRUE;                        //         W on (z)

  ////
  // Process Communication Lattice

  MPI_Comm  C_comm;                                               // process comm latice communicator
  int       C_dims     [N_DIMS];                               // number of ranks per dimension in process communication lattice
  int       C_periodicity [N_DIMS] = {FALSE, FALSE, FALSE};       // periodicity of comm lattice axes
  int       reorder                = FALSE;                       // allow reordering

  ////
  // Init MPI

  int n_procs;                                                    // Number of processes
  MPI_Init       ( & argc, & argv );                              // MPI init call
  MPI_Comm_size  ( MPI_COMM_WORLD, & n_procs );                   // Retreive number of processes

  ////
  // Init C_comm                                                  // C_comm ~ Process Communication Lattice MPI Communicator Handel

  int exact = TRUE;                                               // Don't allow unused processes : exact = TRUE
  set_C_dims  ( n_procs, C_dims, exact );                   // set number of processes per axis in process comm lattice
  MPI_Dims_create( n_procs, N_DIMS, C_dims );                  // not positive what this is doing...
  MPI_Cart_create( MPI_COMM_WORLD, N_DIMS, C_dims,
                   C_periodicity, reorder, & C_comm );            // Set up Processes Communication Lattice

  ////
  // Process Specific Data

  int                    rank;                                     // Process Rank / ID
  int [N_DIMS][N_FACES]  neighbor_ranks;                           // Neighboring Processes Ranks / ID's (2 per axis)
  int [N_DIMS]           C_coords;                             // Coordinates of Process in Process Communication Lattice
  int [N_DIMS]           local_dims;                              // Number of Mesh Points Assigned to Process
  int [N_FEILDS][N_DIMS] local_memory_dims;                       // Memory Allocation for Process

  // Init:    rank,  C_coords,  neighbor_ranks
  set_comm_lattice_vars( & MPI_Comm C_comm, & rank, & C_coords, & neighbor_ranks );

  // Init:    local dims
  set_local_dims       ( global_dims, C_dims, C_coords, & local_dims );

  // Init:    local_memory_dims
  set_local_memory_dims       ( neighbor_ranks, boundry_padding, feild_communications,
                                local_dims, & local_memory_dims );

  ////
  // Allocate memory for all 16 local feilds :

  //   P, U, V, W,   P_next, U_next, V_next, W_next,   P_g1, U_g1, V_g1, W_g1,   P_g2, U_g2, V_g2, W_g2
  double [local_memory_dims[P_INDEX][X]][local_memory_dims[P_INDEX][Y]][local_memory_dims[P_INDEX][Z]] P, P_next, P_g1, P_g2;
  double [local_memory_dims[U_INDEX][X]][local_memory_dims[U_INDEX][Y]][local_memory_dims[U_INDEX][Z]] U, U_next, U_g1, U_g2;
  double [local_memory_dims[V_INDEX][X]][local_memory_dims[V_INDEX][Y]][local_memory_dims[V_INDEX][Z]] V, V_next, V_g1, V_g2;
  double [local_memory_dims[W_INDEX][X]][local_memory_dims[W_INDEX][Y]][local_memory_dims[W_INDEX][Z]] W, W_next, W_g1, W_g2;

  // [TODO]: Fix comm wall mappings
  // [TODO]: To make more flexable and readable initalize all Feilds to true solutions in a single function call
  // [TODO]: Use global mesh indices for true solution functions
  // [TODO]: Hide two stages of sending and recieving corners
  // [TODO]: Define problem specs/settings in seperate file
  // [TODO]: Assign 1 processes to I/O
  // [TODO]: Put enums and consts in correct places
  // [TODO]: Add Timers

  ////
  // Set time to MIN_T

  ////
  // Set Feilds To True Solution For First Step

  ////
  // Begin RK4

  //      - sync at each k
  //      - parallel graphs

  ////
  // Output Error

  ////
  // Shutdown MPI
  MPI_Finalize();
}
