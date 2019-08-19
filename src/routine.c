/*
This includes overlapping benchmarks and rouitnes separately implemented.
The number of MPI processes can be controlled.
The routine and tuning knobs can be configured using config.mk.
It outputs average communication, computation and overlap time agaist the increasing delays
*/

#include "includes/routine.h"

#include <stdio.h>

#include <mpi.h>

#include <stdlib.h>

#include <string.h>

#include "includes/overlap_functions.h"

void Ibcast(int argc, char * argv[]) {

  MPI_Init(NULL, NULL);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  int noOfProcs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, & noOfProcs);
  MPI_Status status;
  MPI_Request request;

  double delay = 0.0, delay_min = 0.0, delay_max = 0.0;
  int buffer = 0;
  int i = 0, j = 0, ITER = 0, N_readings = 1;
  SetValuesOverlap(argc, argv, & ITER, & buffer, & delay_min, & delay_max, & N_readings);

  if (rank == 0) output(argc, argv, ITER, buffer, noOfProcs);

  double t_pure = 0.0, t_comp = 0.0, t_overlap = 0.0; //Pure communication / Pure computation / Overlap time
  int root = 0; //for collective operations root is the sender

  double * pure_time, * comp_time, * overall_time; //Stores values of each run and then used for averaging
  pure_time = (double * ) malloc(sizeof(double) * N_readings);
  comp_time = (double * ) malloc(sizeof(double) * N_readings);
  overall_time = (double * ) malloc(sizeof(double) * N_readings);

  double * vec;
  vec = (double * ) malloc(buffer * sizeof(double));

  if (rank != -1) {
    for (j = 0; j < buffer; j++) {
      vec[j] = (double) j;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /////////////////////////////////////////////////////////////////////////////
  /***********************IBROADCAST******************************************/
  for (delay = delay_min; delay < delay_max; delay += 0.005) {
    /////////////////// PURE COMMUNICATION /////////////////////////////////////
    for (j = 0; j < N_readings; j++) {
      for (i = 0; i < ITER; i++) {
        root = i % noOfProcs;
        t_pure -= MPI_Wtime();
        MPI_Ibcast(vec, buffer, MPI_DOUBLE, root, MPI_COMM_WORLD, & request);
        MPI_Wait( & request, & status);
        t_pure += MPI_Wtime();
      }
      t_pure /= ITER;

      ///////////////////////////////////////////////////////////////////////////

      //////////////// OVERLAP REGION   ////////////////////////////////////////////
      MPI_Barrier(MPI_COMM_WORLD);
      for (i = 0; i < ITER; i++) {
        root = i % noOfProcs;
        t_overlap -= MPI_Wtime();
        MPI_Ibcast(vec, buffer, MPI_DOUBLE, root, MPI_COMM_WORLD, & request);
        t_comp -= MPI_Wtime();
        Computation(vec, delay, buffer);
        t_comp += MPI_Wtime();
        MPI_Wait( & request, & status);
        t_overlap += MPI_Wtime();
      }
      t_comp /= ITER;
      t_overlap /= ITER;

      if (rank == 0) {
        pure_time[j] = t_pure;
        comp_time[j] = t_comp;
        overall_time[j] = t_overlap;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) reading_output(pure_time, comp_time, overall_time, N_readings);
  }
  ////////////////////////////////////////////////////////////////////////////

  free(vec);
  free(pure_time);
  free(comp_time);
  free(overall_time);
  MPI_Finalize();
}

void Igather(int argc, char * argv[]) {
  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  int noOfProcs;
  MPI_Comm_size(MPI_COMM_WORLD, & noOfProcs);
  MPI_Status status;
  MPI_Request request;

  double delay = 0.0, delay_min = 0.0, delay_max = 0.0;
  int buffer = 0;
  int i = 0, j, ITER = 0, N_readings;

  SetValuesOverlap(argc, argv, & ITER, & buffer, & delay_min, & delay_max, & N_readings);

  if (rank == 0) output(argc, argv, ITER, buffer, noOfProcs);

  double t_pure = 0.0, t_comp = 0.0, t_overlap = 0.0;
  int root = 0; //for collective operations root is the sender

  int send_size, recv_size; //for scatter and gather
  send_size = buffer / noOfProcs;
  recv_size = buffer / noOfProcs;

  double * pure_time, * comp_time, * overall_time;
  pure_time = (double * ) malloc(sizeof(double) * N_readings);
  comp_time = (double * ) malloc(sizeof(double) * N_readings);
  overall_time = (double * ) malloc(sizeof(double) * N_readings);

  double * send_vec;
  send_vec = (double * ) malloc(send_size * sizeof(double));

  double * vec; // for scatter this is the sender, for gather this is the receiver
  vec = (double * ) malloc(buffer * sizeof(double));

  if (rank != -1) {
    for (j = 0; j < send_size; j++) {
      send_vec[j] = (double) j;
    }
  }
  if (rank != -1) {
    for (j = 0; j < buffer; j++) {
      vec[j] = (double) j;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /////////////////////////////////////////////////////////////////////////////
  /***********************IGATHER******************************************/
  for (delay = delay_min; delay < delay_max; delay += 0.005) {
    /////////////////// PURE COMMUNICATION /////////////////////////////////////
    for (j = 0; j < N_readings; j++) {
      for (i = 0; i < ITER; i++) {
        root = i % noOfProcs;
        t_pure -= MPI_Wtime();
        MPI_Igather(send_vec, send_size, MPI_DOUBLE, vec, recv_size, MPI_DOUBLE, root, MPI_COMM_WORLD, & request);
        MPI_Wait( & request, & status);
        t_pure += MPI_Wtime();
      }
      t_pure /= ITER;

      ///////////////////////////////////////////////////////////////////////////
      //////////////// OVERLAP REGION   ////////////////////////////////////////////
      MPI_Barrier(MPI_COMM_WORLD);
      for (i = 0; i < ITER; i++) {
        root = i % noOfProcs;
        t_overlap -= MPI_Wtime();
        MPI_Igather(send_vec, send_size, MPI_DOUBLE, vec, recv_size, MPI_DOUBLE, root, MPI_COMM_WORLD, & request);

        t_comp -= MPI_Wtime();
        Computation(vec, delay, buffer);
        t_comp += MPI_Wtime();
        MPI_Wait( & request, & status);
        t_overlap += MPI_Wtime();
      }

      t_comp /= ITER;
      t_overlap /= ITER;
      if (rank == 0) {
        pure_time[j] = t_pure;
        comp_time[j] = t_comp;
        overall_time[j] = t_overlap;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) reading_output(pure_time, comp_time, overall_time, N_readings);
  }
  ////////////////////////////////////////////////////////////////////////////
  free(vec);
  free(pure_time);
  free(comp_time);
  free(overall_time);
  MPI_Finalize();
}

void Ireduce(int argc, char * argv[]) {
  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  int noOfProcs;
  MPI_Comm_size(MPI_COMM_WORLD, & noOfProcs);

  MPI_Status status;
  MPI_Request request;

  double delay = 0.0, delay_min = 0.0, delay_max = 0.0;
  int buffer = 0;
  int i = 0, j, ITER = 0, N_readings;

  SetValuesOverlap(argc, argv, & ITER, & buffer, & delay_min, & delay_max, & N_readings);

  if (rank == 0) output(argc, argv, ITER, buffer, noOfProcs);

  double t_pure = 0.0, t_comp = 0.0, t_overlap = 0.0;
  int root = 0; //for collective operations root is the sender

  double * pure_time, * comp_time, * overall_time;
  pure_time = (double * ) malloc(sizeof(double) * N_readings);
  comp_time = (double * ) malloc(sizeof(double) * N_readings);
  overall_time = (double * ) malloc(sizeof(double) * N_readings);

  double * send_vec, * recv_vec;
  send_vec = (double * ) malloc(buffer * sizeof(double));
  recv_vec = (double * ) malloc(buffer * sizeof(double));

  double * vec; // for scatter this is the sender, for gather this is the receiver
  vec = (double * ) malloc(buffer * sizeof(double));

  if (rank == root) {
    for (j = 0; j < buffer; j++) {
      vec[j] = (double) 1.0;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /////////////////////////////////////////////////////////////////////////////
  /***********************IREDUCE******************************************/
  for (delay = delay_min; delay < delay_max; delay += 0.005) {
    /////////////////// PURE COMMUNICATION /////////////////////////////////////
    for (j = 0; j < N_readings; j++) {
      for (i = 0; i < ITER; i++) {
        root = i % noOfProcs;
        t_pure -= MPI_Wtime();
        MPI_Ireduce(send_vec, recv_vec, buffer, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, & request);
        MPI_Wait( & request, & status);
        t_pure += MPI_Wtime();
      }
      t_pure /= ITER;

      ///////////////////////////////////////////////////////////////////////////
      //////////////// OVERLAP REGION   ////////////////////////////////////////////
      MPI_Barrier(MPI_COMM_WORLD);
      for (i = 0; i < ITER; i++) {
        root = i % noOfProcs;
        t_overlap -= MPI_Wtime();
        MPI_Ireduce(send_vec, recv_vec, buffer, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, & request);
        t_comp -= MPI_Wtime();
        Computation(vec, delay, buffer);
        t_comp += MPI_Wtime();
        MPI_Wait( & request, & status);
        t_overlap += MPI_Wtime();
      }
      t_comp /= ITER;
      t_overlap /= ITER;
      if (rank == 0) {
        pure_time[j] = t_pure;
        comp_time[j] = t_comp;
        overall_time[j] = t_overlap;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) reading_output(pure_time, comp_time, overall_time, N_readings);
  }
  ////////////////////////////////////////////////////////////////////////////
  free(vec);
  free(pure_time);
  free(comp_time);
  free(overall_time);
  MPI_Finalize();
}

void Iscatter(int argc, char * argv[]) {
  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  int noOfProcs;
  MPI_Comm_size(MPI_COMM_WORLD, & noOfProcs);

  MPI_Status status;
  MPI_Request request;

  double delay = 0.0, delay_min = 0.0, delay_max = 0.0;
  int buffer = 0;
  int i = 0, j, ITER = 0, N_readings;

  SetValuesOverlap(argc, argv, & ITER, & buffer, & delay_min, & delay_max, & N_readings);

  if (rank == 0) output(argc, argv, ITER, buffer, noOfProcs);

  double t_pure = 0.0, t_comp = 0.0, t_overlap = 0.0;
  int root = 0; //for collective operations root is the sender

  int send_size, recv_size; //for scatter and gather
  send_size = buffer / noOfProcs;
  recv_size = buffer / noOfProcs;

  double * pure_time, * comp_time, * overall_time;
  pure_time = (double * ) malloc(sizeof(double) * N_readings);
  comp_time = (double * ) malloc(sizeof(double) * N_readings);
  overall_time = (double * ) malloc(sizeof(double) * N_readings);

  double * recv_vec;
  recv_vec = (double * ) malloc(recv_size * sizeof(double));

  double * vec; // for scatter this is the sender, for gather this is the receiver
  vec = (double * ) malloc(buffer * sizeof(double));

  if (rank == root) {
    for (j = 0; j < buffer; j++) {
      vec[j] = (double) j;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /////////////////////////////////////////////////////////////////////////////
  /***********************ISCATTER******************************************/
  for (delay = delay_min; delay < delay_max; delay += 0.005) {

    /////////////////// PURE COMMUNICATION /////////////////////////////////////
    for (j = 0; j < N_readings; j++) {
      for (i = 0; i < ITER; i++) {
        root = i % noOfProcs;
        t_pure -= MPI_Wtime();
        MPI_Iscatter(vec, send_size, MPI_DOUBLE, recv_vec, recv_size, MPI_DOUBLE, root, MPI_COMM_WORLD, & request);
        MPI_Wait( & request, & status);
        t_pure += MPI_Wtime();
      } 
      t_pure /= ITER;
      ///////////////////////////////////////////////////////////////////////////
      //////////////// OVERLAP REGION   ////////////////////////////////////////////
      MPI_Barrier(MPI_COMM_WORLD);
      for (i = 0; i < ITER; i++) {
        root = i % noOfProcs;
        t_overlap -= MPI_Wtime();
        MPI_Iscatter(vec, send_size, MPI_DOUBLE, recv_vec, recv_size, MPI_DOUBLE, root, MPI_COMM_WORLD, & request);
        t_comp -= MPI_Wtime();
        Computation(vec, delay, buffer);
        t_comp += MPI_Wtime();
        MPI_Wait( & request, & status);
        t_overlap += MPI_Wtime();
      }
      t_comp /= ITER;
      t_overlap /= ITER;
      if (rank == 0) {
        pure_time[j] = t_pure;
        comp_time[j] = t_comp;
        overall_time[j] = t_overlap;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) reading_output(pure_time, comp_time, overall_time, N_readings);
  }
  free(vec);
  free(pure_time);
  free(comp_time);
  free(overall_time);
  MPI_Finalize();
}

void Isend(int argc, char * argv[]) {
  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  int noOfProcs;
  MPI_Comm_size(MPI_COMM_WORLD, & noOfProcs);
  MPI_Status status;
  MPI_Request request;

  double delay = 0.0, delay_min = 0.0, delay_max = 0.0;
  int buffer = 0;
  int i = 0, j, ITER = 0, N_readings;

  SetValuesOverlap(argc, argv, & ITER, & buffer, & delay_min, & delay_max, & N_readings);

  if (rank == 0) output(argc, argv, ITER, buffer, noOfProcs);

  double t_pure = 0.0, t_comp = 0.0, t_overlap = 0.0;
  double * pure_time, * comp_time, * overall_time;
  pure_time = (double * ) malloc(sizeof(double) * N_readings);
  comp_time = (double * ) malloc(sizeof(double) * N_readings);
  overall_time = (double * ) malloc(sizeof(double) * N_readings);

  double * vec;
  vec = (double * ) malloc(buffer * sizeof(double));

  if (rank == 0) {
    for (j = 0; j < buffer; j++) {
      vec[j] = (char) j;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  for (delay = delay_min; delay < delay_max; delay += 0.005) {
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////PURE///////////////////////////////////////////////////
    for (j = 0; j < N_readings; j++) {
      for (i = 0; i < ITER; i++) {
        if (rank == 1) {
          MPI_Recv(vec, buffer, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, & status);
        } else if (rank == 0) {
          t_pure -= MPI_Wtime();
          MPI_Isend(vec, buffer, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, & request);
          MPI_Wait( & request, & status);
          t_pure += MPI_Wtime();
        }

      }
      t_pure /= ITER;
      ////////////////////////////////////////////////////////////////////////
      //////////////////////Computation///////////////////////////////////////
      MPI_Barrier(MPI_COMM_WORLD);
      for (i = 0; i < ITER; i++) {
        if (rank == 1) {
          MPI_Recv(vec, buffer, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, & status);
        } else if (rank == 0) {
          t_overlap -= MPI_Wtime();
          MPI_Isend(vec, buffer, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, & request);
          t_comp -= MPI_Wtime();
          Computation(vec, delay, buffer);
          t_comp += MPI_Wtime();
          MPI_Wait( & request, & status);
          t_overlap += MPI_Wtime();
        }

      }
      t_overlap /= ITER;
      t_comp /= ITER;
      ///////////////////////////////////////////////////////////////////////////
      if (rank == 0) {
        pure_time[j] = t_pure;
        comp_time[j] = t_comp;
        overall_time[j] = t_overlap;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) reading_output(pure_time, comp_time, overall_time, N_readings);
  }
  free(vec);
  free(pure_time);
  free(comp_time);
  free(overall_time);
  MPI_Finalize();
}

void Irecv(int argc, char * argv[]) {

  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  int noOfProcs;
  MPI_Comm_size(MPI_COMM_WORLD, & noOfProcs);
  MPI_Status status;
  MPI_Request request;

  double delay = 0.0, delay_min = 0.0, delay_max = 0.0;
  int buffer = 0;
  int i = 0, j, ITER = 0, N_readings;

  SetValuesOverlap(argc, argv, & ITER, & buffer, & delay_min, & delay_max, & N_readings);

  if (rank == 0) output(argc, argv, ITER, buffer, noOfProcs);

  double t_pure = 0.0, t_comp = 0.0, t_overlap = 0.0;
  double * pure_time, * comp_time, * overall_time;
  pure_time = (double * ) malloc(sizeof(double) * N_readings);
  comp_time = (double * ) malloc(sizeof(double) * N_readings);
  overall_time = (double * ) malloc(sizeof(double) * N_readings);

  double * vec;
  vec = (double * ) malloc(buffer * sizeof(double));

  if (rank == 0) {
    for (j = 0; j < buffer; j++) {
      vec[j] = (char) j;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  for (delay = delay_min; delay < delay_max; delay += 0.005) {
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////PURE///////////////////////////////////////////////////
    for (j = 0; j < N_readings; j++) {
      for (i = 0; i < ITER; i++) {
        if (rank == 1) {
          MPI_Send(vec, buffer, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else if (rank == 0) {
          t_pure -= MPI_Wtime();
          MPI_Irecv(vec, buffer, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, & request);
          MPI_Wait( & request, & status);
          t_pure += MPI_Wtime();
        }
      }
      t_pure /= ITER;
      ////////////////////////////////////////////////////////////////////////
      //////////////////////Computation///////////////////////////////////////
      MPI_Barrier(MPI_COMM_WORLD);
      for (i = 0; i < ITER; i++) {

        if (rank == 1) {
          MPI_Send(vec, buffer, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else if (rank == 0) {
          t_overlap -= MPI_Wtime();
          MPI_Irecv(vec, buffer, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, & request);
          t_comp -= MPI_Wtime();
          Computation(vec, delay, buffer);
          t_comp += MPI_Wtime();
          MPI_Wait( & request, & status);
          t_overlap += MPI_Wtime();
        }
      }
      t_overlap /= ITER;
      t_comp /= ITER;
      ///////////////////////////////////////////////////////////////////////////
      if (rank == 0) {
        pure_time[j] = t_pure;
        comp_time[j] = t_comp;
        overall_time[j] = t_overlap;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) reading_output(pure_time, comp_time, overall_time, N_readings);
  }
  free(vec);
  free(pure_time);
  free(comp_time);
  free(overall_time);
  MPI_Finalize();
}
