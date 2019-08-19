#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>

#include <stdlib.h>

#include <unistd.h>

#include <string.h>

#include <mpi.h>

#include "includes/functions.h"

#define MIN(x, y)((x) < (y) ? (x) : (y))
#define MAX(x, y)((x) > (y) ? (x) : (y))

void swap(char ** a, char ** b) {
  char * temp = * a;
  * a = * b;
  * b = temp;
}

void SetValues(int argc, char * argv[], int * ITER, int * Max_size, int * size, int * N_readings) {
  int i;
  for (i = 1; i < argc; ++i) {
    if ((!strcmp(argv[i], "-Max_size"))) {
      * Max_size = atoi(argv[++i]);
    }
    if ((!strcmp(argv[i], "-ITER"))) {
      * ITER = atoi(argv[++i]);
    }
    if ((!strcmp(argv[i], "-size"))) {
      * size = atoi(argv[++i]);
    }
    if ((!strcmp(argv[i], "-N_readings"))) {
      * N_readings = atoi(argv[++i]);
    }
  }
}

void print(double vec_time[], double vec_bwidth[], int N_readings, int size) {
  int i;
  double avg_time, min_time = vec_time[0], max_time = vec_time[0], avg_bwidth, min_bwidht = vec_bwidth[0], max_bwidth = vec_bwidth[0], sum_time = vec_time[0], sum_bwidth = vec_bwidth[0];

  for (i = 1; i < N_readings; i++) {
    sum_time += vec_time[i];
    sum_bwidth += vec_bwidth[i];
    min_time = MIN(vec_time[i], min_time);
    max_time = MAX(vec_time[i], max_time);
    min_bwidht = MIN(vec_bwidth[i], min_bwidht);
    max_bwidth = MAX(vec_bwidth[i], max_bwidth);
}
  avg_time = sum_time / N_readings;
  avg_bwidth = sum_bwidth / N_readings;
  if (size < 16000000)
    printf("\t%d\t\t%.6f\t\t  %.6f\t  %.6f\t\t    %.2f\t\t\t      %.2f\t\t     %.2f\n", size, avg_time, min_time, max_time, avg_bwidth, min_bwidht, max_bwidth);
  else
    printf("\t%d\t%.6f\t\t  %.6f\t  %.6f\t\t    %.2f\t\t\t      %.2f\t\t     %.2f\n", size, avg_time, min_time, max_time, avg_bwidth, min_bwidht, max_bwidth);
}

void pingping(int argc, char * argv[]) {

  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  int noOfProcs;
  MPI_Comm_size(MPI_COMM_WORLD, & noOfProcs);
  MPI_Status status;
  MPI_Request request;

  int j, i, ITER, size, Max_size, N_readings;
  SetValues(argc, argv, & ITER, & Max_size, & size, & N_readings);
  if (rank == 0) {
    double tick = MPI_Wtick();
    printf("************************************************\n");
    printf("Routine \t\t PingPing \nNo of outer iterations:  %d \nNo of runs:\t\t %d \nMin message size: \t %d \nMax message size: \t %d \nNo of Processes: \t %d \nTime Resolution (s): \t %.15f\n", ITER, N_readings, size, Max_size, noOfProcs, tick);
    printf("************************************************\n\n");
    printf("Message size [B]\tAverage time[μs]\tMin time[μs]\tMax time[μs]\t\tAverage Bandwidth[MB/s]\t\tMin Bandwidth[MB/s]\tMax Bandwidth[MB/s] \n");
  }

  double time_start, time_end, time_reported, BWIDTH;
  double * vec_time, * vec_bwidth;
  vec_time = (double * ) malloc(sizeof(double) * N_readings);
  vec_bwidth = (double * ) malloc(sizeof(double) * N_readings);

  for (size = size; size <= Max_size; size *= 2) {
    char * vec, * vec1;
    vec = (char * ) malloc(size * sizeof(char));
    vec1 = (char * ) malloc(size * sizeof(char));
    if (rank < (noOfProcs / 2)) {
      for (j = 0; j < size; j++) {
        vec[j] = (char) j;
      }
    }
    if (rank >= (noOfProcs / 2)) {
      for (j = 0; j < size; j++) {
        vec1[j] = (char) j;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (j = 0; j < N_readings; j++) {
      if (rank == 0) time_start = MPI_Wtime();

      for (i = 0; i < ITER; i++){
        if (rank < (noOfProcs / 2)) {
          MPI_Isend(vec, size, MPI_CHAR, (noOfProcs / 2) + rank, rank, MPI_COMM_WORLD, & request);
          MPI_Recv(vec1, size, MPI_CHAR, rank + (noOfProcs / 2), rank + (noOfProcs / 2), MPI_COMM_WORLD, & status);
          MPI_Wait( & request, & status);
        } else if (rank >= (noOfProcs / 2)) {
          MPI_Isend(vec1, size, MPI_CHAR, rank - (noOfProcs / 2), rank, MPI_COMM_WORLD, & request);
          MPI_Recv(vec, size, MPI_CHAR, rank - (noOfProcs / 2), rank - (noOfProcs / 2), MPI_COMM_WORLD, & status);
          MPI_Wait( & request, & status);
        }
        swap( & vec, & vec1);
				#ifdef BARRIER
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
      }
      if (rank == 0) {
        time_end = MPI_Wtime();

        time_reported = (time_end - time_start) * 1000000.0 / (double) ITER;
        BWIDTH = (double)(noOfProcs / 2) * (double) sizeof(char) * (double) size * (double) ITER * 2.0 / (time_end - time_start) / 1000000.0;
        vec_time[j] = time_reported;
        vec_bwidth[j] = BWIDTH;
      }
    }
    if (rank == 0) {
      print(vec_time, vec_bwidth, N_readings, size);
    }

    free(vec);
    free(vec1);
  }
  free(vec_time);
  free(vec_bwidth);
  MPI_Finalize();
}

void pingpong(int argc, char * argv[]) {

  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  int noOfProcs;
  MPI_Comm_size(MPI_COMM_WORLD, & noOfProcs);
  MPI_Status status;

  int j, i, ITER, size, Max_size, N_readings;
  SetValues(argc, argv, & ITER, & Max_size, & size, & N_readings);
  if (rank == 0) {
    double tick = MPI_Wtick();
    printf("************************************************\n");
    printf("Routine \t\t PingPong \nNo of outer iterations:  %d \nNo of runs:\t\t %d \nMin message size: \t %d \nMax message size: \t %d \nNo of Processes: \t %d \nTime Resolution (s): \t %.15f\n", ITER, N_readings, size, Max_size, noOfProcs, tick);
    printf("************************************************\n\n");
    printf("Message size [B]\tAverage time[μs]\tMin time[μs]\tMax time[μs]\t\tAverage Bandwidth[MB/s]\t\tMin Bandwidth[MB/s]\tMax Bandwidth[MB/s] \n");
  }

  double time_start, time_end, time_reported, BWIDTH;
  double * vec_time, * vec_bwidth;
  vec_time = (double * ) malloc(sizeof(double) * N_readings);
  vec_bwidth = (double * ) malloc(sizeof(double) * N_readings);

  for (size = size; size <= Max_size; size *= 2) {
    char * vec, * vec1;
    vec = (char * ) malloc(size * sizeof(char));
    vec1 = (char * ) malloc(size * sizeof(char));
    if (rank < (noOfProcs / 2)) {
      for (j = 0; j < size; j++) {
        vec[j] = (char) j;
      } 
    }
    if (rank >= (noOfProcs / 2)) {
      for (j = 0; j < size; j++) {
        vec1[j] = (char) j;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (j = 0; j < N_readings; j++) {
      if (rank == 0) time_start = MPI_Wtime();

      for (i = 0; i < ITER; i++) {
        if (rank < (noOfProcs / 2)) {
          MPI_Send((char * ) vec, size, MPI_CHAR, (noOfProcs / 2) + rank, rank, MPI_COMM_WORLD);
          MPI_Recv((char * ) vec1, size, MPI_CHAR, (noOfProcs / 2) + rank, (noOfProcs / 2) + rank, MPI_COMM_WORLD, & status);
        } else if (rank >= (noOfProcs / 2)) {
          MPI_Recv((char * ) vec, size, MPI_CHAR, rank - (noOfProcs / 2), rank - (noOfProcs / 2), MPI_COMM_WORLD, & status);
          MPI_Send((char * ) vec1, size, MPI_CHAR, rank - (noOfProcs / 2), rank, MPI_COMM_WORLD);
        }
        swap( & vec, & vec1);
  			#ifdef BARRIER
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
      }
      if (rank == 0) {
        time_end = MPI_Wtime();
        time_reported = (time_end - time_start) / 2.0 * 1000000.0 / (double) ITER;
        BWIDTH = (double) sizeof(char) * (double) size * (double) ITER * 2.0 / (time_end - time_start) / 1000000.0;
        vec_time[j] = time_reported;
        vec_bwidth[j] = BWIDTH;
      }
    }
    if (rank == 0) {
      print(vec_time, vec_bwidth, N_readings, size);
    }
		free(vec);
    free(vec1);
  }
  free(vec_time);
  free(vec_bwidth);
  MPI_Finalize();
}

#endif
