#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "routine.h"
#include "routine_overlap2.h"

void output (int argc, char* argv[], int ITER, int buffer, int noOfProcs);

void outputOverlap2 (int argc, char* argv[], int ITER, int noOfProcs);

void resolution ();

void reading_output (double pure_time[], double comp_time[], double overall_time[], int N_readings);

void Computation (double arr[], double delay, int buffer);

void SetValuesOverlap (int argc, char* argv[], int *ITER, int *buffer, double *delay_min, double *delay_max, int *N_readings);

void SetValuesOverlap2 (int argc, char *argv[], int* ITER, int* size, int* Max_size, int *N_readings, int *c_switch, double* c_time_small, double* c_time_large);

int division (double comp[], double c_time);

void reading_outputOverlap2 (double pure_time[], double comp_time[], double overall_time[], int N_readings, int size);

void ComputationOverlap2 (double arr[], int divides);

#endif // FUNCTIONS_H
