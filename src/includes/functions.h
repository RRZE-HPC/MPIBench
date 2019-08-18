#ifndef FUNCTIONS
#define FUNCTIONS
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

void pingpong (int argc, char* argv[]);
void pingping (int argc, char* argv[]);
void swap (char** a, char** b);
void SetValues (int argc, char* argv[], int *ITER, int *Max_size, int *size, int *N_readings);
void print (double vec_time[], double vec_bwidth[], int N_readings,int size);
#endif


