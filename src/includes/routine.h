#ifndef ROUTINE_H
#define ROUTINE_H

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "overlap_functions.h"

void Ibcast (int argc, char* argv[]);

void Igather (int argc, char* argv[]);

void Ireduce (int argc, char* argv[]);

void Iscatter (int argc, char* argv[]);

void Isend (int argc, char* argv[]);

void Irecv (int argc, char* argv[]);

#endif // ROUTINE_H
