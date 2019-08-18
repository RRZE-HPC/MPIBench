#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "includes/functions.h"
#include "includes/routine.h"
#include "includes/overlap_functions.h"

#include "includes/routine_overlap2.h"



int main (int argc, char* argv[])
{
	#ifdef PINGPONG
	pingpong (argc, argv);
	#endif

	#ifdef PINGPING
	pingping (argc, argv);
	#endif

	#ifdef OVERLAP
	if (strcmp(argv[1] , "Ibcast") == 0) {
        Ibcast (argc, argv);
    }
    else if (strcmp(argv[1] , "Igather") == 0) {
        Igather (argc, argv);
    }
    else if (strcmp(argv[1] , "Ireduce") == 0) {
        Ireduce (argc, argv);
    }
    else if (strcmp(argv[1] , "Iscatter") == 0) {
        Iscatter (argc, argv);
    }
    else if (strcmp(argv[1] , "Irecv") == 0) {
        Irecv (argc, argv);
    }
    else if (strcmp(argv[1] , "Isend") == 0) {
        Isend (argc, argv);
    }
	#endif

}
