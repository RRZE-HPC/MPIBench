//for compilation:

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "includes/routine.h"


int main (int argc, char* argv[])
{


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

}


