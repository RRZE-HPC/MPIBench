#include "includes/overlap_functions.h"
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "includes/routine.h"
#include "includes/routine_overlap2.h"
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

void output (int argc, char* argv[], int ITER, int buffer, int noOfProcs)
{
    double tick = MPI_Wtick();
    printf ("************************************************\n");
    printf ("Requested routine: \t %s \nNo of outer iterations:  %d \nBuffer size[B]: \t %ld \nNo of Processes: \t %d \nTime Resolution (s) \t %.15f\n", argv[1], ITER, buffer*sizeof(double), noOfProcs, tick);
    printf ("************************************************\n\n");
    printf("\tCommunication Time [s]\t\t\t\t\t Computational Time [s]\t\t\t\t\t Overall Time [s]\n");
    printf ("Average \t   Min \t\t   Max \t\t\t\t   Average \t\t\t\t\t    Average \t   Min \t\t Max\n\n");
}

void outputOverlap2 (int argc, char* argv[], int ITER, int noOfProcs)
{
    double tick = MPI_Wtick();
    printf ("************************************************\n");
    printf ("Requested routine: \t %s \nNo of outer iterations:  %d  \nNo of Processes: \t %d \nTime Resolution (s) \t %.15f\n", argv[1], ITER, noOfProcs, tick);
    printf ("************************************************\n\n");
    printf("Message \t\tCommunication Time [s]\t\t\t\t Computational Time [s]\t\t\t\t\t Overall Time [s]\n");
    printf ("size [B]\t Average \t   Min \t\t   Max \t\t\t\t   Average \t\t\t\t    Average \t   Min \t\t Max\n\n");
}

void reading_output (double pure_time[], double comp_time[], double overall_time[], int N_readings)
{
    int i;
    double avg_ptime, min_ptime=pure_time[0], max_ptime=pure_time[0], avg_ctime, sum_ptime=pure_time[0], sum_ctime=comp_time[0], avg_otime, min_otime=overall_time[0], max_otime=overall_time[0], sum_otime=overall_time[0];
    /*
    for (i=0; i< N_readings; i++)
    {
      printf ("\n Total Time is %.10f  and bandwith is %.2f \n", vec_time[i], vec_bwidth[i]);
    }
    */
    for (i=1; i< N_readings; i++)
    {
        sum_ptime += pure_time[i];
        sum_ctime += comp_time[i];
        sum_otime += overall_time[i];

        min_ptime = MIN (pure_time[i], min_ptime);
        max_ptime = MAX (pure_time[i], max_ptime);

        min_otime = MIN (overall_time[i], min_otime);
        max_otime = MAX (overall_time[i], max_otime);
        //printf ("\n Total Time is %.10f  and bandwith is %.2f \n", vec_time[i], vec_bwidth[i]);
    }
    avg_ptime = sum_ptime / N_readings;
    avg_ctime = sum_ctime / N_readings;
    avg_otime = sum_otime / N_readings;
    printf("%.10f\t %.10f\t %.10f\t\t\t %.10f\t\t\t\t\t %.10f\t %.10f\t %.10f\t\n", avg_ptime, min_ptime, max_ptime, avg_ctime,  avg_otime, min_otime, max_otime );
    /*
    if (size <16000000)
    printf ("\t%d\t\t%.6f\t\t  %.6f\t  %.6f\t\t    %.2f\t\t\t      %.2f\t\t     %.2f\n",size, avg_time, min_time, max_time, avg_bwidth, min_bwidht, max_bwidth);
    else
    printf ("\t%d\t%.6f\t\t  %.6f\t  %.6f\t\t    %.2f\t\t\t      %.2f\t\t     %.2f\n",size, avg_time, min_time, max_time, avg_bwidth, min_bwidht, max_bwidth);
    */

}


void Computation (double arr[], double delay, int buffer)
{
    double t1 = 0.0, t2 = 0.0;
    int j;

    t1 = MPI_Wtime();
    for (j=1; j<=buffer; j++)
    {
        t2 = MPI_Wtime();
        if ( t2 -t1 < delay)
        {
            arr[j] = arr[j] / j;
        }
        else {
            break;
        }
    }
}


void SetValuesOverlap (int argc, char* argv[], int *ITER, int *buffer, double *delay_min, double *delay_max, int *N_readings)
{
    int i;
    for (i =1; i< argc ; ++i)
    {
        if((!strcmp(argv[i], "-ITER"))) {
            *ITER = atoi (argv[++i]);
        }
        if((!strcmp(argv[i], "-buffer"))) {
            *buffer = atoi (argv[++i]);
        }
        if((!strcmp(argv[i], "-delay_min"))) {
            *delay_min = atof (argv[++i]);
        }
        if((!strcmp(argv[i], "-delay_max"))) {
            *delay_max = atof (argv[++i]);
        }
        if((!strcmp(argv[i], "-N_readings"))) {
            *N_readings = atoi (argv[++i]);
        }
    }
}

void SetValuesOverlap2 (int argc, char *argv[], int* ITER, int* size, int* Max_size, int *N_readings, int* c_switch, double* c_time_small, double* c_time_large) //c_switch to change computation time between smaller and larger message sizes
{
    int i;
    //printf("ARGC %d\n", argc);
    for (i =1; i< argc ; ++i)
    {
        if((!strcmp(argv[i], "-ITER"))) {
            *ITER = atoi (argv[++i]);
        }
        if((!strcmp(argv[i], "-size"))) {
            *size = atoi (argv[++i]);
        }
        if((!strcmp(argv[i], "-Max_size"))) {
            *Max_size = atoi (argv[++i]);
        }
        if((!strcmp(argv[i], "-N_readings"))) {
            *N_readings = atoi (argv[++i]);
        }
        if((!strcmp(argv[i], "-c_switch"))) {
            *c_switch = atoi (argv[++i]);
        }
        if((!strcmp(argv[i], "-c_time_small"))) {
            *c_time_small = atof (argv[++i]);
        }
        if((!strcmp(argv[i], "-c_time_large"))) {
            *c_time_large = atof (argv[++i]);
        }
    }

}

int division(double comp[], double c_time)
{
  int j=0;
  double  t1=0.0, t2=0.0,t=0.0;
  t1= MPI_Wtime();
  for (j=0; j<100000000; j++)
  {
      comp[j] /= 1.23;
  }
  t2 = MPI_Wtime();
  t=(t2-t1)/100000000.0;
  //printf ("Time for 1 divide: %.15f\n", t);
  double divides = c_time / t;
  //printf("c_time:%.10f : Divides:   %.15f\n",c_time, divides);
  return divides;
}



void reading_outputOverlap2 (double pure_time[], double comp_time[], double overall_time[], int N_readings, int size)
{
    int i;
    double avg_ptime, min_ptime=pure_time[0], max_ptime=pure_time[0], avg_ctime, sum_ptime=pure_time[0], sum_ctime=comp_time[0], avg_otime, min_otime=overall_time[0], max_otime=overall_time[0], sum_otime=overall_time[0];
    for (i=1; i< N_readings; i++)
    {
        sum_ptime += pure_time[i];
        sum_ctime += comp_time[i];
        sum_otime += overall_time[i];

        min_ptime = MIN (pure_time[i], min_ptime);
        max_ptime = MAX (pure_time[i], max_ptime);

        min_otime = MIN (overall_time[i], min_otime);
        max_otime = MAX (overall_time[i], max_otime);
        //printf ("\n Total Time is %.10f  and bandwith is %.2f \n", vec_time[i], vec_bwidth[i]);
    }
    avg_ptime = sum_ptime / N_readings;
    avg_ctime = sum_ctime / N_readings;
    avg_otime = sum_otime / N_readings;
    printf("%ld \t\t %.10f\t %.10f\t %.10f\t\t\t %.10f\t\t\t\t %.10f\t %.10f\t %.10f\t\n",size*sizeof(double), avg_ptime, min_ptime, max_ptime, avg_ctime,  avg_otime, min_otime, max_otime );
}


void ComputationOverlap2 (double arr[], int divides)
{
    int j;
    for (j=1; j<=divides; j++)
    {
            arr[j] = arr[j] / 1.23;
//		if (j>=divides) printf("J is %d\n", j);
    }
}
