#include "includes/routine_overlap2.h"
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "includes/overlap_functions.h"

void Ibcast (int argc, char* argv[])
{

    MPI_Init(NULL, NULL);

    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int noOfProcs=0;
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcs);
    MPI_Status status;
    MPI_Request request;

    int size=0, Max_size=0;
    int i=0,j=0, ITER=0, N_readings=1, c_switch=1000000;
    double c_time_small=0., c_time_large=0.;
    //resolution();
    SetValuesOverlap2 (argc, argv, &ITER, &size, &Max_size, &N_readings, &c_switch, &c_time_small, &c_time_large);

    if(rank ==0) outputOverlap2(argc, argv, ITER, noOfProcs);

    double t_pure =0.0, t_comp=0.0, t_overlap=0.0;
    int root =0; //for collective operations root is the sender

    double *comp;
    comp = (double*) malloc (100000000*sizeof(double));

    double *pure_time, *comp_time, *overall_time;
    pure_time = (double*) malloc (sizeof(double) * N_readings);
    comp_time = (double*) malloc (sizeof(double) * N_readings);
    overall_time = (double*) malloc (sizeof(double) * N_readings);
    int divides =0, k=0;

    for(size=size; size <= Max_size; size*=2)
    {
	if (rank!=-1)
    {
        for (j=0; j<100000000; j++)
        {
            comp[j] = (double) j;
        }
    }
        double *vec;
        vec = (double*) malloc (size*sizeof(double));

        if (rank!=-1)
        {
            for (j=0; j<size; j++)
            {
                vec[j] = (double) j;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        /////////////////// PURE COMMUNICATION /////////////////////////////////////
        for (j=0; j< N_readings; j++)
        {
	    t_pure -= MPI_Wtime();
	    for (i = 0; i< ITER; i++)
            {
                root = i % noOfProcs;

                MPI_Ibcast(vec, size, MPI_DOUBLE, root, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
            }
            t_pure += MPI_Wtime();
            t_pure /=(double) ITER;
        divides = division (comp, t_pure);
	 MPI_Barrier(MPI_COMM_WORLD);
                t_comp -= MPI_Wtime();
	    for (i=0; i<ITER; i++)
            {
                for (k=0; k<divides; k++)
		{
			comp[k] /=1.01;
		}
            }
                t_comp += MPI_Wtime();
            t_comp /= (double)ITER;

            //////////////// OVERLAP REGION   ////////////////////////////////////////////
            MPI_Barrier (MPI_COMM_WORLD);
                t_overlap -= MPI_Wtime();
            for (i=0; i<ITER; i++)
            {
                root = i% noOfProcs;
                MPI_Ibcast(vec, size, MPI_DOUBLE, root, MPI_COMM_WORLD, &request);
                ComputationOverlap2(comp, divides);
                MPI_Wait(&request, &status);
            }
                t_overlap += MPI_Wtime();
            t_overlap /= (double) ITER;

            if(rank ==0) {
                pure_time[j] = t_pure;
                comp_time[j] = t_comp;
                overall_time[j] = t_overlap;
                //printf("Pure time is: %.10f, computation time is: %.10f, overall time is: %.10f, message_size: %ld\n", t_pure, t_comp, t_overlap, size*sizeof(double));
            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        if(rank==0) reading_outputOverlap2(pure_time, comp_time, overall_time, N_readings, size);
        free (vec);
    }
    ////////////////////////////////////////////////////////////////////////////



    free (comp);
    free (pure_time);
    free (comp_time);
    free (overall_time);
    MPI_Finalize();
}

void Igather (int argc, char* argv[])
{
    MPI_Init(NULL, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int noOfProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcs);
    MPI_Status status;
    MPI_Request request;

    int size=0, Max_size=0;
    int i=0,j=0, k=0, divides =0, ITER=0, N_readings=1, c_switch=1000000;

    double c_time_small=0., c_time_large=0.;
    //resolution();
    SetValuesOverlap2 (argc, argv, &ITER, &size, &Max_size, &N_readings, &c_switch, &c_time_small, &c_time_large);

    if(rank ==0) outputOverlap2(argc, argv, ITER, noOfProcs);

    double t_pure =0.0, t_comp=0.0, t_overlap=0.0;
    int root =0; //for collective operations root is the sender

    double *comp;
    comp = (double*) malloc (100000000*sizeof(double));

    double *pure_time, *comp_time, *overall_time;
    pure_time = (double*) malloc (sizeof(double) * N_readings);
    comp_time = (double*) malloc (sizeof(double) * N_readings);
    overall_time = (double*) malloc (sizeof(double) * N_readings);
    //long long int divides =0;


    /////////////////////////////////////////////////////////////////////////////
    /***********************IGATHER******************************************/
    for(size=size; size <= Max_size; size*=2)
    {
      if (rank!=-1)
    {
        for (j=0; j<100000000; j++)
        {
            comp[j] = (double) j;
        }
    }
        int send_size, recv_size; //for scatter and gather
        send_size = size / noOfProcs;
        recv_size = size / noOfProcs;

        double *send_vec;
        send_vec = (double*) malloc (send_size* sizeof(double));

        double *vec; // for scatter this is the sender, for gather this is the receiver
        vec = (double*) malloc (size*sizeof(double));

        //printf ("size out= %s\n", argv[1]);
        if (rank!=-1)
        {
            for (j=0; j<send_size; j++)
            {
                send_vec[j] = (double) j;
            }
        }
        if (rank!=-1)
        {
            for (j=0; j<size; j++)
            {
                vec[j] = (double) j;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        /////////////////// PURE COMMUNICATION /////////////////////////////////////
        for (j=0; j< N_readings; j++)
        {
t_pure -= MPI_Wtime();
            for (i = 0; i< ITER; i++)
            {
                root = i % noOfProcs;

                MPI_Igather (send_vec, send_size, MPI_DOUBLE, vec, recv_size, MPI_DOUBLE, root, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

            }
t_pure += MPI_Wtime();
            t_pure /= ITER;
		divides = division (comp, t_pure);
	 MPI_Barrier(MPI_COMM_WORLD);
                t_comp -= MPI_Wtime();
	    for (i=0; i<ITER; i++)
            {
                for (k=0; k<divides; k++)
		{
			comp[k] /=1.01;
		}
            }
                t_comp += MPI_Wtime();
            t_comp /= (double)ITER;

            ///////////////////////////////////////////////////////////////////////////
            //
            //////////////// OVERLAP REGION   ////////////////////////////////////////////
            MPI_Barrier (MPI_COMM_WORLD);
            t_overlap -= MPI_Wtime();
	    for (i=0; i<ITER; i++)
            {
                root = i % noOfProcs;

                MPI_Igather (send_vec, send_size, MPI_DOUBLE, vec, recv_size, MPI_DOUBLE, root, MPI_COMM_WORLD, &request);


                ComputationOverlap2(comp, divides);

                MPI_Wait(&request, &status);

            }
	    t_overlap += MPI_Wtime();

            t_overlap /= ITER;
            if(rank ==0) {
                pure_time[j] = t_pure;
                comp_time[j] = t_comp;
                overall_time[j] = t_overlap;
                //printf("Pure time is: %.10f, computation time is: %.10f, overall time is: %.10f, message_size: %ld\n", t_pure, t_comp, t_overlap, size*sizeof(double));

            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        if(rank==0) reading_outputOverlap2(pure_time, comp_time, overall_time, N_readings, size);
    }
    ////////////////////////////////////////////////////////////////////////////




    free (pure_time);
    free (comp_time);
    free (overall_time);
    MPI_Finalize();
}


void Ireduce (int argc, char* argv[])
{
    MPI_Init(NULL, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int noOfProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcs);

    MPI_Status status;
    MPI_Request request;

    int size=0, Max_size=0;
    int i=0,j=0, k=0, divides =0, ITER=0, N_readings=1, c_switch=1000000;

    double c_time_small=0., c_time_large=0.;
    //resolution();
    SetValuesOverlap2 (argc, argv, &ITER, &size, &Max_size, &N_readings, &c_switch, &c_time_small, &c_time_large);

    if(rank ==0) outputOverlap2(argc, argv, ITER, noOfProcs);

    double t_pure =0.0, t_comp=0.0, t_overlap=0.0;
    int root =0; //for collective operations root is the sender

    double *comp;
    comp = (double*) malloc (100000000*sizeof(double));

    //printf("Divides:%lld\n", divides);
    //int s_count=0, l_count=0;
    double *pure_time, *comp_time, *overall_time;
    pure_time = (double*) malloc (sizeof(double) * N_readings);
    comp_time = (double*) malloc (sizeof(double) * N_readings);
    overall_time = (double*) malloc (sizeof(double) * N_readings);
    //long long int divides =0;


    /////////////////////////////////////////////////////////////////////////////
    /***********************IREDUCE******************************************/
    for(size=size; size <= Max_size; size*=2)
    {
      if (rank!=-1)
    {
        for (j=0; j<100000000; j++)
        {
            comp[j] = (double) j;
        }
    }
  double *send_vec, *recv_vec;
        send_vec = (double*) malloc (size* sizeof(double));
        recv_vec = (double*) malloc (size* sizeof(double));

        double *vec; // for scatter this is the sender, for gather this is the receiver
        vec = (double*) malloc (size*sizeof(double));

        if (rank==root)
        {
            for (j=0; j<size; j++)
            {
                vec[j] = (double) 1.0;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        /////////////////// PURE COMMUNICATION /////////////////////////////////////
        for (j=0; j< N_readings; j++)
        {
t_pure -= MPI_Wtime();
            for (i = 0; i< ITER; i++)
            {
                root = i % noOfProcs;

                MPI_Ireduce (send_vec, recv_vec, size,MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

            }
t_pure += MPI_Wtime();
            t_pure /= ITER;
		divides = division (comp, t_pure);
	 MPI_Barrier(MPI_COMM_WORLD);
                t_comp -= MPI_Wtime();
	    for (i=0; i<ITER; i++)
            {
                for (k=0; k<divides; k++)
		{
			comp[k] /=1.01;
		}
            }
                t_comp += MPI_Wtime();
            t_comp /= (double)ITER;

            ///////////////////////////////////////////////////////////////////////////
            //
            //////////////// OVERLAP REGION   ////////////////////////////////////////////
            MPI_Barrier (MPI_COMM_WORLD);

                t_overlap -= MPI_Wtime();
            for (i=0; i<ITER; i++)
            {
                root = i % noOfProcs;
                MPI_Ireduce (send_vec, recv_vec, size,MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, &request);

                ComputationOverlap2(comp, divides);

                MPI_Wait(&request, &status);

            }
t_overlap += MPI_Wtime();
            t_overlap /= ITER;
            if(rank ==0) {
                pure_time[j] = t_pure;
                comp_time[j] = t_comp;
                overall_time[j] = t_overlap;
                //printf("Pure time is: %.10f, computation time is: %.10f, overall time is: %.10f, message_size: %ld\n", t_pure, t_comp, t_overlap, size*sizeof(double));

            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        if(rank==0) reading_outputOverlap2(pure_time, comp_time, overall_time, N_readings, size);
    }
    ////////////////////////////////////////////////////////////////////////////





    free (pure_time);
    free (comp_time);
    free (overall_time);
    MPI_Finalize();
}


void Iscatter (int argc, char* argv[])
{
    MPI_Init(NULL, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int noOfProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcs);

    MPI_Status status;
    MPI_Request request;

    int size=0, Max_size=0;
    int i=0,j=0, k=0, divides =0, ITER=0, N_readings=1, c_switch=1000000;

    double c_time_small=0., c_time_large=0.;
    //resolution();
    SetValuesOverlap2 (argc, argv, &ITER, &size, &Max_size, &N_readings, &c_switch, &c_time_small, &c_time_large);

    if(rank ==0) outputOverlap2(argc, argv, ITER, noOfProcs);

    double t_pure =0.0, t_comp=0.0, t_overlap=0.0;
    int root =0; //for collective operations root is the sender
    double *comp;
    comp = (double*) malloc (100000000*sizeof(double));

    //printf("Divides:%lld\n", divides);
    //int s_count=0, l_count=0;
    double *pure_time, *comp_time, *overall_time;
    pure_time = (double*) malloc (sizeof(double) * N_readings);
    comp_time = (double*) malloc (sizeof(double) * N_readings);
    overall_time = (double*) malloc (sizeof(double) * N_readings);
    //long long int divides =0;


    /////////////////////////////////////////////////////////////////////////////
    /***********************ISCATTER******************************************/
    for(size=size; size <= Max_size; size*=2)
    {
      if (rank!=-1)
    {
        for (j=0; j<100000000; j++)
        {
            comp[j] = (double) j;
        }
    }
        int send_size, recv_size; //for scatter and gather
        send_size = size / noOfProcs;
        recv_size = size / noOfProcs;

        double *recv_vec;
        recv_vec = (double*) malloc (recv_size* sizeof(double));

        double *vec; // for scatter this is the sender, for gather this is the receiver
        vec = (double*) malloc (size*sizeof(double));

        if (rank==root)
        {
            for (j=0; j<size; j++)
            {
                vec[j] = (double) j;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        /////////////////// PURE COMMUNICATION /////////////////////////////////////
        for (j=0; j< N_readings; j++)
        {
t_pure -= MPI_Wtime();
            for (i = 0; i< ITER; i++)
            {
                root = i % noOfProcs;

                MPI_Iscatter (vec, send_size, MPI_DOUBLE, recv_vec, recv_size, MPI_DOUBLE, root, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

            }
                t_pure += MPI_Wtime();
            t_pure /= ITER;
	divides = division (comp, t_pure);
	 MPI_Barrier(MPI_COMM_WORLD);
                t_comp -= MPI_Wtime();
	    for (i=0; i<ITER; i++)
            {
                for (k=0; k<divides; k++)
		{
			comp[k] /=1.01;
		}
            }
                t_comp += MPI_Wtime();
            t_comp /= (double)ITER;

            ///////////////////////////////////////////////////////////////////////////
            //
            //////////////// OVERLAP REGION   ////////////////////////////////////////////
            MPI_Barrier (MPI_COMM_WORLD);
 t_overlap -= MPI_Wtime();
		for (i=0; i<ITER; i++)
            {
                root = i % noOfProcs;

                MPI_Iscatter (vec, send_size, MPI_DOUBLE, recv_vec, recv_size, MPI_DOUBLE, root, MPI_COMM_WORLD, &request);

                ComputationOverlap2(comp, divides);

                MPI_Wait(&request, &status);

            }
t_overlap += MPI_Wtime();

            t_overlap /= ITER;
            if(rank ==0) {
                pure_time[j] = t_pure;
                comp_time[j] = t_comp;
                overall_time[j] = t_overlap;
                //printf("Pure time is: %.10f, computation time is: %.10f, overall time is: %.10f, message_size: %ld\n", t_pure, t_comp, t_overlap, size*sizeof(double));

            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        if(rank==0) reading_outputOverlap2(pure_time, comp_time, overall_time, N_readings, size);
    }
    ////////////////////////////////////////////////////////////////////////////


    free (pure_time);
    free (comp_time);
    free (overall_time);
    MPI_Finalize();
}

void Isend (int argc, char* argv[])
{
    MPI_Init(NULL, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int noOfProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcs);
    MPI_Status status;
    MPI_Request request;

    int size=0, Max_size=0;
    int i=0,j=0, k=0, divides =0, ITER=0, N_readings=1, c_switch=1000000;

    double c_time_small=0., c_time_large=0.;
    //resolution();
    SetValuesOverlap2 (argc, argv, &ITER, &size, &Max_size, &N_readings, &c_switch, &c_time_small, &c_time_large);

    if(rank ==0) outputOverlap2(argc, argv, ITER, noOfProcs);
    double *comp;
    comp = (double*) malloc (100000000*sizeof(double));


    //printf("Divides:%lld\n", divides);

    //int s_count=0, l_count=0;
    double t_pure =0.0, t_comp=0.0, t_overlap=0.0;
    double *pure_time, *comp_time, *overall_time;
    pure_time = (double*) malloc (sizeof(double) * N_readings);
    comp_time = (double*) malloc (sizeof(double) * N_readings);
    overall_time = (double*) malloc (sizeof(double) * N_readings);
    //long long int divides =0;


    for(size=size; size <= Max_size; size*=2)
    {
      if (rank!=-1)
    {
        for (j=0; j<100000000; j++)
        {
            comp[j] = (double) j;
        }
    }
        double *vec;
        vec = (double*) malloc (size*sizeof(double));

        if (rank==0)
        {
            for (j=0; j<size; j++)
            {
                vec[j] = (char) j;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        /////////////////////////////////////////////////////////////////////////////
        /////////////////////PURE///////////////////////////////////////////////////
        for (j=0; j< N_readings; j++)
        {
            for (i=0; i<ITER; i++)
            {
                if (rank==1)
                {
                    MPI_Recv(vec, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                }
                else if (rank ==0)
                {
                    t_pure -= MPI_Wtime();
                    MPI_Isend(vec, size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,&request);
                    MPI_Wait (&request, &status);
                    t_pure += MPI_Wtime();
                }

            }
            t_pure /= ITER;
divides = division (comp, t_pure);
	 MPI_Barrier(MPI_COMM_WORLD);
                t_comp -= MPI_Wtime();
	    for (i=0; i<ITER; i++)
            {
                for (k=0; k<divides; k++)
		{
			comp[k] /=1.01;
		}
            }
                t_comp += MPI_Wtime();
            t_comp /= (double)ITER;

            ////////////////////////////////////////////////////////////////////////
            //////////////////////Computation///////////////////////////////////////
            MPI_Barrier (MPI_COMM_WORLD);
            for (i=0; i<ITER; i++)
            {

                if (rank==1)
                {
                    MPI_Recv(vec, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                }

                else if (rank ==0)
                {
                    MPI_Isend(vec, size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,&request);

                    ComputationOverlap2(comp, divides);
                    t_overlap -=MPI_Wtime();
                    MPI_Waitall (1, &request, &status);
                    t_overlap +=MPI_Wtime();
                }

            }
            t_overlap /= ITER;

            ///////////////////////////////////////////////////////////////////////////
            if(rank ==0) {
                pure_time[j] = t_pure;
                comp_time[j] = t_comp;
                overall_time[j] = t_overlap;
                //printf("Pure time is: %.10f, computation time is: %.10f, overall time is: %.10f, message_size: %ld\n", t_pure, t_comp, t_overlap, size*sizeof(double));

            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        if(rank==0) reading_outputOverlap2(pure_time, comp_time, overall_time, N_readings, size);
    }
    //  free (vec);
    free (pure_time);
    free (comp_time);
    free (overall_time);
    MPI_Finalize();
}


void Irecv (int argc, char* argv[])
{
    MPI_Init(NULL, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int noOfProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcs);
    MPI_Status status;
    MPI_Request request;

    int size=0, Max_size=0;
    int i=0,j=0, k=0, divides =0, ITER=0, N_readings=1, c_switch=1000000;

    double c_time_small=0., c_time_large=0.;
    //resolution();
    SetValuesOverlap2 (argc, argv, &ITER, &size, &Max_size, &N_readings, &c_switch, &c_time_small, &c_time_large);

    if(rank ==0) outputOverlap2(argc, argv, ITER, noOfProcs);
    double *comp;
    comp = (double*) malloc (100000000*sizeof(double));



    //printf("Divides:%lld\n", divides);

    //int s_count=0, l_count=0;
    double t_pure =0.0, t_comp=0.0, t_overlap=0.0;
    double *pure_time, *comp_time, *overall_time;
    pure_time = (double*) malloc (sizeof(double) * N_readings);
    comp_time = (double*) malloc (sizeof(double) * N_readings);
    overall_time = (double*) malloc (sizeof(double) * N_readings);
    //long long int divides =0;


    for(size=size; size <= Max_size; size*=2)
    {
      if (rank!=-1)
    {
        for (j=0; j<100000000; j++)
        {
            comp[j] = (double) j;
        }
    }
        double *vec;
        vec = (double*) malloc (size*sizeof(double));

        if (rank==0)
        {
            for (j=0; j<size; j++)
            {
                vec[j] = (char) j;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        /////////////////////////////////////////////////////////////////////////////
        /////////////////////PURE///////////////////////////////////////////////////
        for (j=0; j< N_readings; j++)
        {
            for (i=0; i<ITER; i++)
            {
                if (rank==1)
                {
                    MPI_Send(vec, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                }
                else if (rank ==0)
                {
                    t_pure -= MPI_Wtime();
                    MPI_Irecv(vec, size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,&request);
                    MPI_Wait (&request, &status);
                    t_pure += MPI_Wtime();
                }

            }
            t_pure /= ITER;
divides = division (comp, t_pure);
	 MPI_Barrier(MPI_COMM_WORLD);
                t_comp -= MPI_Wtime();
	    for (i=0; i<ITER; i++)
            {
                for (k=0; k<divides; k++)
		{
			comp[k] /=1.01;
		}
            }
                t_comp += MPI_Wtime();
            t_comp /= (double)ITER;

            ////////////////////////////////////////////////////////////////////////
            //////////////////////Computation///////////////////////////////////////
            MPI_Barrier (MPI_COMM_WORLD);
            for (i=0; i<ITER; i++)
            {

                if (rank==1)
                {
                    MPI_Send(vec, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                }

                else if (rank ==0)
                {
                    t_overlap -=MPI_Wtime();
                    MPI_Irecv(vec, size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,&request);

                    ComputationOverlap2(comp, divides);

                    MPI_Wait (&request, &status);
                    t_overlap +=MPI_Wtime();
                }

            }
            t_overlap /= ITER;

            ///////////////////////////////////////////////////////////////////////////
            if(rank ==0) {
                pure_time[j] = t_pure;
                comp_time[j] = t_comp;
                overall_time[j] = t_overlap;
                //printf("Pure time is: %.10f, computation time is: %.10f, overall time is: %.10f, message_size: %ld\n", t_pure, t_comp, t_overlap, size*sizeof(double));

            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        if(rank==0) reading_outputOverlap2(pure_time, comp_time, overall_time, N_readings, size);
    }
    //  free (vec);
    free (pure_time);
    free (comp_time);
    free (overall_time);
    MPI_Finalize();
}
