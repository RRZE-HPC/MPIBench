***
# MPIBench
MPIBench is a software that contains benchmark suites used for investigating MPI communication behavior. These benchmarks probe point-to-point communication and the communication/computation overlap. It is written in C language and implemented on top of MPI standard routines. We implemented each routine according to our use case. It is portable and flexible and works with several MPI implemetations. 

It contains following benchmarks:
1. PingPong
2. PingPing
3. Overlapping benchmarks: 
* Pair-wise routines  
* Collective routines  

## Build
1. Configure the application, timing, and problem specific parameters in  `config.mk`.
```C
#Feature options
					
BENCHMARK =     overlap2 //pingpong, pingping, overlap1, overlap2
MACRO 	  =    -DOVERLAP //-DOVERLAP, -DPINGPNG, -DPINGPONG
routine_type =  Ibcast //Ibcast, Igather, Ireduce, Isactter, Isend, Irecv
PROCESSES =     2
	
OPTIONS =       -ITER 10 //outer iterations for timing calculations
OPTIONS +=      -N_readings 3 // No of times a benchmark is run (for averaging)
ifeq ($(BENCHMARK),pingpong)
OPTIONS +=      -Max_size 140000 //Maximum data size being transferred between two processes
OPTIONS +=      -size 1 //smallest data sent between processes

else ifeq ($(BENCHMARK),pingping)
OPTIONS +=      -Max_size 140000 //Maximum data size being transferred between two processes
OPTIONS +=      -size 1 //smallest data sent between processes
			
else ifeq ($(BENCHMARK),overlap1)
OPTIONS +=      -buffer 1000 //communication buffer size
OPTIONS +=      -delay_min 0.01 //minimum computation time
OPTIONS +=      -delay_max 0.05 //maximum computation time
			
else ifeq ($(BENCHMARK),overlap2)
OPTIONS +=      -Max_size 140000 //Maximum data size sent between processes
OPTIONS +=      -size 1 // smallest data sent between processes
					
endif
```
The variable **BENCHMARK** takes the the name of the benchmark. For overlapping benchmarks, the variable **routine_type** takes the routine whose overlap needs to be investigated. The variable **PROCESSES** takes the number of MPI processes for the benchmark run. 

2. To build and run:  
``` 
make
```
It will automatically run the benchmark according the configuration file and provides the output.

3. Clean up with:
``` 
make clean
```
For pinning MPI threads, either [Intel](https://software.intel.com/en-us/mpi-developer-reference-linux-environment-variables-for-process-pinning) and [Likwid](https://github.com/RRZE-HPC/likwid/wiki/Likwid-Pin) can be used. 

