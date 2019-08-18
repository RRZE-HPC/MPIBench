***
# MPIBench
MPIBench is a software that contains benchmark suites used for investigating MPI communication behavior. Benchmarks probe point-to-point communication and the communication/computation overlap. It is written in C language and implemented on top of MPI standard routines. Each routine is implemented in our way. It is portable and flexible and works with several MPI implemetations. 

It contains following benchmarks:
1. PingPong
2. PingPing
3. Overlapping benchmarks: 
* Pair-wise routines  
* Collective routines  

##Build
1. Configure the application, timing, and problem specific parameters in  `config.mk`.
```C
#Feature options
					
					BENCHMARK =     overlap2
					MACRO =         -DOVERLAP
					routine_type =  Ibcast
					PROCESSES =     2
					
					OPTIONS =       -ITER 10
					OPTIONS +=      -N_readings 3
					ifeq ($(BENCHMARK),pingpong)
					OPTIONS +=      -Max_size 140000
					OPTIONS +=      -size 1
					
					else ifeq ($(BENCHMARK),pingping)
					OPTIONS +=      -Max_size 140000
					OPTIONS +=      -size 1
					
					else ifeq ($(BENCHMARK),overlap1)
					OPTIONS +=      -buffer 1000
					OPTIONS +=      -delay_min 0.01
					OPTIONS +=      -delay_max 0.05
					
					else ifeq ($(BENCHMARK),overlap2)
					OPTIONS +=      -Max_size 140000
					OPTIONS +=      -size 1
					
					endif
```
