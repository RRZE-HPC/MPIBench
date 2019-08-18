***
# MPIBench
MPIBench is a software that contains benchmark suites used for investigating MPI communication behavior. Benchmarks probe point-to-point communication and the communication/computation overlap. It is written in C language and implemented on top of MPI standard routines. Each routine is implemented in our way. It is portable and flexible and works with several MPI implemetations. 

It contains following benchmarks:
1. PingPong
2. PingPing
3. Overlapping benchmarks:  
..* Pair-wise routines  
..* Collective routines  
