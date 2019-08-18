#declare the variables
CC=mpicc
include config.mk
CFLAGS=-c -Wall -g $(MACRO)

MAINOBJ = main.o
ifeq ($(BENCHMARK),pingpong)
	MAINOBJ += functions.o
else ifeq ($(BENCHMARK),pingping)
	MAINOBJ += functions.o
else ifeq ($(BENCHMARK),overlap1)
	MAINOBJ += overlap_functions.o routine.o
else ifeq ($(BENCHMARK),overlap2)
        MAINOBJ +=  routine_overlap2.o overlap_functions.o
endif


.PHONY: run

run: main
ifeq ($(BENCHMARK),pingpong)
	mpiexec -n $(PROCESSES) ./main $(OPTIONS)
else ifeq ($(BENCHMARK),pingping)
	mpiexec -n $(PROCESSES) ./main $(OPTIONS)
else ifeq ($(BENCHMARK),overlap1)
	mpiexec -n $(PROCESSES) ./main $(routine_type) $(OPTIONS)
else ifeq ($(BENCHMARK),overlap2)
	mpiexec -n $(PROCESSES) ./main $(routine_type) $(OPTIONS)	
endif




main: $(MAINOBJ)
ifeq ($(BENCHMARK),pingpong)
	$(CC) $(MAINOBJ) -o main
else ifeq ($(BENCHMARK),pingping)
	$(CC) $(MAINOBJ) -o main
else ifeq ($(BENCHMARK),overlap1)
	$(CC) $(MAINOBJ) -o main
else ifeq ($(BENCHMARK),overlap2)
	$(CC) $(MAINOBJ) -o main
endif




main.o: src/main.c
	$(CC) $(CFLAGS) src/main.c

functions.o: src/functions.c
	$(CC) $(CFLAGS) src/functions.c

overlap_functions.o: src/overlap_functions.c
	$(CC) $(CFLAGS) src/overlap_functions.c

routine.o: src/routine.c
	$(CC) $(CFLAGS) src/routine.c

routine_overlap2.o: src/routine_overlap2.c
	$(CC) $(CFLAGS) src/routine_overlap2.c
clean:
	rm -rf *.o main
	
