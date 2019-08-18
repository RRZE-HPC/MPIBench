#Feature options


BENCHMARK =	overlap1
MACRO =		-DOVERLAP
routine_type =  Ibcast
PROCESSES = 	2

OPTIONS =	-ITER 10
OPTIONS +=      -N_readings 3
ifeq ($(BENCHMARK),pingpong)
OPTIONS +=	-Max_size 140000
OPTIONS +=	-size 1

else ifeq ($(BENCHMARK),pingping)
OPTIONS +=      -Max_size 140000
OPTIONS +=      -size 1

else ifeq ($(BENCHMARK),overlap1)
OPTIONS +=	-buffer 1000
OPTIONS +=	-delay_min 0.01
OPTIONS +=	-delay_max 0.05

else ifeq ($(BENCHMARK),overlap2)
OPTIONS +=      -Max_size 140000
OPTIONS +=      -size 1

endif
