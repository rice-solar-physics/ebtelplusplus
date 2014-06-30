CC=gcc
CFLAGS=-Wall -g -std=c99

all: ebtel

ebtel: ebtel_functions_loop.o ebtel_functions_param.o ebtel_functions_util.o ebtel_functions_solvers.o ebtel_main.o 
	$(CC) $(CFLAGS) ebtel_functions_loop.o ebtel_functions_param.o ebtel_functions_util.o ebtel_functions_solvers.o ebtel_main.o -o ebtel -lm

ebtel_functions_loop.o: ebtel_functions_loop.c
	$(CC) -c  $(CFLAGS) ebtel_functions_loop.c 

ebtel_functions_param.o: ebtel_functions_param.c
	$(CC) -c  $(CFLAGS) ebtel_functions_param.c

ebtel_functions_util.o: ebtel_functions_util.c
	$(CC) -c $(CFLAGS) ebtel_functions_util.c

ebtel_functions_solvers.o: ebtel_functions_solvers.c
	$(CC) -c  $(CFLAGS) ebtel_functions_solvers.c

ebtel_main.o: ebtel_main.c
	$(CC) -c $(CFLAGS) ebtel_main.c

clean: 
	rm *.o ebtel
