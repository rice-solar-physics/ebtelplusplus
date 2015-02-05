CC=gcc
CFLAGS=-Wall -g -std=c99

all: ebtel-2fl

ebtel-2fl: ebtel-2fl_functions_loop.o ebtel-2fl_functions_param.o ebtel-2fl_functions_util.o ebtel-2fl_functions_solvers.o ebtel-2fl_functions_heating.o ebtel-2fl_main.o 
	$(CC) $(CFLAGS) ebtel-2fl_functions_loop.o ebtel-2fl_functions_param.o ebtel-2fl_functions_util.o ebtel-2fl_functions_solvers.o ebtel-2fl_functions_heating.o ebtel-2fl_main.o -o ebtel-2fl -lm

ebtel-2fl_functions_loop.o: ebtel-2fl_functions_loop.c
	$(CC) -c  $(CFLAGS) ebtel-2fl_functions_loop.c 

ebtel-2fl_functions_param.o: ebtel-2fl_functions_param.c
	$(CC) -c  $(CFLAGS) ebtel-2fl_functions_param.c

ebtel-2fl_functions_util.o: ebtel-2fl_functions_util.c
	$(CC) -c $(CFLAGS) ebtel-2fl_functions_util.c

ebtel-2fl_functions_solvers.o: ebtel-2fl_functions_solvers.c
	$(CC) -c  $(CFLAGS) ebtel-2fl_functions_solvers.c

ebtel-2fl_functions_heating.o: ebtel-2fl_functions_heating.c
	$(CC) -c  $(CFLAGS) ebtel-2fl_functions_heating.c

ebtel-2fl_main.o: ebtel-2fl_main.c
	$(CC) -c $(CFLAGS) ebtel-2fl_main.c

clean: 
	rm ebtel-2fl*.o ebtel-2fl
