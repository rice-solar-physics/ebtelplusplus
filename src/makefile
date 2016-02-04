UNAME = $(shell uname)
CC=gcc
CFLAGS=-Wall -g -std=c99
ifeq ($(UNAME), Darwin)
LFLAGS=-lxml2 -lz -liconv -lm
IFLAGS=-I /opt/local/include/libxml2
endif
ifeq ($(UNAME), Linux)
LFLAGS=-lxml2 -lz -lm
IFLAGS=-I /usr/include/libxml2
endif

all: set_xml_flags ebtel-2fl

set_xml_flags:
	@if test -z "$$LFLAGS"; then echo "Default values for LFLAGS"; else echo "LFLAGS = $$LFLAGS"; fi
	@if test -z "$$IFLAGS"; then echo "Default values for IFLAGS"; else echo "IFLAGS = $$IFLAGS"; fi

ebtel-2fl: ebtel-2fl_functions_loop.o ebtel-2fl_functions_param.o ebtel-2fl_functions_util.o ebtel-2fl_functions_solvers.o ebtel-2fl_functions_heating.o ebtel-2fl_main.o 
	$(CC) $(CFLAGS) ebtel-2fl_functions_loop.o ebtel-2fl_functions_param.o ebtel-2fl_functions_util.o ebtel-2fl_functions_solvers.o ebtel-2fl_functions_heating.o ebtel-2fl_main.o -o ../bin/ebtel $(LFLAGS)

ebtel-2fl_functions_loop.o: ebtel-2fl_functions_loop.c
	$(CC) -c  $(CFLAGS) $(IFLAGS) ebtel-2fl_functions_loop.c 

ebtel-2fl_functions_param.o: ebtel-2fl_functions_param.c
	$(CC) -c  $(CFLAGS) $(IFLAGS) ebtel-2fl_functions_param.c

ebtel-2fl_functions_util.o: ebtel-2fl_functions_util.c
	$(CC) -c $(CFLAGS) $(IFLAGS) ebtel-2fl_functions_util.c

ebtel-2fl_functions_solvers.o: ebtel-2fl_functions_solvers.c
	$(CC) -c  $(CFLAGS) $(IFLAGS) ebtel-2fl_functions_solvers.c

ebtel-2fl_functions_heating.o: ebtel-2fl_functions_heating.c
	$(CC) -c  $(CFLAGS) $(IFLAGS) ebtel-2fl_functions_heating.c

ebtel-2fl_main.o: ebtel-2fl_main.c
	$(CC) -c $(CFLAGS) $(IFLAGS) ebtel-2fl_main.c

clean: 
	rm *.o ../bin/ebtel
