PEDANTIC_PARANOID_FREAK =       -g -lm -O0 -Wall -Wshadow -Wcast-align \
				-Waggregate-return -Wstrict-prototypes \
				-Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions -Werror
CC = gcc
CPPFLAGS = -DDEBUG
LDFLAGS = -g -lm

include Makefile.config

ifdef FILTER_SIZE
CONFIG_FLAGS += -DFILTER_SIZ=$(FILTER_SIZE)
endif

CFLAGS = $(PEDANTIC_PARANOID_FREAK) $(CONFIG_FLAGS)

PROG = main

.PHONY: all
all : $(PROG)

OBJS = main.o imgoper.o convol.o edge.o filter.o hist.o matrix.o median.o pgmfil.o arr.o lst.o kmeans.o offfil.o
main : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(PROG) $(LDFLAGS)

clean :
	@rm -f $(PROG) *.o

cleanall : clean
	@rm -f $(PROG)
