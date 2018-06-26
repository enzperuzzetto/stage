ROOT=../../bcl

CC = gcc
CPPFLAGS = -I. -I$(ROOT)/include
#CFLAGS   = -Wall -Wextra -Werror -pedantic -std=c99 -g
CFLAGS   = -Wall -Wextra -pedantic -std=c99  -O3
LDFLAGS  = -L$(ROOT)/lib
LDLIBS   = -lbcl -lm

BIN= analogy

.PHONY: all
all: $(BIN)

analogy : main.o pyramid.o conversion.o BestMatch.o patchMatch.o stats.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o:%.c
	$(CC) $(CFLAGS) -c  $< -I$(ROOT)/include


.PHONY: clean cleanall
clean:
	$(RM) *.o *~ *.ppm
cleanall: clean
	$(RM) $(BIN)
