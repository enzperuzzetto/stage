ROOT=../bcl

CC = gcc
CPPFLAGS = -I. -I$(ROOT)/include 
#CFLAGS   = -Wall -Wextra -Werror -pedantic -std=c99 -g
CFLAGS   = -Wall -Wextra -pedantic -std=c99 -O
LDFLAGS  = -L$(ROOT)/lib 
LDLIBS   = -lbcl -lm

BIN= analogy

.PHONY: all
all: $(BIN)

analogy : main.o transfert.o synthese.o filtre.o features.o pyramid.o conversion.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o:%.c
	$(CC) $(CFLAGS) -c  $< -I../bcl/include 


.PHONY: clean cleanall
clean:
	$(RM) *.o *~ *.ppm
cleanall: clean
	$(RM) $(BIN) 
