ROOT=../../bcl

CPPFLAGS = -I$(ROOT)/include
CFLAGS   = -Wall -Wextra -Werror -pedantic -std=c99
LDFLAGS  = -L$(ROOT)/lib 
LDLIBS   = -lbcl -lm

BIN=\
	yiq\
	pyramid_gaussian

.PHONY: all
all: $(OBJ) $(BIN)


.PHONY: clean cleanall
clean:
	$(RM) *.o *~ *.ppm
cleanall: clean
	$(RM) $(BIN) 
