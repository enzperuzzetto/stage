ROOT=../../bcl

CPPFLAGS = -I$(ROOT)/include
#CFLAGS   = -Wall -Wextra -Werror -pedantic -std=c99 -g
CFLAGS   = -Wall -Wextra  -pedantic -std=c99 -O
LDFLAGS  = -L$(ROOT)/lib 
LDLIBS   = -lbcl -lm

BIN=\
	yiq\
	yiq2\
	pyramid_gaussian\
	patch_match\
	random\
	cannaux\
	option\
	texture\
	texture_pyramid

.PHONY: all
all: $(OBJ) $(BIN)


.PHONY: clean cleanall
clean:
	$(RM) *.o *~ *.ppm
cleanall: clean
	$(RM) $(BIN) 
