#COMPILER
CC = cc
#OPTIONS
CFLAGS =  -O3 

#DEFINITIONS
#GDEFS  +=  -DVERBOSE  # activate extra checks in verbose mode
#GDEFS  +=  -DDEBUG    #fix the seed of the random number generator
GDEFS  +=  -DDUMP     #print configuration in a ASCII file for each realization
#GDEFS   += -DBINARY #with DUMP active print in BINARY   

OS      = $(shell uname -s)
ARCH    = $(shell uname -m)
HNAME   = $(shell uname -n)

ifeq ($(OS),Darwin) 
GDEFS +=-DMACOSX
FLAG = MACOSX
endif

ifeq ($(OS),Linux) 
FLAG = LINUX
LDFLAGS = -lm
endif


#FILES
INCLUDES =
OBJS = evolution.o lookuptable.o moving.o statistics.o  rnd.o debug.o 
SRCS = evolution.c lookuptable.c moving.c statistics.c rnd.c debug.c voter2d.c
HDRS =  common.h defs.h gvar.h functionsdef.h


default: v2d post

os:
	@echo "running on $(FLAG)"


v2d: voter2d.o ${OBJS}
	${CC} ${GDEFS}   ${INCLUDES} -o $@ voter2d.o ${OBJS} ${CFLAGS} ${LDFLAGS}


.c.o: 
	${CC}  ${GDEFS}  ${INCLUDES} -c $< ${CFLAGS}

clean:
	rm *.o 
	@echo "Cleaned."


post:
	@echo "compiled completed"
	@echo "output is:  v2d"

