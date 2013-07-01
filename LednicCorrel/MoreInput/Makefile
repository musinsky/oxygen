NAME = lednic_input
SRCS = $(NAME).F
OBJS = $(NAME).o
OUT = $(NAME).exe

FCC = g77
CERNFLAGS = -fno-automatic -fno-second-underscore -fugly-complex \
 -DCERNLIB_LINUX -DCERNLIB_UNIX -DCERNLIB_LNX -DCERNLIB_QMGLIBC

CERNLIBS = `cernlib mathlib`

all: ${OUT}

${OUT}: ${OBJS}
	${FCC} -o $@ ${OBJS} ${CERNLIBS}

${OBJS}: ${SRCS}
	${FCC} -c ${CERNFLAGS} ${SRCS}
clean:
	rm -f ${NAME}.o
