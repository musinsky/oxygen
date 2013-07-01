NAME = lednic_input
SRCS = $(NAME).F
OUT = $(NAME).exe

FCC = gfortran
CERNLIBS = `cernlib mathlib`

all:	${OUT}

${OUT}:	${SRCS}
	${FCC} -o $@ $^ ${CERNLIBS}

clean:
	rm -f ${OUT}
