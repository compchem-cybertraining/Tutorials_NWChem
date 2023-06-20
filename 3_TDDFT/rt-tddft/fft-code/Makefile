#
## Makefile
#
EXE     = do_fft.x
SRCS = \
do_fft.f90
OBJS =$(SRCS:.f=.o)
FC	 = gfortran
${EXE}: $(OBJS) ;
	${FC} -o ${EXE} ${OBJS} ./dfftpack/libdfftpack.a -lm
