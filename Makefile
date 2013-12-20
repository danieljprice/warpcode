.KEEP_STATE:
FC = gfortran
FFLAGS = -Wall -Wextra -O3 -g -ffast-math -fopenmp
DEBUGFLAG = -fcheck=all -ffpe-trap=invalid -finit-real=nan -finit-integer=nan -fbacktrace

SOURCES= getalpha2_gordon.f warpeddisc.f90 main.f90
OBJ1=$(SOURCES:.f90=.o)
OBJ=$(OBJ1:.f=.o)

ifeq ($(DEBUG), yes)
    FFLAGS += ${DEBUGFLAG}
    FFLAGS := $(FFLAGS:-O3=-O0)
endif

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

default: disc

.PHONY: disc

disc: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	\rm -f *.o *.mod disc
