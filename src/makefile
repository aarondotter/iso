#for MESA pieces
include $(MESA_DIR)/utils/makefile_header

COLOR_DIR=/home/dotter/science/colors
LOAD_COLOR=-L$(COLOR_DIR)/make -lcolor
INCLUDE=-I$(MESASDK_ROOT)/include -I$(MESA_DIR)/include -I$(COLOR_DIR)/make
LOAD_MESA=-L$(MESA_LIB_DIR) -linterp_1d -lnum -lmtx -lutils

FC=gfortran 
DEBUG=-g -fbacktrace -fcheck=all 
WARN=-Wall -Wextra -pedantic
OMP=-fopenmp
OPT=-O3
FFLAGS= -ffree-form $(OPT) $(WARN) $(DEBUG)

ISO=make_iso
EEP=make_eeps
BLEND=blend_eeps

all: $(EEP) $(ISO) $(BLEND)

$(ISO): make_iso.o iso_eep_support.o iso_color.o
	$(FC) $(OMP) -o ../$(ISO) $^ $(LOAD_MESA) $(LOAD_COLOR)

$(EEP): make_eeps.o phases.o eep.o iso_eep_support.o 
	$(FC) $(OMP) -o ../$(EEP) $^ $(LOAD_MESA)

$(BLEND): blend_eeps.o iso_eep_support.o
	$(FC) $(OMP) -o ../$(BLEND) $^ $(LOAD_MESA)

make_iso.o: make_iso.f90 iso_eep_support.o iso_color.o

make_eeps.o: make_eeps.f90 phases.o eep.o iso_eep_support.o

blend_eeps.o: blend_eeps.f90 iso_eep_support.o

eep.o: eep.f90 iso_eep_support.o

phases.o: phases.f90 iso_eep_support.o

iso_eep_support.o: iso_eep_support.f90

%.o:%.f90
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

check: $(ISO) $(EEP)
	cd ..; ./make_iso input.041214

clean:
	rm -f *.o *.mod *~ ../$(EEP) ../$(ISO) ../$(BLEND)
