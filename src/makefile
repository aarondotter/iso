#for MESA pieces
include $(MESA_DIR)/utils/makefile_header

INCLUDE=-I$(MESASDK_ROOT)/include -I$(MESA_DIR)/include
LOAD_MESA=-L$(MESA_LIB_DIR) -linterp_1d -lnum -lmtx -lutils $(LOAD_BLAS) $(LOAD_LAPACK)

FC=gfortran 
DEBUG=-g -fbacktrace -fcheck=all 
WARN=-Wall -Wextra -pedantic
OMP=-fopenmp
OPT=-O0
FFLAGS= -ffree-form $(OPT) $(WARN) $(DEBUG)

ISO=make_iso
EEP=make_eeps
BLEND=blend_eeps

all: $(EEP) $(ISO) $(BLEND)

$(ISO): make_iso.o iso_eep_support.o 
	$(FC) $(OMP) -o ../$(ISO) $^ $(LOAD_MESA)

$(EEP): make_eeps.o phases.o eep.o iso_eep_support.o 
	$(FC) $(OMP) -o ../$(EEP) $^ $(LOAD_MESA)

$(BLEND): blend_eeps.o iso_eep_support.o
	$(FC) $(OMP) $(FFLAGS) -o ../$(BLEND) $^ $(LOAD_MESA)

make_iso.o: make_iso.f iso_eep_support.o
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

make_eeps.o: make_eeps.f phases.o eep.o iso_eep_support.o
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

blend_eeps.o: blend_eeps.f iso_eep_support.o
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

eep.o: eep.f iso_eep_support.o
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

phases.o: phases.f iso_eep_support.o
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

iso_eep_support.o: iso_eep_support.f
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

check: $(ISO) $(EEP)
	cd ..; ./make_iso input.041214

clean:
	rm -f *.o *.mod *~ ../$(EEP) ../$(ISO) ../$(BLEND)
