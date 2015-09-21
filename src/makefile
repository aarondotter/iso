#for MESA pieces
include $(MESA_DIR)/utils/makefile_header

INCLUDE=-I$(MESASDK_ROOT)/include -I$(MESA_DIR)/include
LOAD_MESA=-L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)

FC=gfortran 
DEBUG=-g -fbacktrace -fcheck=all 
WARN=-Wall -Wextra -pedantic
OMP=-fopenmp
OPT=-O3
FFLAGS= -ffree-form $(OPT) $(WARN) $(DEBUG)

ISO=make_iso
EEP=make_eep
BLEND=blend_eeps

all: $(EEP) $(ISO) $(BLEND)

$(ISO): make_iso.o iso_eep_support.o iso_color.o BC_tables.o
	$(FC) $(OMP) -o ../$(ISO) $^ $(LOAD_MESA)

$(EEP): make_eep.o phases.o eep.o iso_eep_support.o 
	$(FC) $(OMP) -o ../$(EEP) $^ $(LOAD_MESA)

$(BLEND): blend_eeps.o iso_eep_support.o
	$(FC) $(OMP) -o ../$(BLEND) $^ $(LOAD_MESA)

make_iso.o: make_iso.f90 iso_eep_support.o iso_color.o BC_tables.o

make_eep.o: make_eep.f90 phases.o eep.o iso_eep_support.o

blend_eeps.o: blend_eeps.f90 iso_eep_support.o

eep.o: eep.f90 iso_eep_support.o

phases.o: phases.f90 iso_eep_support.o

iso_eep_support.o: iso_eep_support.f90

iso_color.o: iso_eep_support.o BC_tables.o

BC_tables.o: BC_tables.f90

%.o:%.f90
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

check: $(ISO) $(EEP)
	cd ..; ./make_iso input.041214

clean:
	rm -f *.o *.mod *~ ../$(EEP) ../$(ISO) ../$(BLEND)
