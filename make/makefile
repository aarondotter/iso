#for MESA pieces, you need to have MESA_DIR environment variable set
include $(MESA_DIR)/utils/makefile_header

INCLUDE=-I$(MESASDK_ROOT)/include -I$(MESA_DIR)/include
LOAD_MESA=-L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)

ROOT_DIR=/home/dotter/science/iso
SRC_DIR=$(ROOT_DIR)/src

FC=gfortran 
DEBUG=-g -fbacktrace -fcheck=all 
WARN=-Wall -Wextra -pedantic
OMP=-fopenmp
OPT=-O3
FFLAGS= -ffree-form $(OPT) $(WARN) $(DEBUG)

ISO=make_iso
EEP=make_eep
TRK=make_track
BLEND=blend_eeps
INTERP=iso_interp

all: $(EEP) $(ISO) $(BLEND) $(TRK) $(INTERP)

$(INTERP): iso_interp.o iso_eep_support.o iso_color.o BC_tables.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(INTERP) $^ $(LOAD_MESA)

$(ISO): make_iso.o iso_eep_support.o iso_color.o BC_tables.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(ISO) $^ $(LOAD_MESA)

$(EEP): make_eep.o phases.o eep.o iso_eep_support.o 
	$(FC) $(OMP) -o $(ROOT_DIR)/$(EEP) $^ $(LOAD_MESA)

$(BLEND): blend_eeps.o iso_eep_support.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(BLEND) $^ $(LOAD_MESA)

$(TRK): make_track.o iso_eep_support.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(TRK) $^ $(LOAD_MESA)

iso_interp.o: $(SRC_DIR)/iso_interp.f90 iso_eep_support.o

make_track.o: $(SRC_DIR)/make_track.f90 iso_eep_support.o

make_iso.o: $(SRC_DIR)/make_iso.f90 iso_eep_support.o iso_color.o BC_tables.o

make_eep.o: $(SRC_DIR)/make_eep.f90 phases.o eep.o iso_eep_support.o

blend_eeps.o: $(SRC_DIR)/blend_eeps.f90 iso_eep_support.o

eep.o: $(SRC_DIR)/eep.f90 iso_eep_support.o

phases.o: $(SRC_DIR)/phases.f90 iso_eep_support.o

iso_eep_support.o: $(SRC_DIR)/iso_eep_support.f90

iso_color.o: iso_eep_support.o BC_tables.o

BC_tables.o: $(SRC_DIR)/BC_tables.f90

%.o : $(SRC_DIR)/%.f90
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

check: $(ISO) $(EEP)
	cd ..; ./make_iso input.041214

clean:
	rm -f $(ROOT_DIR)/$(EEP) $(ROOT_DIR)/$(ISO) $(ROOT_DIR)/$(BLEND) $(ROOT_DIR)/$(TRK) $(ROOT_DIR)/$(INTERP) *.mod *.o
