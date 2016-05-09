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

CMD=make_cmd
ISO=make_iso
EEP=make_eep
TRK=make_track
BLEND=blend_eeps
INTMET=iso_interp_met
INTAGE=iso_interp_age
INTEEP=eep_interp_met
CHA=charlies_angel

basic: $(EEP) $(ISO) 
all: $(EEP) $(ISO) $(BLEND) $(TRK) $(INTMET) $(INTAGE) $(CMD) $(CHA) $(INTEEP)

$(CMD): make_cmd.o iso_eep_color.o iso_eep_support.o BC_tables.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(CMD) $^ $(LOAD_MESA)

$(INTEEP): eep_interp_met.o iso_eep_support.o iso_eep_color.o BC_tables.o iso_interp_support.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(INTEEP) $^ $(LOAD_MESA)

$(INTMET): iso_interp_met.o iso_eep_support.o iso_interp_support.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(INTMET) $^ $(LOAD_MESA)

$(INTAGE): iso_interp_age.o iso_eep_support.o iso_interp_support.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(INTAGE) $^ $(LOAD_MESA)

$(ISO): make_iso.o iso_eep_support.o 
	$(FC) $(OMP) -o $(ROOT_DIR)/$(ISO) $^ $(LOAD_MESA)

$(EEP): make_eep.o phases.o eep.o iso_eep_support.o 
	$(FC) $(OMP) -o $(ROOT_DIR)/$(EEP) $^ $(LOAD_MESA)

$(BLEND): blend_eeps.o iso_eep_support.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(BLEND) $^ $(LOAD_MESA)

$(TRK): make_track.o iso_eep_support.o iso_eep_color.o BC_tables.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(TRK) $^ $(LOAD_MESA)

$(CHA): charlies_angel.o iso_eep_support.o
	$(FC) $(OMP) -o $(ROOT_DIR)/$(CHA) $^ $(LOAD_MESA)

eep_interp_met.o: $(SRC_DIR)/eep_interp_met.f90 iso_eep_support.o iso_interp_support.o iso_eep_color.o BC_tables.o

charlies_angel.o: $(SRC_DIR)/charlies_angel.f90 iso_eep_support.o

iso_interp_support.o: $(SRC_DIR)/iso_interp_support.f90 iso_eep_support.o

iso_interp_met.o: $(SRC_DIR)/iso_interp_met.f90 iso_eep_support.o iso_interp_support.o

iso_interp_age.o: $(SRC_DIR)/iso_interp_age.f90 iso_eep_support.o iso_interp_support.o

make_track.o: $(SRC_DIR)/make_track.f90 iso_eep_support.o iso_eep_color.o

make_cmd.o: $(SRC_DIR)/make_cmd.f90 iso_eep_support.o iso_eep_color.o BC_tables.o

make_iso.o: $(SRC_DIR)/make_iso.f90 iso_eep_support.o

make_eep.o: $(SRC_DIR)/make_eep.f90 phases.o eep.o iso_eep_support.o

blend_eeps.o: $(SRC_DIR)/blend_eeps.f90 iso_eep_support.o

eep.o: $(SRC_DIR)/eep.f90 iso_eep_support.o

phases.o: $(SRC_DIR)/phases.f90 iso_eep_support.o

iso_eep_support.o: $(SRC_DIR)/iso_eep_support.f90

iso_eep_color.o: iso_eep_support.o BC_tables.o

BC_tables.o: $(SRC_DIR)/BC_tables.f90

%.o : $(SRC_DIR)/%.f90
	$(FC) $(OMP) $(FFLAGS) $(INCLUDE) -c $<

clean:
	rm -f $(ROOT_DIR)/$(EEP) $(ROOT_DIR)/$(ISO) $(ROOT_DIR)/$(BLEND) $(ROOT_DIR)/$(TRK) $(ROOT_DIR)/$(INTMET) $(ROOT_DIR)/$(INTAGE) $(ROOT_DIR)/$(CMD)
	rm -f *.mod *.o
