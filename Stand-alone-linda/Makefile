# Makefile

SHELL := /bin/bash
FORT?=gfortran
FORTFLAGS?=-fdefault-real-8 -O0 -g -ffpe-trap=zero,invalid,overflow,underflow -fimplicit-none -Wunused -Wall
EXE = soil_model
INC = /usr/local/netcdf-4.3.2-4.2-gfortran/include/
LIB = /usr/local/netcdf-4.3.2-4.2-gfortran/lib64/ -lnetcdf -lnetcdff
#compilation
.PHONY: all info
all: $(EXE)

#main executable
$(EXE): soil_model.o
	$(FORT) $(FORTFLAGS) data_parameters.o data_soil.o data_namelist.o soil_model.o -L$(LIB) $(OBJ) -o $(EXE)
# modules
data_parameters.o: data_parameters.f90
	$(FORT) $(FORTFLAGS) -c data_parameters.f90

data_soil.o: data_soil.f90 data_parameters.o
	$(FORT) $(FORTFLAGS) -c data_soil.f90

data_namelist.o: data_namelist.f90 data_parameters.o
	$(FORT) $(FORTFLAGS) -c data_namelist.f90
# main program
soil_model.o: soil_model.f90 data_namelist.o data_soil.o
	$(FORT) $(FORTFLAGS) -c soil_model.f90 -I$(INC)
info :
	@echo "generating compile information"
	@-rm -f .fconfig
	@echo "Exe              : $(EXE)" >> .fconfig
	@echo "Compiler command : $(FORT)" >> .fconfig
	@echo "Compiler includes: $(INC)" >> .fconfig
	@echo "Compiler flags   : $(FORTFLAGS)" >> .fconfig
	@echo "Linker command   : $(LD)" >> .fconfig
	@echo "Linker libraries : $(LIB)" >> .fconfig

.PHONY: clean
clean:
	-rm $(EXE)
	-rm *.o
	-rm *.mod
	-rm *~
del_temp:
	-rm *~
