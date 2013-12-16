##### SELECT ARCHITECTURE#######################
#PERSONAL COMPUTER
arc=mod_mpi_dummy.o
#compiler=ifort -convert big_endian -openmp -fast -parallel
#compiler=ifort -convert big_endian -fast -parallel
compiler=ifort -convert big_endian -check all -traceback -g
##### END SELECT ARCHITECTURE###################

core=n_grid.o
thermalmodel=sub_chem.o mod_chem.o func_chem.o LU.o
all: driver

#dummy
n_grid.o: n_grid.f90
	$(compiler) -c $<
mod_mpi_dummy.o: mod_mpi_dummy.f90 n_grid.o
	$(compiler) -c $<

#chemlib
LU.o: LU.f90
	$(compiler) -c $<
mod_chem.o:       mod_chem.f90
	$(compiler) -c $<
func_chem.o:      func_chem.f90 mod_chem.o
	$(compiler) -c $<
sub_chem.o:       sub_chem.f90 func_chem.o LU.o
	$(compiler) -c $<
thermal_model.o:  thermal_model.f90 sub_chem.o
	$(compiler) -c $<
driver : driver.f90 $(core) $(arc) $(thermalmodel)
	$(compiler) -o $@ $^

clean:
	rm -rf driver *.mod *.o *.dat *.eps *.vtk *.bin temp.* result/*.*

