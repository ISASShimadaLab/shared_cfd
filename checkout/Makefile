##### SELECT ARCHITECTURE#######################
#PERSONAL COMPUTER
arc=mod_mpi_dummy.o
#compiler=ifort -convert big_endian -openmp -fast -parallel
#compiler=ifort -convert big_endian -fast -parallel
compiler=ifort -convert big_endian -check all -traceback -g

##### END SELECT ARCHITECTURE###################
program=driver
all: $(program)

#dummy
n_grid_dummy.o: n_grid_dummy.f90
	$(compiler) -c $<
mod_mpi_dummy.o: mod_mpi_dummy.f90 n_grid_dummy.o
	$(compiler) -c $<

#chemlib
LU.o: LU.f90
	$(compiler) -c $<
mod_chem.o:  mod_chem.f90
	$(compiler) -c $<
func_chem.o:  func_chem.f90 mod_chem.o
	$(compiler) -c $<
sub_chem.o:  sub_chem.f90  func_chem.o LU.o mod_chem.o
	$(compiler) -c $<
chemlib=sub_chem.o func_chem.o mod_chem.o LU.o

driver : driver.f90 n_grid_dummy.o mod_mpi_dummy.o $(chemlib)
	$(compiler) -o $@ $^

clean:
	rm -rf $(program) *.mod *.o *.dat *.eps *.vtk *.bin temp.* result/*.* 
tar:
	tar cvf  result.tar.gz result/*.bin geometry.bin RMS.dat
