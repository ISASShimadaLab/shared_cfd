##### SELECT ARCHITECTURE#######################
#PERSONAL COMPUTER
arc=mod_mpi_dummy.o
compiler=ifort -convert big_endian -openmp -fast -parallel
#compiler=ifort -convert big_endian -fast -parallel
#compiler=ifort -convert big_endian -check all -traceback -g
##### END SELECT ARCHITECTURE###################

core =  n_grid.o read_control.o scheme.o variable.o prmtr.o geometry.o set_dt.o check_convergence.o inout.o 
time_scheme= time_precon_viscous.o sch_precon.o var_precon.o
highorder = muscl.o
viscous = sch_viscous.o
thermalmodel=thermal_model.o

all: main

#Core
n_grid.o:n_grid.f90
	$(compiler) -c $<
prmtr.o: prmtr.f90
	$(compiler) -c $<

variable.o:variable.f90 n_grid.o
	$(compiler) -c $<
mod_mpi.o:mod_mpi.f90 n_grid.o
	$(compiler) -c $<
mod_mpi_dummy.o:mod_mpi_dummy.f90 n_grid.o
	$(compiler) -c $<

geometry.o: geometry.f90 variable.o $(arc)
	$(compiler) -c $<

scheme.o: scheme.f90     variable.o $(arc)
	$(compiler) -c $<
inout.o:  inout.f90                        variable.o $(arc)
	$(compiler) -c $<
read_control.o:read_control.f90            variable.o $(arc)
	$(compiler) -c $<
check_convergence.o: check_convergence.f90 variable.o $(arc)
	$(compiler) -c $<
set_dt.o: set_dt.f90                       variable.o $(arc)
	$(compiler) -c $<

#Time Scheme
var_precon.o: var_precon.f90 n_grid.o
	$(compiler) -c $<
sch_precon.o: sch_precon.f90 var_precon.o variable.o $(arc)
	$(compiler) -c $<
time_precon_viscous.o: time_precon_viscous.f90 var_precon.o
	$(compiler) -c $<


#High order scheme
muscl.o: muscl.f90 $(arc)
	$(compiler) -c $<


#Viscosity
sch_viscous.o: sch_viscous.f90 variable.o $(arc)
	$(compiler) -c $<

#thermal model
thermal_model.o:  thermal_model.f90 variable.o $(arc)
	$(compiler) -c $<


#user definition
condition.o: condition.f90 $(core) $(arc) $(time_scheme) $(viscous) $(thermalmodel) $(highorder)
	$(compiler) -c $<

main : main.f90 $(core) $(arc) $(time_scheme) $(viscous) $(thermalmodel) $(highorder) condition.o
	$(compiler) -o $@ $^

clean: remove
	rm -rf main *.mod *.o *.dat *.eps *.vtk *.bin temp.* result/*.*
remove:
	rm -rf *.o* *.e* result/* result.tar.gz

tar:
	tar cvf  result.tar result/*.bin geometry.bin RMS.dat
	gzip -f result.tar

