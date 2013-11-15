compiler=ifort -fast -parallel -convert big_endian -assume byterecl

remote=/large/y/y535/2dcfd/share_code
directory=../checkout/

#program=postprocess change_num_step flow_history
#program=postprocess
program = postprocess
all: $(program)

#main
n_grid.o:$(directory)n_grid.f90
	$(compiler) -c $<
module.o:$(directory)module.f90  n_grid.o
	$(compiler) -c $<

postprocess.o:postprocess.f90 n_grid.o
	$(compiler) -c $<
postprocess  : postprocess.o n_grid.o
	$(compiler) -o $@ $^

change_num_step  : change_num_step.f90 n_grid.o module.o
	$(compiler) -o $@ $^

flow_history     : flow_history.f90 n_grid.o module.o
	$(compiler) -o $@ $^

modify           : modify.f90 n_grid.o module.o
	$(compiler) -o $@ $^

clean: 
	rm -rf $(program) *.o *.mod result/* result.tar.gz geometry.bin RMS.dat temp*

pull:
	scp maja:$(remote)/result.tar.gz .
	tar xzvf result.tar.gz
