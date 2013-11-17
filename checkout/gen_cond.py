#!/usr/bin/python
from mod_cond import *
from read_plot3d import *
from set_cond_var import *
from fetchMPI import *
from out_cond import *
import os

def gen_cond(filename):
	[Nplane,nijk,BD] = read_plot3d(filename)
	
	touch = set_touch(BD)
	toSetManually = set_toSetManually(Nplane,nijk,touch)
	
	out_cond(toSetManually)
	
	[bmpi,Nproc]=fetch_myid(Nplane,nijk)
	[touch,MPIcomm]= modify_touch(touch,bmpi,nijk,Nproc)
	
	out_cut(touch,nijk,Nproc)
	out_MPI(MPIcomm,nijk,Nproc)
	os.system('rm -f checkout/condition.head.f90 checkout/condition.tail.f90 checkout/fetchMPI.py checkout/mod_cond.py checkout/out_cond.py checkout/read_plot3d.py checkout/set_cond_var.py')

