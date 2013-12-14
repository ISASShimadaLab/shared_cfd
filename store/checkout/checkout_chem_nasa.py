#!/usr/bin/python
import os
import sys
from checkout_utility import *

def checkout_chem_nasa():
	if(os.path.exists("checkout_chem")):
		os.system("rm -r checkout_chem")
	
	os.system("mkdir checkout_chem")
	
	#initialize arr_engage
	arr_engage = [
		"checkout_chem/control.part.inp",\
		"checkout_chem/Makefile.main",\
		"checkout_chem/Makefile.var"]
	for i,v in enumerate(arr_engage): arr_engage[i]=[v,""]
	
	####################################################################################
	fp        = open("checkout_chem.inp","r")
	
	### CORE PARTS ###
	engage("chem_drive","checkout_chem",arr_engage)
	os.system("cp store/core/n_grid.raw.f90 checkout_chem/")
	
	####################################################################################
	print "***architecture***"
	print "\tPersonal Computer is selected."
	engage("architecture/PC","checkout_chem",arr_engage)

	####################################################################################
	print "***thermal model***"
	val = read_control_next_int(fp)
	if(val == 0):
		print "\tflame sheet model is selected."
		engage("therm_lib/NASA/core","checkout_chem",arr_engage)
		engage("therm_lib/NASA/driver","checkout_chem",arr_engage)
		
		# process mod_chem.f90
		if not os.path.exists("chem.inp"):
			print "ERROR:Can't find 'chem.inp'. Please try again."
			sys.exit(1)
		os.system("cp chem.inp checkout/")
		import store.checkout.ckinterp
		val = store.checkout.ckinterp.ckinterp()
		val = map(str,val)
		fromto = [ \
			["NE",val[0]],\
			["NS",val[1]]]
		raw2pro("checkout_chem/mod_chem.raw.f90","checkout_chem/mod_chem.f90",fromto)
		nY  = 2
		nV  = 3
		ABOUTNV="with-nV"
	else:
		print "\tOdd Input at thermal model! value is ",val
		sys.exit(1)
	####################################################################################
	# close checkout.inp
	fp.close()

	#expand variables
	part_control         = arr_engage[0][1]
	part_Makefile_main   = arr_engage[1][1]
	part_Makefile_var    = arr_engage[2][1]
	
	# generate control.raw.inp
	fcontrol  = open("checkout_chem/control.raw.inp","w")
	fcontrol.write(part_control)
	fcontrol.close()
	
	# generate n_grid.f90
	fromto = [\
		["NPLANE","1"],\
		["NumI","1"],\
		["NumJ","1"],\
		["NIMAX","1"],\
		["NJMAX","1"],\
		["NumY",str(nY)],\
		["NV",  str(nV)],\
		["GridFileName","dummy"]]
	raw2pro("checkout_chem/n_grid.raw.f90","checkout_chem/n_grid.f90",fromto)
	
	# generate Makefile
	fMakefile = open("checkout_chem/Makefile","w")
	append_file_to_file(fMakefile,"checkout_chem/Makefile.head")
	fMakefile.write(part_Makefile_var)
	fMakefile.write(part_Makefile_main)
	append_file_to_file(fMakefile,"checkout_chem/Makefile.foot")
	fMakefile.close()
	
	#remove unnecessary files
	if os.path.exists("checkout_chem/__init__.py"): os.remove("checkout_chem/__init__.py")
	os.system("rm -f checkout_chem/*.pyc")
