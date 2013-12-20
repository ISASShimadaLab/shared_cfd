#!/usr/bin/python
import os
import sys
from checkout_utility import *

def checkout_flow():
	if(os.path.exists("checkout")):
	#	print "There is already checkout directory."
	#	print "Do you want to Overwrite? (Y/N)"
	#	line = sys.stdin.readline()
	#	if( line != "y\n" and line != "Y\n"):
	#		sys.exit(0)
		if     os.path.exists("checkout/condition.f90"):os.system("cp checkout/condition.f90 condition.bak.f90")
		if     os.path.exists("checkout/control.inp"):	os.system("cp checkout/control.inp   control.bak.inp")
		if not os.path.exists("grid_bak"):		os.system("mkdir grid_bak")
		os.system("cp checkout/grid/* grid_bak/")
		os.system("rm -r checkout")
	
	os.system("mkdir checkout")
	os.system("mkdir checkout/grid")
	os.system("mkdir checkout/result")
	
	#initialize arr_engage
	arr_engage = [
		"checkout/control.part.inp",\
		"checkout/Makefile.main",\
		"checkout/Makefile.var",\
		"checkout/main.part_top.f90",\
		"checkout/main.part_variable.f90",\
		"checkout/main.part_init.f90",\
		"checkout/main.part_point_implicit.f90",\
		"checkout/main.part_primitive.f90",\
		"checkout/main.part_main.f90"]
	for i,v in enumerate(arr_engage): arr_engage[i]=[v,""]
	
	
	####################################################################################
	fp        = open("checkout.inp","r")
	
	### CORE PARTS ###
	engage("core","checkout",arr_engage)
	
	# read grid file
	filename = read_control_next(fp)
	fgrid = open(filename,'r')
	Nplane = int(fgrid.readline().strip().split().pop())
	nijk=[]
	nijk_str=["","",""]
	nijk_max=[0,0,0]
	for i in range(Nplane):
		nijk_mono = fgrid.readline().strip().split()
		for i in range(0,3):
			nijk_mono[i] = int(nijk_mono[i])-1
			nijk_str[i]+=','+str(nijk_mono[i])
			nijk_max[i] = max(nijk_max[i],nijk_mono[i])
		if(nijk_mono[2] != 0):
			print "This program cannot be applied yet to 3D problem."
			print nijk_mono[2]
			sys.exit(1)
		nijk_mono.pop()
		nijk.append(nijk_mono)
	fgrid.close()
	os.system("cp "+filename+" checkout/grid/")
	
	####################################################################################
	print "***architecture***"
	val = read_control_next_int(fp)
	if(val == 0):
		print "\tPersonal Computer is selected."
		engage("architecture/PC","checkout",arr_engage)
	elif(val >= 1):
		print "\tM-System[Super Computer] is selected."
		engage("architecture/SC","checkout",arr_engage)
	
		import checkout.generate_proc_grid
		[MPINumGrid,bwmax]=checkout.generate_proc_grid.split_grid_from_Nproc(val,nijk)
		Nproc=0
		MPINumGridLine=['','']
		MPINumGridMax=[0,0]
		for MPINumGridMono in MPINumGrid:
			Nproc += MPINumGridMono[0]*MPINumGridMono[1]
			for i in range(2):
				MPINumGridMax[i]   = max(MPINumGridMax[i],MPINumGridMono[i])
				MPINumGridLine[i] += ','+str(MPINumGridMono[i])
		import shutil
		shutil.move("grid_separation.inp","checkout/")
		raw2pro("checkout/mpirun.raw.sh","checkout/mpirun.sh",[["Nproc",str(Nproc)]])
		fromto = [\
			 ["NGXMAX",str(MPINumGridMax[0])],\
			 ["NGYMAX",str(MPINumGridMax[1])],\
			 ["NGX",MPINumGridLine[0][1:]],\
			 ["NGY",MPINumGridLine[1][1:]],\
			 ["BWMAX",str(bwmax)]]
		raw2pro("checkout/mod_mpi.raw.f90","checkout/mod_mpi.f90",fromto)
		print "\tActual NumProc is set to"
		for i,MPINumGridMono in enumerate(MPINumGrid):
			print "\t\tPlane%3i : %3i x %3i = %3i" % (i+1,MPINumGridMono[0],MPINumGridMono[1],MPINumGridMono[0]*MPINumGridMono[1])
		print "\t\t\t\tSum = %3i" % (Nproc)
	else:
		print "\tOdd Input at architecture! value is ",val
		sys.exit(1)
	####################################################################################
	print "***dimension***"
	val = read_control_next_int(fp)
	if(val == 0):
		print "\ttwo dimension is selected."
		engage("dim/2d","checkout",arr_engage)
		DIMENSION = "2d"
	elif(val == 1):
		print "\taxi-symmetry is selected."
		engage("dim/q2d","checkout",arr_engage)
		DIMENSION = "q2d"
	else:
		print "\tOdd Input at dimension! value is ",val
		sys.exit(1)
	
	####################################################################################
	print "***time scale***"
	val_dt = read_control_next_int(fp)
	if(val_dt == 0):
		print "\tGlobal Time Step is selected."
		engage("time_step/global","checkout",arr_engage)
		DT_LOCAL_GLOBAL = "dt_grbl"
	elif(val_dt == 1):
		print "\tLocal Time Step is selected."
		engage("time_step/local","checkout",arr_engage)
		DT_LOCAL_GLOBAL = "dt_mat(i,j,plane)"
	else:
		print "\tOdd Input at time step! value is ",val_dt
		sys.exit(1)
	
	####################################################################################
	print "***time scheme***"
	val = read_control_next_int(fp)
	if(val == 0):
		print "\tEuler is selected."
		engage("time_scheme/euler","checkout",arr_engage)
		raw2pro("checkout/time_euler_viscous.raw.f90","checkout/time_euler_viscous.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
	elif(val == 1):
		print "\ttwo-step Runge-Kutta is selected."
		engage("time_scheme/RK2","checkout",arr_engage)
		raw2pro("checkout/time_RK2_viscous.raw.f90","checkout/time_RK2_viscous.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
	elif(val == 2):
		print "\tLU-SGS is selected."
		engage("time_scheme/LU-SGS/"+DIMENSION,"checkout",arr_engage)
		raw2pro("checkout/sch_lusgs.raw.f90","checkout/sch_lusgs.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
	elif(val == 3):
		print "\tNR is selected."
		engage("time_scheme/NR/"+DIMENSION,"checkout",arr_engage)
		raw2pro("checkout/sch_NR.raw.f90","checkout/sch_NR.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
	elif(val == 4):
		print "\tPre-conditioner is selected."
		engage("time_scheme/precon/"+DIMENSION,"checkout",arr_engage)
		raw2pro("checkout/sch_precon.raw.f90","checkout/sch_precon.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
		if(val_dt == 0):
			print "ERROR: This scheme cannot be used at global time step."
			sys.exit(1)
	elif(val == 5):
		if(DIMENSION == 'q2d'):
			print "ERROR: This scheme cannot be used at axi-symmetric problem now."
			sys.exit(1)
		if(val_dt == 1):
			print "ERROR: This scheme cannot be used at local time step."
			sys.exit(1)
		print "\tDual Time is selected."
		engage("time_scheme/dual/"+DIMENSION,"checkout",arr_engage)
		raw2pro("checkout/sch_dual.raw.f90","checkout/sch_dual.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
	elif(val == 6):
		print "\tPre-conditioner using LU-SGS is selected."
		engage("time_scheme/preconLU-SGS/"+DIMENSION,"checkout",arr_engage)
		raw2pro("checkout/sch_precon.raw.f90","checkout/sch_precon.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
		if(val_dt == 0):
			print "ERROR: This scheme cannot be used at global time step."
			sys.exit(1)
	else:
		print "\tOdd Input at time scheme! value is ",val
		sys.exit(1)
	
	####################################################################################
	print "***high order scheme***"
	val = read_control_next_int(fp)
	if(val == 0):
		print "\tUp Wind is selected."
		engage("high_order/upwind","checkout",arr_engage)
	elif(val == 1):
		print "\tMUSCL is selected."
		engage("high_order/muscl","checkout",arr_engage)
	else:
		print "\tOdd Input at high order scheme! value is ",val
		sys.exit(1)
	
	####################################################################################
	print "***viscosity***"
	val_vis = read_control_next_int(fp)
	if(val_vis == 0):
		print "\tnon-viscous flow is selected."
	elif(val_vis == 1):
		print "\tviscous flow is selected."
	else:
		print "\tOdd Input at viscosity! value is ",val_vis
		sys.exit(1)
	
	####################################################################################
	print "***thermal model***"
	val = read_control_next_int(fp)
	if(val == 0):
		print "\tideal gas is selected."
		engage("therm_lib/ideal","checkout",arr_engage)
		ABOUTNV="no-nV"
		nV  = 1
		nY  = 1
	
		# process thermal_model.f90
		fromto = [ \
			["KAPPA" ,"1.4d0"],\
			["RGAS"  ,"287d0"],\
			["NU"    ,"1.6d-5"]]
		raw2pro("checkout/thermal_model.raw.f90","checkout/thermal_model.f90",fromto)
	elif(val == 1 or val == 2):
		engage("therm_lib/chemkin/core","checkout",arr_engage)
		if(val == 1) :
			print "\tCold Flow using chemical kinetics database is selected."
			engage("therm_lib/chemkin/flow/cold","checkout",arr_engage)
		else:
			print "\tChemical Kinetics Model is selected."
			if(val_dt == 1):
				print "Error : Chemical Kinetics Model cannot be used at local time step."
				sys.exit(1)
			engage("therm_lib/chemkin/flow/reactive","checkout",arr_engage)

		# process mod_chem.f90
		if not os.path.exists("chem.inp"):
			print "ERROR:Can't find 'chem.inp'. Please try again."
			sys.exit(1)
		os.system("cp chem.inp checkout/")
		import store.checkout.ckinterp
		val = store.checkout.ckinterp.ckinterp()
		val = map(str,val)
		val[2]=max(val[2],1)
		fromto = [ \
			["NE" ,val[0]],\
			["NS" ,val[1]],\
			["NR" ,val[2]]]
		raw2pro("checkout/mod_chem.raw.f90","checkout/mod_chem.f90",fromto)
	
		nY  = val[1]
		nV  = val[1]
		ABOUTNV="no-nV"
	elif(val == 3):
		print "\tflame sheet model is selected."
		engage("therm_lib/NASA/core","checkout",arr_engage)
		engage("therm_lib/NASA/flame_sheet/flow","checkout",arr_engage)
	
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
		raw2pro("checkout/mod_chem.raw.f90","checkout/mod_chem.f90",fromto)
		nY  = 2
		nV  = 3
		ABOUTNV="with-nV"
	elif(val == 4):
		print "\tPerfect Equilibrium model is selected."
		engage("therm_lib/NASA/core","checkout",arr_engage)
		engage("therm_lib/NASA/cea/flow","checkout",arr_engage)
	
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
		raw2pro("checkout/mod_chem.raw.f90","checkout/mod_chem.f90",fromto)
		nY  = 2
		nV  = val[1]
		ABOUTNV="with-nV"
	else:
		print "\tOdd Input at thermal model! value is ",val
		sys.exit(1)
	####################################################################################
	# close checkout.inp
	fp.close()
	
	#viscosity selection
	if(val_vis == 0):
		engage("viscosity/non-viscous","checkout",arr_engage)
	elif(val_vis == 1):
		engage("viscosity/viscous/"+ABOUTNV+"/"+DIMENSION,"checkout",arr_engage)
	
	# generation of condition.raw.f90
	print "***generation of condition.raw.f90***"
	engage("cond/core","checkout",arr_engage)
	engage("cond/"+ABOUTNV,"checkout",arr_engage)
	import checkout.gen_cond
	checkout.gen_cond.gen_cond(filename,ABOUTNV)
	print "\tDone."

	#expand variables
	part_control         = arr_engage[0][1]
	part_Makefile_main   = arr_engage[1][1]
	part_Makefile_var    = arr_engage[2][1]
	part_top             = arr_engage[3][1]
	part_variable        = arr_engage[4][1]
	part_init            = arr_engage[5][1]
	part_point_implicit  = arr_engage[6][1]
	part_primitive       = arr_engage[7][1]
	part_main            = arr_engage[8][1]
	
	# generate control.raw.inp
	fcontrol  = open("checkout/control.raw.inp","w")
	fcontrol.write(part_control)
	fcontrol.close()
	
	# generate n_grid.f90
	fromto = [\
		["NPLANE",str(Nplane)],\
		["NumI",nijk_str[0][1:]],\
		["NumJ",nijk_str[1][1:]],\
		["NIMAX",str(nijk_max[0])],\
		["NJMAX",str(nijk_max[1])],\
		["NumY",str(nY)],\
		["NV",  str(nV)],\
		["GridFileName",os.path.basename(filename)]]
	raw2pro("checkout/n_grid.raw.f90","checkout/n_grid.f90",fromto)
	
	# generate Makefile
	fMakefile = open("checkout/Makefile","w")
	append_file_to_file(fMakefile,"checkout/Makefile.head")
	fMakefile.write(part_Makefile_var)
	fMakefile.write(part_Makefile_main)
	append_file_to_file(fMakefile,"checkout/Makefile.foot")
	fMakefile.close()
	
	# generate main.f90
	fmain = open("checkout/main.f90","w")
	fmain.write(part_top)
	append_file_to_file(fmain,"checkout/main.head.f90")
	fmain.write(part_variable)
	fmain.write(part_init)
	fromto = [\
	["part_top\n"           , part_top],\
	["part_init\n"          , part_init],\
	["part_point_implicit\n", part_point_implicit],\
	["part_primitive\n"     , part_primitive],\
	["part_main\n"          , part_main]]
	raw2pro("checkout/main.body.raw.f90","checkout/main.body.f90",fromto)
	append_file_to_file(fmain,"checkout/main.body.f90")
	if(part_point_implicit == ""):
		raw2pro("checkout/main.main.raw.f90","checkout/main.main.f90",fromto)
		if(os.path.exists("checkout/main.main.point_implicit.f90")):
			os.system("rm checkout/main.main.point_implicit.f90")
	else:
		if(not os.path.exists("checkout/main.main.point_implicit.f90")):
			print "This time scheme cannnot be applied to point implicit scheme."
		raw2pro("checkout/main.main.point_implicit.f90","checkout/main.main.f90",fromto)
		os.system("rm checkout/main.main.raw.f90")
	append_file_to_file(fmain,"checkout/main.main.f90")
	append_file_to_file(fmain,"checkout/main.foot.f90")
	fmain.close()
	
	#remove unnecessary files
	if os.path.exists("checkout/__init__.py"): os.remove("checkout/__init__.py")
	os.system("rm -f checkout/*.pyc")
