#!/usr/bin/python
import os
import sys

#########################################################################################
############################## utilities ################################################
#########################################################################################

def append_file_to_file(fp,filename):
	ftmp = open(filename)
	fp.write(ftmp.read())
	ftmp.close()
	os.system("rm "+filename)

def file_to_str(filename):
	ftmp = open(filename)
	string = ftmp.read()
	ftmp.close()
	os.system("rm "+filename)
	return string

def raw2pro(name_raw,name_pro,fromto):
	fraw = open(name_raw,"r")
	fpro = open(name_pro,"w")
	for line in fraw:
		for element in fromto:
			line = line.replace(element[0],element[1])
		fpro.write(line)
	fraw.close()
	fpro.close()
	os.system("rm "+name_raw)

def read_control_next(fp):
	string = fp.readline()
	string = string[string.index(':')+1:string.index('#')-1].strip()
	return string

def read_control_next_int(fp):
	return int(read_control_next(fp))

def read_control_next_split(fp,num,place):
	val = read_control_next(fp).split()
	if(len(val) != num):
		print "Error: The number of parameters for "+place+" is odd."
		sys.exit(1)
	return val

def engage(place):
	global fcontrol, part_Makefile_main, part_Makefile_var,\
		part_top,part_variable,part_init,part_point_implicit,part_primitive,part_main
	os.system("cp store/"+place+"/* checkout/")

	filename = "checkout/control.part.inp"
	if(os.path.exists(filename)):
		append_file_to_file(fcontrol,filename)
	filename = "checkout/Makefile.main"
	if(os.path.exists(filename)):
		part_Makefile_main += file_to_str(filename)
	filename = "checkout/Makefile.var"
	if(os.path.exists(filename)):
		part_Makefile_var  += file_to_str(filename)
	filename = "checkout/main.part_top.f90"
	if(os.path.exists(filename)):
		part_top           += file_to_str(filename)
	filename = "checkout/main.part_variable.f90"
	if(os.path.exists(filename)):
		part_variable      += file_to_str(filename)
	filename = "checkout/main.part_init.f90"
	if(os.path.exists(filename)):
		part_init          += file_to_str(filename)
	filename = "checkout/main.part_point_implicit.f90"
	if(os.path.exists(filename)):
		part_point_implicit+= file_to_str(filename)
	filename = "checkout/main.part_primitive.f90"
	if(os.path.exists(filename)):
		part_primitive     += file_to_str(filename)
	filename = "checkout/main.part_main.f90"
	if(os.path.exists(filename)):
		part_main          += file_to_str(filename)


#########################################################################################
############################## start of main program ####################################
#########################################################################################

if(os.path.exists("checkout")):
#	print "There is already checkout directory."
#	print "Do you want to Overwrite? (Y/N)"
#	line = sys.stdin.readline()
#	if( line != "y\n" and line != "Y\n"):
#		sys.exit(0)
	os.system("rm -r checkout")

os.system("mkdir checkout")
os.system("mkdir checkout/grid")
os.system("mkdir checkout/result")

fp        = open("checkout.inp","r")

fcontrol  = open("checkout/control.raw.inp","w")
[part_top,part_variable,part_init,part_point_implicit,part_primitive,part_main] = 6*[""]
[part_Makefile_main,part_Makefile_var] = 2*[""]

####################################################################################
### CORE PARTS ###
engage("core")

# process n_grid.f90
filename = read_control_next(fp)

fgrid = open(filename,'r')
fgrid.readline()
nijk = fgrid.readline().strip().split()
fgrid.close()
for i in range(0,3):
	nijk[i] = int(nijk[i])-1
if(nijk[2] != 0):
	print "This program cannot be applied yet to 3D problem."
	print nijk[2]
	sys.exit(1)
os.system("cp "+filename+" checkout/grid/")
nY  = read_control_next(fp)
fromto = [\
	["NumI",str(nijk[0])],\
	["NumJ",str(nijk[1])],\
	["NumY",nY],\
	["GridFileName",os.path.basename(filename)]]
raw2pro("checkout/n_grid.raw.f90","checkout/n_grid.f90",fromto)
####################################################################################
print "***architecture***"
val = read_control_next_int(fp)
if(val == 0):
	print "\tPersonal Computer is selected."
	engage("architecture/PC")

elif(val >= 1):
	print "\tM-System[Super Computer] is selected."
	engage("architecture/SC")

	import checkout.generate_proc_grid
	[MPINumGrid,bwmax]=checkout.generate_proc_grid.split_grid_from_Nproc(val,nijk)
	Nproc = MPINumGrid[0]*MPINumGrid[1]
	import shutil
	shutil.move("grid_separation.inp","checkout/")
	raw2pro("checkout/mpirun.raw.sh","checkout/mpirun.sh",[["Nproc",str(Nproc)]])
	fromto = [\
		 ["NGX",str(MPINumGrid[0])],\
		 ["NGY",str(MPINumGrid[1])],\
		 ["BWMAX",str(bwmax)]]
	raw2pro("checkout/mod_mpi.raw.f90","checkout/mod_mpi.f90",fromto)
	print "\tActual NumProc is set to %3i x %3i = %3i" % (MPINumGrid[0],MPINumGrid[1],Nproc)
else:
	print "\tOdd Input at architecture! value is ",val
	sys.exit(1)

####################################################################################
print "***time scale***"
val_dt = read_control_next_int(fp)
if(val_dt == 0):
	print "\tGlobal Time Step is selected."
	engage("time_step/global")
	DT_LOCAL_GLOBAL = "dt_grbl"
elif(val_dt == 1):
	print "\tLocal Time Step is selected."
	engage("time_step/local")
	DT_LOCAL_GLOBAL = "dt_mat(i,j)"
else:
	print "\tOdd Input at time step! value is ",val_dt
	sys.exit(1)

####################################################################################
print "***time scheme***"
val = read_control_next_int(fp)
if(val == 0):
	print "\tEuler is selected."
	engage("time_scheme/euler")
	raw2pro("checkout/time_euler_viscous.raw.f90","checkout/time_euler_viscous.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
elif(val == 1):
	print "\ttwo-step Runge-Kutta is selected."
	engage("time_scheme/RK2")
	raw2pro("checkout/time_RK2_viscous.raw.f90","checkout/time_RK2_viscous.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
elif(val == 2):
	print "\tLU-SGS is selected."
	engage("time_scheme/LU-SGS")
	raw2pro("checkout/sch_lusgs.raw.f90","checkout/sch_lusgs.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
elif(val == 3):
	print "\tNR is selected."
	engage("time_scheme/NR")
	raw2pro("checkout/sch_NR.raw.f90","checkout/sch_NR.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
elif(val == 4):
	print "\tPre-conditioner is selected."
	engage("time_scheme/precon")
	raw2pro("checkout/sch_precon.raw.f90","checkout/sch_precon.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
	if(val_dt == 0):
		print "ERROR: This scheme cannot be used at global time step."
		sys.exit(1)
elif(val == 5):
	print "\tDual Time is selected."
	engage("time_scheme/dual")
	raw2pro("checkout/sch_dual.raw.f90","checkout/sch_dual.f90",[["DT_LOCAL_GLOBAL",DT_LOCAL_GLOBAL]])
	if(val_dt == 1):
		print "ERROR: This scheme cannot be used at local time step."
		sys.exit(1)
elif(val == 6):
	print "\tPre-conditioner using LU-SGS is selected."
	engage("time_scheme/preconLU-SGS")
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
	engage("high_order/upwind")
elif(val == 1):
	print "\tMUSCL is selected."
	engage("high_order/muscl")
else:
	print "\tOdd Input at high order scheme! value is ",val
	sys.exit(1)

####################################################################################
print "***viscosity***"
val = read_control_next_int(fp)
if(val == 0):
	print "\tnon-viscous flow is selected."
	engage("viscosity/non-viscous")
elif(val == 1):
	print "\tviscous flow is selected."
	engage("viscosity/viscous")
else:
	print "\tOdd Input at viscosity! value is ",val
	sys.exit(1)

####################################################################################
print "***thermal model***"
val = read_control_next_int(fp)
if(val == 0):
	print "\tideal gas is selected."
	engage("therm_lib/ideal")

	# process thermal_model.f90
	val = read_control_next_split(fp,3,"thermal model")
	fromto = [ \
		["KAPPA" ,val[0]],\
		["RGAS"  ,val[1]],\
		["NU"    ,val[2]]]
	raw2pro("checkout/thermal_model.raw.f90","checkout/thermal_model.f90",fromto)
elif(val == 1 or val == 2):
	if(val == 1) :
		print "\tCold Flow using chemical kinetics database is selected."
	else:
		print "\tChemical Kinetics Model is selected."
		if(val_dt == 1):
			print "Error : Chemical Kinetics Model cannot be used at local time step."
			sys.exit(1)
		engage("therm_lib/kinetics_reaction")

	engage("therm_lib/kinetics")


	# process mod_chem.f90
	val = read_control_next_split(fp,2,"thermal model")
	fromto = [ \
		["NumOfElements " ,val[0]],\
		["NumOfSpecies  " ,nY],\
		["NumOfReactions" ,val[1]]]
	raw2pro("checkout/mod_chem.raw.f90","checkout/mod_chem.f90",fromto)
else:
	print "\tOdd Input at thermal model! value is ",val
	sys.exit(1)


####################################################################################
# close checkout.inp
fp.close()

# close control.raw.inp
fcontrol.close()

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
append_file_to_file(fmain,"checkout/main.body.f90")
fromto = [\
["part_top\n"           , part_top],\
["part_init\n"          , part_init],\
["part_point_implicit\n", part_point_implicit],\
["part_primitive\n"     , part_primitive],\
["part_main\n"          , part_main]]
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
if os.path.exists("checkout/__init__.py"):
	os.remove("checkout/__init__.py")
os.system("rm -f checkout/*.pyc")

