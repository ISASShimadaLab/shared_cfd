#!/usr/bin/python
def plot_grid(nijk):
	space=20 #space between planes
	#file read
	fi = open("grid_separation.inp")
	fo = open("vis_grid.dat","w")

	isum = 0
	for nijk_mono in nijk:
		ni=nijk_mono.pop(0)
		nj=nijk_mono.pop(0)

		line  = fi.readline()
		arr = line.split()
		for var in arr:
		   var=int(var)
		   fo.write("%3.3i %3.3i\n" % (var+isum,1))
		   fo.write("%3.3i %3.3i\n" % (var+isum,nj))
		   fo.write("\n")
		fo.write("%3.3i %3.3i\n" % (ni+isum,1))
		fo.write("%3.3i %3.3i\n" % (ni+isum,nj))
		fo.write("\n")

		line  = fi.readline()
		arr = line.split()
		for var in arr:
		   var=int(var)
		   fo.write("%3.3i %3.3i\n" % (1 +isum,var))
		   fo.write("%3.3i %3.3i\n" % (ni+isum,var))
		   fo.write("\n")
		fo.write("%3.3i %3.3i\n" % (1 +isum,nj))
		fo.write("%3.3i %3.3i\n" % (ni+isum,nj))
		fo.write("\n")

		isum +=ni+space
	fi.close()
	fo.close()

def makelist(arr,ng):
	sum_r=0
	sum_a=0
	for i in range(len(arr)):
		if arr[i]<0:
			sum_r -= arr[i]
		else:
			sum_a += arr[i]
	
	unit = (ng-sum_a)/sum_r
	res  = (ng-sum_a)%sum_r
	
	now=1
	line=""
	max_bw = 0
	for i in range(len(arr)):
		line+=("%3.3i " % now)
		if (i+1)%5 == 0:
			line+="   "
		
		if arr[i]<0:
			delta = -arr[i]*unit
		else:
			delta =  arr[i]
		
		if i<res:
			delta+=1

		if delta>max_bw:
			max_bw = delta

		now+=delta
	return [line,max_bw]

def generate_grid_separation(nijk,MPINumGrid):
	max_bw=0
	f=open("grid_separation.inp","w")

	for plane,nijk_mono in enumerate(nijk):
		MPINumGrid_mono = MPINumGrid[plane]
		for cind,var_n in enumerate(nijk_mono):
			var_g = MPINumGrid_mono[cind]
			arr=var_g*[-1]

			import os
			if os.path.exists("split_ratio.py"):
				import split_ratio
				arr=split_ratio.ratio(plane,cind)
				if len(arr) != var_g:
					print "Odd array at split_ratio.py"
					sys.exit(1)

			[line,max_bw_candidate] = makelist(arr,var_n)
			max_bw = max(max_bw,max_bw_candidate)
			f.write(line+"\n")
	f.close()
	
	plot_grid(nijk)

	return max_bw

def split_grid(bw,nijk):
	MPINumGrid=[]
	for plane in nijk:
		MPINumGrid_mono = []
		for var in plane:
			MPINumGrid_mono.append(int((var-1)/bw+1))
		MPINumGrid.append(MPINumGrid_mono)
	return MPINumGrid

def split_grid_from_Nproc(Nproc,nijk):
	if Nproc <len(nijk):
		print "Error: At least more than the number of plane is necessary for processors."
		import sys
		sys.exit(1)
	bw = 1
	while 1:
		MPINumGrid = split_grid(bw,nijk)
		Nproc_now = 0

		for arr in MPINumGrid:
			Nproc_plane = 1
			for var in arr:
				Nproc_plane *= var
			Nproc_now += Nproc_plane

		if(Nproc_now<=Nproc):
			break
		bw+=1
	max_bw = generate_grid_separation(nijk,MPINumGrid)
	return [MPINumGrid,max_bw]

if __name__ == "__main__":
	# set nijk
	nijk =[]

#	# MB_test00.x
#	# 1st plane
#	nijk_plane=[]
#	nijk_plane.append(2)
#	nijk_plane.append(2)
#	nijk.append(nijk_plane)
#
#	# 2nd plane
#	nijk_plane=[]
#	nijk_plane.append(2)
#	nijk_plane.append(2)
#	nijk.append(nijk_plane)

	# MB_test2.x
	# 1st plane
	nijk_plane=[]
	nijk_plane.append(50)
	nijk_plane.append(100)
	nijk.append(nijk_plane)

	# 2nd plane
	nijk_plane=[]
	nijk_plane.append(50)
	nijk_plane.append(50)
	nijk.append(nijk_plane)

#	# MB_test4.x
#	# 1st plane
#	nijk_plane=[]
#	nijk_plane.append(40)
#	nijk_plane.append(10)
#	nijk.append(nijk_plane)
#
#	# 2nd plane
#	nijk_plane=[]
#	nijk_plane.append(40)
#	nijk_plane.append(40)
#	nijk.append(nijk_plane)
#
#	# 3rd plane
#	nijk_plane=[]
#	nijk_plane.append(40)
#	nijk_plane.append(70)
#	nijk.append(nijk_plane)


	Nproc = 20
	print split_grid_from_Nproc(Nproc,nijk)
