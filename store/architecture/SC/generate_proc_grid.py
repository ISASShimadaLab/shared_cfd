#!/usr/bin/python
def plot_grid(ni,nj):
	#file read
	fi = open("grid_separation.inp")
	fo = open("vis_grid.dat","w")
	
	line  = fi.readline()
	arr = line.split()
	ng = len(arr)
	for i in range(ng):
	   tmp=int(arr[i])
	   fo.write("%3.3i %3.3i\n" % (tmp,1))
	   fo.write("%3.3i %3.3i\n" % (tmp,nj))
	   fo.write("\n")
	fo.write("%3.3i %3.3i\n" % (ni,1))
	fo.write("%3.3i %3.3i\n" % (ni,nj))
	fo.write("\n")
	
	line  = fi.readline()
	arr = line.split()
	ng = len(arr)
	for j in range(ng):
	   tmp=int(arr[j])
	   fo.write("%3.3i %3.3i\n" % (1 ,tmp))
	   fo.write("%3.3i %3.3i\n" % (ni,tmp))
	   fo.write("\n")
	fo.write("%3.3i %3.3i\n" % (1 ,nj))
	fo.write("%3.3i %3.3i\n" % (ni,nj))
	fo.write("\n")
	
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

def generate_grid_separation(ni,nj,ngx,ngy):
	arr_x=ngx*[-1]
	arr_y=ngy*[-1]

	import os
	if os.path.exists("split_ratio.py"):
		import split_ratio
		arr_x=split_ratio.ratio_x()
		arr_y=split_ratio.ratio_y()
		if len(arr_x) != ngx or len(arr_y) != ngy:
			print "Odd array at split_ratio.py"
			sys.exit(1)

	f=open("grid_separation.inp","w")
	[line,max_bw_x] = makelist(arr_x,ni)
	f.write(line+"\n")
	[line,max_bw_y] = makelist(arr_y,nj)
	f.write(line+"\n")
	f.close()
	
	plot_grid(ni,nj)

	return max(max_bw_x,max_bw_y)

def split_grid(bw,nijk):
	MPINumGrid = 3*[0]
	for i in range(0,3):
		MPINumGrid[i]  = int((nijk[i]-1)/bw+1)
	return MPINumGrid

def split_grid_from_Nproc(Nproc,nijk):
	bw = 1
	while 1:
		MPINumGrid = split_grid(bw,nijk)
		Nproc_now = 1
		for i in range(0,3):
			if(MPINumGrid[i] != 0):
				Nproc_now *= MPINumGrid[i]
		if(Nproc_now<=Nproc):
			break
		bw+=1
	max_bw = generate_grid_separation(nijk[0],nijk[1],MPINumGrid[0],MPINumGrid[1])
	return [MPINumGrid,max_bw]

if __name__ == "__main__":
	import sys
	array = sys.argv
	del array[0]
	if   len(array) == 2:
		print "plot_grid"
		[ni,nj] = map(int,array)
		plot_grid(ni,nj)
	if   len(array) == 3:
		print "split_grid_from_Nproc"
		nijk=3*[0]
		[nijk[0],nijk[1],Nproc] = map(int,array)
		MPINumGrid = split_grid_from_Nproc(Nproc,nijk)
		print MPINumGrid
	elif len(array) == 4:
		print "generate_grid_separation"
		[ni,nj,ngx,ngy] = map(int,array)
		print generate_grid_separation(ni,nj,ngx,ngy)
	else:
		print "Odd Command Parameters."
		sys.exit(1)

