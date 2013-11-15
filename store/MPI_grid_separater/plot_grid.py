#!/usr/bin/python

#parameter
ni = 99
nj = 49 

#file read
fi = open("grid_separation.dat")
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

