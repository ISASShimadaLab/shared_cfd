#!/usr/bin/python
from mod_cond import *
import sys
def fetch_myid(Nplane,nijk):
#	print 
	Nproc = 0
	gsep =[]
	import os
	if not os.path.exists("checkout/grid_separation.inp"):
		Nproc = 1
		bmpi=[]
		for var in nijk:
			bmpi_mono=[]
			ngx = var[0] -1
			ngy = var[1] -1
			bmpi_mono.append([[[1,ngy,0]],[[1,ngy,0]]])
			bmpi_mono.append([[[1,ngx,0]],[[1,ngx,0]]])
			bmpi.append(bmpi_mono)
#			print bmpi_mono[0][0]
#			print bmpi_mono[0][1]
#			print bmpi_mono[1][0]
#			print bmpi_mono[1][1]
#		sys.exit(1)
		return [bmpi,Nproc]

	fp = open("checkout/grid_separation.inp","r")
	for var in nijk:
		gsep_mono=[]
		for n in var[0:2]:
			arr = map(int,fp.readline().strip().split())
			arr.append(n)
			gsep_mono.append(arr)
		gsep.append(gsep_mono)
	fp.close()

	bmpi=[]
	for gsep_mono in gsep:
		bmpi_mono=[[[],[]],[[],[]]]
		gsepx = gsep_mono[0]
		gsepy = gsep_mono[1]
		ngx = len(gsepx)-1
		ngy = len(gsepy)-1

		for i in range(ngx):
			bmpi_mono[1][0].append([gsepx[i],gsepx[i+1]-1,Nproc+i*ngy])
			bmpi_mono[1][1].append([gsepx[i],gsepx[i+1]-1,Nproc+i*ngy+ngy-1])

		for i in range(ngy):
			bmpi_mono[0][0].append([gsepy[i],gsepy[i+1]-1,Nproc+i])
			bmpi_mono[0][1].append([gsepy[i],gsepy[i+1]-1,Nproc+i+ngy*(ngx-1)])
		bmpi.append(bmpi_mono)
		Nproc+=ngx*ngy
#		print bmpi_mono[0][0]
#		print bmpi_mono[0][1]
#		print bmpi_mono[1][0]
#		print bmpi_mono[1][1]
#	sys.exit(1)
	return [bmpi,Nproc]

def modify_touch(touch_old,bmpi,nijk,Nproc):
	touch   = []
	MPIcomm = []
	for touch_elm in touch_old:
		order = touch_elm.order
		width = touch_elm.width
		p1    = touch_elm.p1   
		p2    = touch_elm.p2   
		cind1 = touch_elm.cind1
		cind2 = touch_elm.cind2
		lu1   = touch_elm.lu1  
		lu2   = touch_elm.lu2  
		s1    = touch_elm.s1   
		s2    = touch_elm.s2   
		pm1   = touch_elm.pm1  
		pm2   = touch_elm.pm2  
#		print "WIDTH=",width
		bmpi_local1 = bmpi[p1][cind1][lu1]
		bmpi_local2 = bmpi[p2][cind2][lu2]
		myid1 =  search_myid(bmpi_local1,s1)
		myid2 =  search_myid(bmpi_local2,s2)

		width_before=0
		inc = 0
#		print "s1=",s1+inc,"s2=",s2+inc*order,"inc",inc,"myidnow1=",myid1,"myidnow2=",myid2
#		sys.exit(1)
		for inc in range(1,width):
			myid_now1 =  search_myid(bmpi_local1,s1+inc)
			myid_now2 =  search_myid(bmpi_local2,s2+inc*order)
			if(myid_now1 != myid1 or myid_now2 != myid2):
				recordMPI(	inc,\
						myid1,\
						myid2,\
						touch_elm,\
						width_before,\
						touch,\
						MPIcomm)
				myid1=myid_now1
				myid2=myid_now2
				width_before=inc
#			print "s1=",s1+inc,"s2=",s2+inc*order,"inc",inc,"myidnow1=",myid_now1,"myidnow2=",myid_now2

		recordMPI(	width,\
				myid1,\
				myid2,\
				touch_elm,\
				width_before,\
				touch,\
				MPIcomm)
#		print

	return [touch,MPIcomm]

def recordMPI(inc,myid1,myid2,touch_elm,width_before,touch,MPIcomm):
	MPIelm = TOUCH_elm(touch_elm)
	MPIelm.s1    = touch_elm.s1+width_before
	MPIelm.s2    = touch_elm.s2+width_before*touch_elm.order
	MPIelm.width = inc-width_before
	if myid1 == myid2:
		MPIelm.myid  = myid1
		touch.append(MPIelm)
	else:
		MPIelm.myid1 = myid1
		MPIelm.myid2 = myid2
		MPIcomm.append(MPIelm)


def search_myid(bmpi_local,s1):
	for arr in bmpi_local:
		if(arr[0]<=s1 and s1<= arr[1]):
			return arr[2]
	print "Error at search_myid!!"
	print s1
	print bmpi_local
	sys.exit(1)
