#!/usr/bin/python
import sys
import os
from numpy import *
from mod_cond import *

def read_plot3d(filename):
	fp = open(filename,'r')
	
	buf=fp.readline().split()
	Nplane=int(buf.pop(0))
	nijk=[]
	BD=[]
	
	for i in range(Nplane):
	   buf=fp.readline().split()
	   for i,var in enumerate(buf):
	      buf[i]=int(var)
	   nijk.append(buf)
	
	for plane in range(Nplane):
	   ng=1
	   for var in nijk[plane]:
	      ng*=var
	
	   ### read data ###
	   prearr=[]
	   for i in range(int((ng+3)/4)):
	      buf=fp.readline().split()
	      for j,var in enumerate(buf):
	         buf[j]=float(var)
	      prearr.extend(buf)
	   prearr = array(prearr)
	   gridx = prearr.reshape(nijk[plane][1],nijk[plane][0]).transpose()
	
	   prearr=[]
	   for i in range(int((ng+3)/4)):
	      buf=fp.readline().split()
	      for j,var in enumerate(buf):
	         buf[j]=float(var)
	      prearr.extend(buf)
	   prearr = array(prearr)
	   gridy = prearr.reshape(nijk[plane][1],nijk[plane][0]).transpose()
	
	   prearr=[]
	   for i in range(int((ng+3)/4)):
	      buf=fp.readline().split()
	      for j,var in enumerate(buf):
	         buf[j]=float(var)
	      prearr.extend(buf)
	   prearr = array(prearr)
	   buf = prearr.reshape(nijk[plane][1],nijk[plane][0]).transpose()
	
	   ### get only boundaries ###
	   ni=nijk[plane][0]
	   nj=nijk[plane][1]
	   # i-plane
	   for j in range(nijk[plane][1]):
	        BD.append(BD_elm([gridx[   0][j],gridy[   0][j]],plane,0,0,j))
	   for j in range(nijk[plane][1]):
	        BD.append(BD_elm([gridx[ni-1][j],gridy[ni-1][j]],plane,0,1,j))
	
	   # j-plane
	   for i in range(nijk[plane][0]):
	        BD.append(BD_elm([gridx[i][   0],gridy[i][   0]],plane,1,0,i))
	   for i in range(nijk[plane][0]):
	        BD.append(BD_elm([gridx[i][nj-1],gridy[i][nj-1]],plane,1,1,i))
	
	fp.close()
	return [Nplane,nijk,BD]
