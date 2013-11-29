#!/usr/bin/python
import sys
from mod_cond import *

def set_touch(BD):
	touch=[]
	i=0
	while i <len(BD)-1:
	   BDi=BD[i]
	   BDI=BD[i+1]
	   if not BDi.isSameLine(BDI):
	      i+=1
	      continue
	
	   for j in range(i+1,len(BD)):
	      BDj=BD[j]
	      if BDi.grid == BDj.grid:
	         pm=0
	         if j+1<len(BD) and BDj.isSameLine(BD[j+1]) and BDI.grid   == BD[j+1].grid:
	            pm=+1
	         elif               BDj.isSameLine(BD[j-1]) and BDI.grid   == BD[j-1].grid:
	            pm=-1
	
	         if pm != 0:
	            touch_elm = TOUCH_elm(BDi,BDj,pm)
	            width=1
	            while True:
	               if(j+width*pm >=len(BD) or i+width >=len(BD)):
	                  break
	
	               BDI=BD[i+width]
	               BDJ=BD[j+width*pm]
	               if BDi.isSameLine(BDI) and BDj.isSameLine(BDJ) and BDI.grid   == BDJ.grid:
	                  width+=1
	               else:
	                  break
	            touch_elm.width=width-1
	            i+=width-2
	            touch.append(touch_elm)
	   i+=1

	return touch


def set_toSetManually(Nplane,nijk,touch):
	# Initialize
	toSetManually=[]
	for plane in range(Nplane):
	   toSetManually.append([])
	   for cind in range(2):
	      if cind ==0:
	         n=nijk[plane][1]-1
	      if cind ==1:
	         n=nijk[plane][0]-1
	      toSetManually[plane].append([[[1,n]],[[1,n]]])
	
	# generate
	for var in touch:
	   for lu in range(2):
	      if lu ==0:
	         toRemove = [var.s1,var.s1+var.width-1]
                 plane = var.p1
                 cind  = var.cind1
                 lu    = var.lu1
	      else:
	         toRemove = [var.s2,var.s2+var.order*(var.width-1)]
	         toRemove.sort()
                 plane = var.p2
                 cind  = var.cind2
                 lu    = var.lu2
	
	      now_work=toSetManually[plane][cind][lu]
	
	      for i,elm in enumerate(now_work):
	         if elm[0] <= toRemove[0] and toRemove[1] <= elm[1]:
	            tmp =[]
	            if elm[0]<toRemove[0]:
	               tmp.append([elm[0],toRemove[0]-1])
	            if elm[1]>toRemove[1]:
	               tmp.append([toRemove[1]+1,elm[1]])
	            now_work[i:i+1]=tmp
	            toRemove=None
	            break
	
	      if toRemove!= None:
	         print "Unexpected Situation Occured!"
	         print "STOP!"
	         sys.exit(1)
	
#	print "toSetManually=",toSetManually

	return toSetManually
