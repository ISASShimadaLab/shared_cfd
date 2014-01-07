#!/usr/bin/python
import os
import sys
from store.checkout.checkout_flow import *
from store.checkout.checkout_chem_nasa import *
from store.checkout.checkout_model import *

flag=[0,0,0]
if   os.path.exists("checkout.inp"):		flag[0]=1
elif os.path.exists("checkout_chem.inp"):	flag[1]=1
elif os.path.exists("checkout_model.inp"):	flag[2]=1

s=sum(flag)
if s ==0:
	print "You have to make input file."
elif s>1:
	print "More than 1 input file."
	print "Remove unnecessary file."

if   flag[0]==1:checkout_flow()	
if   flag[1]==1:checkout_chem_nasa()
if   flag[2]==1:checkout_model()
os.system("rm -f store/*.pyc store/checkout/*.pyc")
