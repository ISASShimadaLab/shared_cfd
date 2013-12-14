#!/usr/bin/python
import os
import sys
from checkout_utility import *

def checkout_model():
	fp        = open("checkout_model.inp","r")
	
	val = read_control_next_int(fp)

	if(val == 0):
		import preCEA
		preCEA.preCEA(0,fp)
	else:
		print "\tOdd Input at thermal model! value is ",val
		sys.exit(1)

	fp.close()
