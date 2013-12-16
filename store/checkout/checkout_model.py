#!/usr/bin/python
import os
import sys
from checkout_utility import *
from preCEA import *

def checkout_model():
	fp        = open("checkout_model.inp","r")
	
	val = read_control_next_int(fp)

	if(val == 0 or val == 1):
		preCEA(val,fp)
	else:
		print "\tOdd Input at thermal model! value is ",val
		sys.exit(1)

	fp.close()
