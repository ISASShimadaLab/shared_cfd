#!/usr/bin/python
import os
import sys
from store.checkout.checkout_flow import *
from store.checkout.checkout_nasa import *

if os.path.exists("checkout.inp"):
	checkout_flow()	
elif os.path.exists("checkout_chem.inp"):
	checkout_nasa()
else:
	print "You have to make input file."
os.system("rm -f store/*.pyc store/checkout/*.pyc")
