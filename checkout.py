#!/usr/bin/python
import os
import sys
from store.checkout.checkout_flow import *
from store.checkout.checkout_nasa import *

argv = sys.argv
var=""
if len(argv) >1: var = argv[1]

if var == "nasa":
	checkout_nasa()
else:
	checkout_flow()	
os.system("rm -f store/*.pyc store/checkout/*.pyc")
