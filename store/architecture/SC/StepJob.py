#/usr/bin/python
import os
import sys

arr= sys.argv
if(len(arr)<2):
	num=24
else:
	num=int(arr[1])

for var in range(num):
	os.system("""
./mpirun.sh <<'''

'''
""")

