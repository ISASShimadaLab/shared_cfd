#/usr/bin/python
import os
import sys
import subprocess

arr= sys.argv
if(len(arr)<2):
	p = subprocess.Popen(["qstat"],stdout=subprocess.PIPE)
	i=0
	for line in p.stdout:
		if "y535" in line:i+=1
	num=10-i
	print num,"job to throw"
else:
	num=int(arr[1])

for var in range(num):
	os.system("sh vecrun.sh ")

