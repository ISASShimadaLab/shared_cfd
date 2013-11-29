#!/usr/bin/python
import sys
import os
def ckinterp(): 
	if not os.path.exists("chem.inp"):
		print "Can't find 'chem.inp'. Please try again."
		sys.exit(1)

	fp = open("chem.inp")
	#element
	ELM=[]
	while True:
		line = fp.readline()
		if line.startswith("!"):continue
		if line.startswith("elements"):continue
		if line.startswith("end"):break
		ELM.extend(line.split())
	#print ELM
	
	#species
	SPC=[]
	while True:
		line = fp.readline()
		if line.startswith("!"):continue
		if line.startswith("species"):continue
		if line.startswith("end"):break
		SPC.extend(line.split())
	#print SPC
	
	#thermo
	while True:
		line = fp.readline()
		if line.startswith("end"):break
	
	#reactions
	RCT=[]
	while True:
		line = fp.readline()
		if line.startswith("!"):continue
		if line.startswith("reactions"):continue
		if line.startswith("end"):break
		if "/" in line:
			RCT[-1]+=line
		elif "duplicate" in line.lower():
			RCT[-1]+=line
		else:
			RCT.append(line)
	#print RCT

	fp.close()
	return [len(ELM),len(SPC),len(RCT)]
if __name__ == "__main__":
	print ckinterp()
