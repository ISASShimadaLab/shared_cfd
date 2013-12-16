#!/usr/bin/python
import os
import sys

def append_file_to_file(fp,filename):
	ftmp = open(filename)
	fp.write(ftmp.read())
	ftmp.close()
	os.system("rm "+filename)

def file_to_str(filename):
	ftmp = open(filename)
	string = ftmp.read()
	ftmp.close()
	os.system("rm "+filename)
	return string

def raw2pro(name_raw,name_pro,fromto):
	fraw = open(name_raw,"r")
	fpro = open(name_pro,"w")
	for line in fraw:
		for element in fromto:
			line = line.replace(element[0],element[1])
		fpro.write(line)
	fraw.close()
	fpro.close()
	os.system("rm "+name_raw)

def read_control_next(fp):
	string = fp.readline()
	if '#' in string:
		string = string[string.index(':')+1:string.index('#')-1].strip()
	else:
		string = string[string.index(':')+1:].strip()
	return string

def read_control_next_int(fp):
	return int(read_control_next(fp))

def read_control_next_split(fp,num,place):
	val = read_control_next(fp).split()
	if(len(val) != num):
		print "Error: The number of parameters for "+place+" is odd."
		sys.exit(1)
	return val

def engage(place,to,arr_engage):
	os.system("cp store/"+place+"/* "+to+"/")
	for var in arr_engage:
		if(os.path.exists(var[0])):var[1]+= file_to_str(var[0])
