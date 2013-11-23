#!/usr/bin/python
import sys

class SPCelm:
	def __init__(self,name,elm_dic):
		self.name    = name
		self.elm_dic = elm_dic
class CONTROL:
	def __init__(self,fuel,oxid,prod):
		self.fuel = fuel
		self.oxid = oxid
		self.prod = prod
	def species(self):
		arr = []
		arr.extend(self.fuel.keys())
		arr.extend(self.oxid.keys())
		arr.extend(self.prod.keys())
		return list(set(arr))
	def set_species(self,species):
		self.species = species
	def elements(self):
		arr = []
		for SPC in self.species:
			arr.extend(SPC.elm_dic.keys())
		return list(set(arr))
	def check_consistency():

def read_thermo():
	fp = open("thermo.inp")
	while True:
		if fp.readline()[:1] != '!': break
	fp.readline()
	
	SPClist = []
	while True:
		name = fp.readline()[:18].strip()
		if name[:3] == "END": break
		line = fp.readline()
		Nsctn     = int(line[:2])
		phase     = int(line[51:52])
		elem_line = line[10:50]
		elm_dic={}
		for i in range(5):
			elm       = elem_line[ :2]
			amount    = elem_line[2:8].strip()
			elem_line = elem_line[8:]
			if(elm == '  '):break
			elm_dic[elm]=float(amount)
		for i in range(Nsctn*3):fp.readline()
		if phase ==0: SPClist.append(SPCelm(name,elm_dic))
	fp.close()
	return SPClist

def read_control_chem():
	fp = open("control_chem.inp")

	for i in range(4): fp.readline()
	fuel = {}
	while True:
		line = fp.readline().strip()
		if(line == 'end'):break
		[key,val] = line.split()
		fuel[key] = float(val)

	for i in range(3): fp.readline()
	oxid = {}
	while True:
		line = fp.readline().strip()
		if(line == 'end'):break
		[key,val] = line.split()
		oxid[key] = float(val)

	fp.readline()
	prod = {}
	while True:
		line = fp.readline().strip()
		if(line == 'end'):break
		[key,val] = line.split()
		prod[key] = float(val)

	fp.close()
	return CONTROL(fuel,oxid,prod)

def fetch_elm_info(SPClist,control):
	controlSPC = control.species()
	arr = []
	for SPCmono in SPClist:
		for species in controlSPC:
			if(SPCmono.name == species):
				arr.append(SPCmono)
				break
	control.set_species(arr)
	return control
	


if __name__ == "__main__":
	SPClist = read_thermo()
	control = read_control_chem()
	control = fetch_elm_info(SPClist,control)
	print control.elements()

