#!/usr/bin/python
#########################################################
#########################################################
#######     preCEA                                 ######
#######     IN:  checkout_model.inp,thermo.inp     ######
#######     OUT: control_chem.inp, chem.inp        ######
#########################################################
#########################################################
import sys
DIR = "" #DIR to output

######################## CLASSES ##################################################
class SPCelm:
	# name:species name
	# elm_dic : dictionary of number of elements
	#	key->element name
	#	val->number of elements
	def __init__(self,name,elm_dic):
		self.name    = name
		self.elm_dic = elm_dic
class CONTROL:
	# species : list of species name used
	# elements list of elements name used
	# dis_species: dictionary
	#	key->species name
	#	val->dictionary
	#		key->element name
	#		key->number of elements
	def __init__(self,species):
		self.species = list(set(species))
	def set_species(self,dic_species):
		# set dic_species and elements from dictionary.
		self.dic_species = dic_species
		arr =[]
		for var in self.dic_species.values(): arr.extend(var.keys())
		self.elements = list(set(arr))
	
######################## FUNCTIONS ################################################
def read_thermo():
	#read thermo.inp to gather species name
	fp = open("store/therm_lib/NASA/core/thermo.inp")
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

def read_checkout_chem(ntype,fp):
	#read checkout_chem.inp and write control_chem.raw.inp

	species=[]
	line  = fp.readline()
	while True:
		line = fp.readline().strip()
		if(line == 'end'):break
		species.append(line)
	fp.close()

	fo = open(DIR+"control_chem.raw.inp","w")
	fo.write("#Control parameters for chemistry\n")
	fo.write("Oxidizer Pressure(Pa)   : 1.e5\n")
	fo.write("Oxidizer Temperature(K) : 300.\n")
	fo.write("Oxygen Composition(mole ratio):\n")
	for var in species:fo.write(var+" 1\n")
	fo.write("end\n")
	fo.write("Fuel Pressure(Pa)       : 1.e5\n")
	fo.write("Fuel Temperature(K)     : 300.\n")
	fo.write("Fuel Composition(mole ratio):\n")
	for var in species:fo.write(var+" 1\n")
	fo.write("end\n")
	if ntype==0:
		fo.write("Products Composition(mole ratio):\n")
		for var in species:fo.write(var+" 2\n")
		fo.write("end\n")
	fo.close()

	return CONTROL(species)

def fetch_elm_info(SPClist,control):
	#species -> elements in species
	controlSPC = control.species
	dic = {}
	for SPCmono in SPClist:
		for species in controlSPC:
			if(SPCmono.name == species):
				dic[species] = SPCmono.elm_dic
				break
	control.set_species(dic)
	return control

def trimSPC(SPClist,elements):
	# gather species existable
	trimed = []
	for SPC in SPClist:
		flag = True
		for elm in SPC.elm_dic.keys():
			if not elm in elements : flag = False
		if flag: trimed.append(SPC.name)
	return trimed

def out_chem(ELM,SPC):
	#write out chem.inp
	fp = open(DIR+"chem.inp","w")

	fp.write("elements\n")
	line=""
	for i,var in enumerate(ELM):
		line += " "+var.strip()
		if len(line)>50 or i+1 == len(ELM):
			fp.write(line+"\n")
			line=""
	fp.write("end\n")

	fp.write("species\n")
	line=""
	for i,var in enumerate(SPC):
		line += " "+var.strip()
		if len(line)>50 or i+1 == len(SPC):
			fp.write(line+"\n")
			line=""
	fp.write("end\n")

	fp.write("thermo\n")
	fp.write("end\n")
	fp.write("reactions cal/mole  moles\n")
	fp.write("end\n")

	fp.close()
######################## MAIN ################################################
def preCEA(var,fp):
	# var =0 : flame sheet
	# var =1 : cea
	# fp :file pointer of checkout_chem.inp
	
	# SPClist: list of SPCelm
	# control: instance of CONTROL
	# trimedSPC: list of species name the combustion model uses.
	
	SPClist = read_thermo()
	control = read_checkout_chem(var,fp)
	control = fetch_elm_info(SPClist,control)
	if var ==0:
		trimedSPC = control.species
	else:
		trimedSPC = trimSPC(SPClist,control.elements)
	out_chem(control.elements,trimedSPC)
def SPCcandidates():
	return read_thermo()
