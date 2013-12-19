#!/usr/bin/python
import wx
import os
from control_buttons import *
from subprocess import Popen,PIPE

# global variables
[NG,RA,cb,frame]=4*[None]

###############################################################
############################ CONTROL_BUTTONS_CHEM #############
###############################################################
class CONTROL_BUTTONS_CHEM(CONTROL_BUTTONS):
	def __init__(self,RightPanel):
		CONTROL_BUTTONS.__init__(self,RightPanel,event_handler)
	def enable_generate_button(self):
		# calc logical product
		lgcsum=True
		for var in RA:
			lgcsum &= var.isValidate
	
		# Enable of disable button_generate
		if(lgcsum):
			self.button_generate.Enable()
		else:
			self.button_generate.Disable()

###############################################################
############################ RA ###############################
###############################################################
class RAelm(wx.RadioBox):
	def __init__(self,i,RightPanel,datRA):
		wx.RadioBox.__init__(	self,\
					RightPanel,\
					i*10+3,datRA[0],\
					choices=datRA[1],\
					style = wx.RA_HORIZONTAL,\
					majorDimension=3)
		self.Bind(wx.EVT_RADIOBOX,event_handler)
		self.isValidate = True
		self.inphead = datRA[2]
		self.inptail = datRA[3]
		self.disabledBy = [[] for ii in range(self.GetCount())]
	def inputstring(self):
		return self.inphead+str(self.GetSelection())+self.inptail+"\n"

###############################################################
####################### I/O ###################################
###############################################################
def data_set_flow():
	[datRA,NG]=[[],[]]
	fp = open("store/gui/chem.inp")

	fp.readline()
	while True:
		line=fp.readline().rstrip()
		if line[:3] =="---": break
		arr=[line]
		arr.append(fp.readline().rstrip().split(","))
		for i in range(2): arr.append(fp.readline().rstrip())
		datRA.append(arr)
		fp.readline()

	for line in fp:
		NG.append(map(int,line.rstrip().split(',')))

	fp.close()

	return [datRA,NG]

###############################################################
#################### For First Panel Generation ###############
###############################################################
def set_panel(mainframe,base_panel):
	global NG,RA,cb,frame
	frame=mainframe
	RightPanel = wx.Panel(base_panel,wx.ID_ANY)
	RightPanel.SetBackgroundColour('#FFFFFF')
	layoutRight = wx.BoxSizer(wx.VERTICAL)

	#data set
	[datRA,NG]=data_set_flow()

	# Radio Buttons
	RA=[]
	for i in range(len(datRA)):
		RA.append(RAelm(i,RightPanel,datRA[i]))
		layoutRight.Add(RA[i],flag=wx.GROW,border=10)

	cb=CONTROL_BUTTONS_CHEM(RightPanel)
	layoutRight.Add(cb.panel,flag=wx.ALIGN_RIGHT)
	RightPanel.SetSizer(layoutRight)

	#initialize
	for i,v in enumerate(RA):
		v.SetSelection(0)
		process_values(3,i)

	return RightPanel

###############################################################
#################### Event Handlers ###########################
###############################################################
def checkvalues(typeID,partID): #used only at process_values
	#remove disability by old (typeID,partID)
	for i,v in enumerate(RA):
		for j in range(RA[i].GetCount()):
			for k in range(len(RA[i].disabledBy[j])-1,-1,-1):
				# if partID disables [i,j]
				if RA[i].disabledBy[j][k] == partID:
					del RA[i].disabledBy[j][k]
					# if no selection denies selecting "item j"
					if len(RA[i].disabledBy[j]) == 0 :
						RA[i].EnableItem(j,True)
						# if selected item now is "item j"
						if(RA[i].GetSelection() == j):
							RA[i].isValidate = True

	#set disability caused by new (typeID,partID)
	RA[partID].isValidate = True
	value = RA[partID].GetSelection()

	for NGelm in NG:
		KeyValue = NGelm[partID]
		if KeyValue == -1:continue
		if KeyValue != value:continue

		for i,BadValue in enumerate(NGelm):
			if partID == i:continue
			if BadValue != -1:
				RA[i].EnableItem(BadValue,False)
				RA[i].disabledBy[BadValue].append(partID)
				if RA[i].GetSelection() == BadValue:
					RA[i].isValidate = False

# general interface
def event_handler(e):
	ID = e.GetId()
	typeID = ID%10;partID = ID/10
	process_values(typeID,partID)

def process_values(typeID,partID):
	if   typeID == 1 : # FileOpen
		fn.OnOpen()
	elif typeID == 3 :
		checkvalues(typeID,partID)
	elif typeID == 9 : # close button
		if   partID == 98:
			output_and_run()
		else:
			frame.Close()

	cb.enable_generate_button()

###############################################################
#################### For Execution ############################
###############################################################
def output_and_run():
	os.system("rm -f checkout.inp")
	os.system("rm -f checkout_model.inp")

	# output checkout.inp
	finp = open('checkout_chem.inp','w')
	for var in RA:
		finp.write(var.inputstring())
	finp.close()

	# run program
	p = Popen("./checkout.py",stdout = PIPE)
	printedout = p.stdout.read()
	value = p.wait()
	wx.MessageBox("Output from checkout.py\n"\
			+"**************************\n"\
			+printedout\
			+"**************************\n"\
			+"exit value: "+str(value),"Message")
	if value==0:frame.Close()

