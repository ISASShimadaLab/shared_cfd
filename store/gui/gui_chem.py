#!/usr/bin/python
#coding:utf-8
import wx
import os
import wx.lib.mixins.listctrl  as  listmix
from control_buttons import *
from subprocess import Popen,PIPE

# global variables
[RA,cb,frame,SPCcand]=4*[None]

###############################################################
############################ CONTROL_BUTTONS_FLOW #############
###############################################################
class CONTROL_BUTTONS_CHEM(CONTROL_BUTTONS):
	def __init__(self,RightPanel):
		CONTROL_BUTTONS.__init__(self,RightPanel,event_handler)
	def enable_generate_button(self,RA):
		# calc logical product
		lgcsum = True
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
					majorDimension=1)
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
def data_set():
	datRA=[]
	fp = open("store/gui/chem.inp")

	fp.readline()
	while True:
		line=fp.readline().rstrip()
		if line[:3] =="---": break
		arr=[line]
		arr.append(fp.readline().rstrip().split(","))
		arr.append(fp.readline().rstrip())
		arr.append(fp.readline().rstrip().replace("\\n","\n"))
		datRA.append(arr)
		fp.readline()

	fp.close()

	return datRA

###############################################################
#################### For First Panel Generation ###############
###############################################################
def set_panel(mainframe,base_panel):
	global RA,cb,frame,SPCcand
	frame=mainframe
	RightPanel = wx.Panel(base_panel,wx.ID_ANY)
	RightPanel.SetBackgroundColour('#FFFFFF')
	layoutRight = wx.BoxSizer(wx.VERTICAL)

	#data set
	datRA=data_set()

	# Radio Buttons
	RA=[]
	for i in range(len(datRA)):
		RA.append(RAelm(i,RightPanel,datRA[i]))
		layoutRight.Add(RA[i],		flag=wx.GROW|wx.ALL,border=10)

	# Control Buttons
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
# general interface
def event_handler(e):
	ID = e.GetId()
	typeID = ID%10;partID = ID/10
	process_values(typeID,partID)

def process_values(typeID,partID):
	if(typeID == 9): # close button
		if   partID == 98:
			output_and_run()
		else:
			frame.Close()

	cb.enable_generate_button(RA)

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

