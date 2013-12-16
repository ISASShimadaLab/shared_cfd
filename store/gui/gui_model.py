#!/usr/bin/python
#coding:utf-8
import wx
import os
import wx.lib.mixins.listctrl  as  listmix
from ..checkout import preCEA
from control_buttons import *
from subprocess import Popen,PIPE
 
data = [
"H2O",
"O2",
"LOX",
]
# global variables
[RA,LB,cb,frame,SPCcand]=5*[None]

###############################################################
############################ CONTROL_BUTTONS_FLOW #############
###############################################################
class CONTROL_BUTTONS_MODEL(CONTROL_BUTTONS):
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
					majorDimension=2)
		self.Bind(wx.EVT_RADIOBOX,event_handler)
		self.isValidate = True
		self.inphead = datRA[2]
		self.inptail = datRA[3]
		self.disabledBy = [[] for ii in range(self.GetCount())]
	def inputstring(self):
		return self.inphead+str(self.GetSelection())+self.inptail+"\n"
  
###############################################################
############################ LB ###############################
###############################################################
class LBplane(wx.Panel):
	def __init__(self,RightPanel):
		size = (150,380)
		wx.Panel.__init__(self,RightPanel,wx.ID_ANY)
		self.SetBackgroundColour('#FFFFFF')

		tmpbox    = wx.StaticBox(self,wx.ID_ANY,"chemical species selection")
		tmplayout = wx.StaticBoxSizer(tmpbox,wx.HORIZONTAL)

		tmppanel1=wx.Panel(self,wx.ID_ANY)
		tmppanel1.SetBackgroundColour("#FFFFFF")
		tmplayout1 = wx.BoxSizer(wx.VERTICAL)
		self.lb_select = wx.ListBox(tmppanel1,wx.ID_ANY,size=size,\
				style=wx.LB_SINGLE|wx.LB_SORT)
		self.lb_select. Bind(wx.EVT_LISTBOX_DCLICK,self.event_out)
		tmplayout1.Add(wx.StaticText(tmppanel1,wx.ID_ANY,"Selected"),\
				flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		tmplayout1.Add(self.lb_select,flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		tmppanel1.SetSizer(tmplayout1)

		tmppanel2=wx.Panel(self,wx.ID_ANY)
		tmppanel2.SetBackgroundColour("#FFFFFF")
		tmplayout2 = wx.BoxSizer(wx.VERTICAL)
		self.inbutton =wx.Button(tmppanel2,14,"<-----")
		self.outbutton=wx.Button(tmppanel2,24,"----->")
		self.inbutton. Bind(wx.EVT_BUTTON,self.event_in)
		self.outbutton.Bind(wx.EVT_BUTTON,self.event_out)
		tmplayout2.Add(self.inbutton, flag=wx.ALL |wx.ALIGN_CENTER,border=5)
		tmplayout2.Add(self.outbutton,flag=wx.ALL |wx.ALIGN_CENTER,border=5)
		tmppanel2.SetSizer(tmplayout2)

		tmppanel3=wx.Panel(self,wx.ID_ANY)
		tmppanel3.SetBackgroundColour("#FFFFFF")
		tmplayout3 = wx.BoxSizer(wx.VERTICAL)
		self.lb_notselect = wx.ListBox(tmppanel3,wx.ID_ANY,choices=SPCcand,\
					size=size,style=wx.LB_SINGLE|wx.LB_SORT)
		self.lb_notselect. Bind(wx.EVT_LISTBOX_DCLICK,self.event_in)
		tmplayout3.Add(wx.StaticText(tmppanel3,wx.ID_ANY,"Not Selected"),\
				flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		tmplayout3.Add(self.lb_notselect,flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		tmppanel3.SetSizer(tmplayout3)

		tmplayout.Add(tmppanel1,flag=wx.ALL |wx.ALIGN_CENTER,border=5)
		tmplayout.Add(tmppanel2,flag=wx.ALL |wx.ALIGN_CENTER,border=5)
		tmplayout.Add(tmppanel3,flag=wx.ALL |wx.ALIGN_CENTER,border=5)
		self.SetSizer(tmplayout)

		self.isValidate = True

	def event_in(self,e):
		select = self.lb_notselect.GetSelection()
		if select==-1: return
		self.lb_select.Append(self.lb_notselect.GetString(select))
		self.lb_notselect.Delete(select)

	def event_out(self,e):
		select = self.lb_select.GetSelection()
		if select==-1: return
		self.lb_notselect.Append(self.lb_select.GetString(select))
		self.lb_select.Delete(select)
	 
###############################################################
####################### I/O ###################################
###############################################################
def data_set():
	datRA=[]
	fp = open("store/gui/model.inp")

	fp.readline()
	while True:
		line=fp.readline().rstrip()
		if line[:3] =="---": break
		arr=[line]
		arr.append(fp.readline().rstrip().split(","))
		for i in range(2): arr.append(fp.readline().rstrip())
		datRA.append(arr)
		fp.readline()

	fp.close()

	SPCcand = preCEA.SPCcandidates()
	for i,v in enumerate(SPCcand):SPCcand[i]=v.name
	return [datRA,SPCcand]

###############################################################
#################### For First Panel Generation ###############
###############################################################
def set_panel(mainframe,base_panel):
	global RA,LB,cb,frame,SPCcand
	frame=mainframe
	RightPanel = wx.Panel(base_panel,wx.ID_ANY)
	RightPanel.SetBackgroundColour('#FFFFFF')
	layoutRight = wx.BoxSizer(wx.VERTICAL)

	#data set
	[datRA,SPCcand]=data_set()

	# Radio Buttons
	RA=[]
	for i in range(len(datRA)):
		RA.append(RAelm(i,RightPanel,datRA[i]))
		layoutRight.Add(RA[i],		flag=wx.GROW|wx.ALL,border=10)

	layoutRight.Add(LBplane(RightPanel),flag=wx.GROW|wx.ALL,border=10)

	# Control Buttons
	cb=CONTROL_BUTTONS_MODEL(RightPanel)
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
def checkvalues(typeID,partID):
	print typeID,partID

def event_handler(e):
	ID = e.GetId()
	typeID = ID%10;partID = ID/10
	process_values(typeID,partID)

def process_values(typeID,partID):
	if(  typeID == 1): # FileOpen
		fn.OnOpen()
	elif(typeID ==3):
		checkvalues(typeID,partID)
	elif(typeID == 9): # close button
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
	os.system("rm -f checkout_chem.inp")

	# output checkout.inp
	finp = open('checkout_model.inp','w')
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

