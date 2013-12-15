#!/usr/bin/python
import wx
import os
from control_buttons import *
from subprocess import Popen,PIPE

# global variables
[NG,fn,TC,RA,cb,frame]=6*[None]

###############################################################
############################ CONTROL_BUTTONS_FLOW #############
###############################################################
class CONTROL_BUTTONS_FLOW(CONTROL_BUTTONS):
	def __init__(self,RightPanel):
		CONTROL_BUTTONS.__init__(self,RightPanel,event_handler)
	def enable_generate_button(self,BtnFileName,TC,RA):
		# calc logical product
		lgcsum = BtnFileName.logical
		for var in TC:
			lgcsum &= var.isValidate
		for var in RA:
			lgcsum &= var.isValidate
	
		# Enable of disable button_generate
		if(lgcsum):
			self.button_generate.Enable()
		else:
			self.button_generate.Disable()

###############################################################
############################ TC ###############################
###############################################################
class TCelm(wx.TextCtrl):
	def __init__(self,i,RightPanel,datTC):
		tmpPanel = wx.Panel(RightPanel,wx.ID_ANY)
		tmpPanel.SetBackgroundColour('#FFFFFF')
		tmpbox    = wx.StaticBox(tmpPanel,wx.ID_ANY,datTC[0])
		tmplayout = wx.StaticBoxSizer(tmpbox,wx.HORIZONTAL)
		wx.TextCtrl.__init__(self,tmpPanel,i*10+2,datTC[2])
		self.Bind(wx.EVT_TEXT,event_handler)
		self.isValidate = True
		self.initial_value = datTC[2]
		self.inphead       = datTC[3]
		self.inptail       = datTC[4]
		tmplayout.Add(self,flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		tmplayout.Add(wx.StaticText(tmpPanel,wx.ID_ANY,datTC[1]),\
			flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		tmpPanel.SetSizer(tmplayout)
		self.panel=tmpPanel
	def inputstring(self):
		return self.inphead+self.GetValue()+self.inptail+"\n"

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
############################ FN ###############################
###############################################################
class FileName:
	def __init__(self,RightPanel):
		tmpPanel = wx.Panel(RightPanel,wx.ID_ANY)
		tmpPanel.SetBackgroundColour('#FFFFFF')
		tmpbox    = wx.StaticBox(tmpPanel,wx.ID_ANY,"Grid File Name")
		tmplayout = wx.StaticBoxSizer(tmpbox,wx.HORIZONTAL)
		self.btn = wx.Button(tmpPanel,1,"Choose File.")
		self.btn.Bind(wx.EVT_BUTTON,event_handler)
		tmplayout.Add(self.btn,flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		self.lbl =wx.StaticText(tmpPanel,wx.ID_ANY,"No File Selected yet.") 
		self.inphead = "grid file name    : "
		self.inptail = " #Grid file name (*.x)"
		tmplayout.Add(self.lbl,flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		tmpPanel.SetSizer(tmplayout)
		self.logical = False
		self.panel=tmpPanel
	def inputstring(self):
		return self.inphead+self.fullpath+self.inptail+"\n"

	def OnOpen(self):
		dirName = os.getcwd()
		dialog = wx.FileDialog(frame, "Choose the grid file", dirName, "", "*.x", wx.OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			fileName = dialog.GetPath()
			self.fullpath = fileName
			self.lbl.SetLabel(os.path.basename(fileName))
			self.logical = True
		dialog.Destroy()

###############################################################
####################### I/O ###################################
###############################################################
def data_set_flow():
	[datTC,datRA,NG]=[[],[],[]]
	fp = open("flow.inp")

	fp.readline()
	while True:
		line=fp.readline().rstrip()
		if line[:3] =="---": break
		arr=[line]
		for i in range(4): arr.append(fp.readline().rstrip())
		datTC.append(arr)
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
		arr=line.rstrip().split(',')
		arr=map(int,arr)
		arr[0]=[arr[0]]
		arr[1]=arr[1:]
		NG.append(arr[0:2])

	fp.close()

	return [datTC,datRA,NG]

###############################################################
#################### For First Panel Generation ###############
###############################################################
def set_panel(mainframe,base_panel):
	global NG,fn,TC,RA,cb,frame
	frame=mainframe
	RightPanel = wx.Panel(base_panel,wx.ID_ANY)
	RightPanel.SetBackgroundColour('#FFFFFF')
	layoutRight = wx.BoxSizer(wx.VERTICAL)

	#data set
	[datTC,datRA,NG]=data_set_flow()

	# Grid File Name
	fn=FileName(RightPanel)
	layoutRight.Add(fn.panel,flag=wx.GROW,border=10)

	# TextCtrl
	TC=[]
	for i,v in enumerate(datTC):
		TC.append(TCelm(i,RightPanel,datTC[i]))
		layoutRight.Add(TC[i].panel,	flag=wx.GROW,border=10)

	# Radio Buttons
	RA=[]
	for i in range(len(datRA)):
		RA.append(RAelm(i,RightPanel,datRA[i]))
		layoutRight.Add(RA[i],		flag=wx.GROW,border=10)

	cb=CONTROL_BUTTONS_FLOW(RightPanel)
	layoutRight.Add(cb.panel,flag=wx.ALIGN_RIGHT)
	RightPanel.SetSizer(layoutRight)

	#initialize
	for i,v in enumerate(TC):
		v.SetValue(v.initial_value)
		process_values(2,i)
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
				# if [typeID,partID] disables [i,j]
				if RA[i].disabledBy[j][k] == [typeID,partID]:
					del RA[i].disabledBy[j][k]
					# if no selection denies selecting "item j"
					if len(RA[i].disabledBy[j]) == 0 :
						RA[i].EnableItem(j,True)
						# if selected item now is "item j"
						if(RA[i].GetSelection() == j):
							RA[i].isValidate = True

	#set disability caused by new (typeID,partID)
	if   typeID==2:# TC
		value = TC[partID].GetValue()
		if(value.isdigit() and int(value)>=0):
			TC[partID].isValidate = True
			value=int(value)
		else:
			TC[partID].isValidate = False
			return
	elif typeID==3:# RA
		RA[partID].isValidate = True
		value = RA[partID].GetSelection()

	for NGelm in NG:
		KeyValue = NGelm[typeID-2][partID]
		if KeyValue == -1:continue
		if typeID==2 and KeyValue == value:continue
		if typeID==3 and KeyValue != value:continue

		for i,AllowedValue in enumerate(NGelm[0]):
			if (typeID == 2) and (partID == i):continue
			if TC[i].isValidate\
				and (AllowedValue != -1)\
				and (AllowedValue != int(TC[i].GetValue())):
				TC[i].SetValue('')
				TC[i].isValidate = False

		for i,BadValue in enumerate(NGelm[1]):
			if (typeID == 3) and (partID == i):continue
			if BadValue != -1:
				RA[i].EnableItem(BadValue,False)
				RA[i].disabledBy[BadValue].append([typeID,partID])
				if RA[i].GetSelection() == BadValue:
					RA[i].isValidate = False

# general interface
def event_handler(e):
	ID = e.GetId()
	typeID = ID%10;partID = ID/10
	process_values(typeID,partID)

def process_values(typeID,partID):
	if(  typeID == 1): # FileOpen
		fn.OnOpen()
	elif(typeID == 2 or typeID ==3):
		checkvalues(typeID,partID)
	elif(typeID == 9): # close button
		if   partID == 98:
			output_and_run()
		else:
			frame.Close()

	cb.enable_generate_button(fn,TC,RA)

###############################################################
#################### For Execution ############################
###############################################################
def output_and_run():
	os.system("rm -f checkout_chem.inp")
	os.system("rm -f checkout_model.inp")

	# output checkout.inp
	finp = open('checkout.inp','w')
	finp.write(fn.inputstring())
	for var in TC:
		finp.write(var.inputstring())
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

