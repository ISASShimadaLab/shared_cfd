#!/usr/bin/python
import wx
import os
from control_buttons import *
###############################################################
############################ CONTROL_BUTTONS_FLOW #############
###############################################################
class CONTROL_BUTTONS_FLOW(CONTROL_BUTTONS):
	def __init__(self,RightPanel,event_handler):
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
	def __init__(self,i,RightPanel,datTC,event_handler):
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
	def __init__(self,i,RightPanel,datRA,event_handler):
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
############################ FN ###############################
###############################################################
class FileName:
	def __init__(self,RightPanel,event_handler):
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

	def OnOpen(self,frame):
		dirName = os.getcwd()
		dialog = wx.FileDialog(frame, "Choose the grid file", dirName, "", "*.x", wx.OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			fileName = dialog.GetPath()
			self.fullpath = fileName
			self.lbl.SetLabel(os.path.basename(fileName))
			self.logical = True
		dialog.Destroy()


###############################################################
############################ UTILITY ##########################
###############################################################
def data_set_flow():
	[datTC,datRA,datOption,NG]=[[],[],[],[]]

	datTC.append([	"Maximum Number of Processors",\
			"0 @ Personal Computer",\
			"0",\
			"NumOfProcessors   : ",\
			"  #The Maximum Number of Processors @ M-System[Super Computer] / 0 @ Personal Computer"])

	datRA.append(["dimension",("two dimension","axi-symmetry"),"dimension         : ","  #two dimension(0) / axi-symmetry(1)"])
	datRA.append(["time scale",("global","local"),"time step         : ","  #global(0)/local(1)"])
	datRA.append(["time scheme",("Euler","two-step Runge-Kutta",\
			"LU-SGS","Newton-Raphson",\
			"Pre-Conditioner","Dual Time Step","Pre-conditioner with LU-SGS"),"time scheme       : ",\
			"  #Euler(0)/Runge-Kutta two-step(1)/LU-SGS(2)/Newton-Raphson(3)/Pre-Conditioner(4)/Dual Time Step(5)/Pre-conditioner with LU-SGS(6)"])
	datRA.append(["high order scheme",("Up Wind","MUSCL"),"high order scheme : ","  #Up Wind(0)/MUSCL(1)"])
	datRA.append(["viscosity",("non-viscous flow","viscous flow"),"viscosity         : ",\
			"  #non-viscous flow(0)/viscous flow(1)"])
	datRA.append(["thermal model",(
				"Ideal Gas",\
				"Cold Flow by chemical kinetics D.B.",\
				"Chemical Kinetics Model",\
				"Flame Sheet Model"),\
			"thermal model     : ",\
			"  #Ideal Gas(0)/Cold flow using chemical kinetics database(1) / Chemical Kinetics Model(2) / Flame Sheet Model(3)"])

	datOption.append([["specific heat ratio","1.4"],\
			["gas constant (J/kg/K)","287"],\
			["coefficient of kinematic viscosity (m^2/s)","1.6e-5"]])
	datOption.append([])
	datOption.append([])
	datOption.append([])

	NG.append([[-1],[-1, 0, 4,-1,-1,-1]])
	NG.append([[-1],[-1, 1, 5,-1,-1,-1]])
	NG.append([[-1],[ 1,-1, 5,-1,-1,-1]])
	NG.append([[-1],[-1, 0, 6,-1,-1,-1]])
	NG.append([[-1],[-1,-1,-1,-1,-1, 0]])
	NG.append([[-1],[-1, 1,-1,-1,-1, 2]])

	return [datTC,datRA,datOption,NG]


def set_flowpanel(frame,event_handler):
	RightPanel = wx.Panel(frame,wx.ID_ANY)
	RightPanel.SetBackgroundColour('#FFFFFF')
	layoutRight = wx.BoxSizer(wx.VERTICAL)

	#data set
	[datTC,datRA,datOption,NG]=data_set_flow()

	# Grid File Name
	fn=FileName(RightPanel,event_handler)
	layoutRight.Add(fn.panel,flag=wx.GROW,border=10)

	# TextCtrl
	TC=[]
	for i,v in enumerate(datTC):
		TC.append(TCelm(i,RightPanel,datTC[i],event_handler))
		layoutRight.Add(TC[i].panel,	flag=wx.GROW,border=10)

	# Radio Buttons
	RA=[]
	for i in range(len(datRA)):
		RA.append(RAelm(i,RightPanel,datRA[i],event_handler))
		layoutRight.Add(RA[i],		flag=wx.GROW,border=10)

	cb=CONTROL_BUTTONS_FLOW(RightPanel,event_handler)
	layoutRight.Add(cb.panel,flag=wx.ALIGN_RIGHT)
	RightPanel.SetSizer(layoutRight)

	return [RightPanel,NG,fn,TC,RA,cb]























