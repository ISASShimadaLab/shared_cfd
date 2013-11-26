#!/usr/bin/python
import wx
import os
from PIL import Image
from subprocess import Popen,PIPE

def data_set():
	global NG,datRA,datTC,datOption

	[datTC,datRA,datOption]=[[],[],[]]

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
	datRA.append(["thermal model",("Ideal Gas",\
			"Cold Flow by chemical kinetics D.B.",\
			"Chemical Kinetics Model"),"thermal model     : ",\
			"  #Ideal Gas(0)/Cold flow using chemical kinetics database(1) / Chemical Kinetics Model(2)"])

	datOption.append([["specific heat ratio","1.4"],\
			["gas constant (J/kg/K)","287"],\
			["coefficient of kinematic viscosity (m^2/s)","1.6e-5"]])
	datOption.append([])
	datOption.append([])

	NG.append([[-1],[-1, 0, 4,-1,-1,-1]])
	NG.append([[-1],[-1, 1, 5,-1,-1,-1]])
	NG.append([[-1],[ 1,-1, 5,-1,-1,-1]])
	NG.append([[-1],[-1, 0, 6,-1,-1,-1]])
	NG.append([[-1],[-1,-1,-1,-1,-1, 0]])
	NG.append([[-1],[-1, 1,-1,-1,-1, 2]])

class OptionPanel:
	inphead = "ParametersForAbove: "
	inptail = " #depends on the model. see text."
	isValidate = True
	panel = None
	Num=0
	NumNow=0
	TC=[]
	ST=[]
	def __init__(self,datOption,RightPanel):
		for var in datOption:
			self.Num = max(len(var),self.Num)
		
		self.panel = wx.Panel(RightPanel,wx.ID_ANY)
		self.panel.SetBackgroundColour('#FFFFFF')
		layout = wx.StaticBoxSizer(\
				wx.StaticBox(\
					self.panel,\
					wx.ID_ANY,\
					"parameters for thermal model"),\
				wx.VERTICAL)
		for i in range(self.Num):
			tmppanel = wx.Panel(self.panel,wx.ID_ANY)
			tmppanel.SetBackgroundColour('#FFFFFF')
			tmplayout = wx.BoxSizer(wx.HORIZONTAL)
			self.TC.append(wx.TextCtrl(  tmppanel,i*10+4))
			self.ST.append(wx.StaticText(tmppanel,wx.ID_ANY))
			tmplayout.Add(self.TC[i],flag=wx.ALL|wx.ALIGN_CENTER,border=5)
			tmplayout.Add(self.ST[i],flag=wx.ALL|wx.ALIGN_CENTER|wx.GROW,border=5)
			self.TC[i].Bind(wx.EVT_TEXT,event_handler)
			tmppanel.SetSizer(tmplayout)
			layout.Add(tmppanel,flag=wx.GROW)
		self.panel.SetSizer(layout)

	def change_plane(self,i):
		self.NumNow = i
		for i,var in enumerate(datOption[self.NumNow]):
			self.ST[i].SetLabel(":\t"+var[0])
			self.TC[i].SetValue(      var[1])
			self.TC[i].Enable()
		for i in range(len(datOption[self.NumNow]),self.Num):
			self.ST[i].SetLabel("")
			self.TC[i].SetValue("")
			self.TC[i].Disable()

	def check_validation(self):
		self.isValidate = True
		for i,var in enumerate(datOption[self.NumNow]):
			if(self.TC[i].GetValue() == ''):
				self.isValidate = False

	def inputstring(self):
		string =""
		string += self.inphead
		for i,var in enumerate(datOption[self.NumNow]):
			string += self.TC[i].GetValue() + " "
		string += self.inptail
		return string

def output_and_run():
	finp = open('checkout.inp','w')
	finp.write(LblFileName.inphead+LblFileName.fullpath+LblFileName.inptail+"\n")
	for var in TC:
		finp.write(var.inphead+var.GetValue()+var.inptail+"\n")
	for var in RA:
		finp.write(var.inphead+str(var.GetSelection())+var.inptail+"\n")
	finp.write(option.inputstring())
	finp.close()

#	os.system("./checkout.py")
#	wx.MessageBox("Output from checkout.py","Message")

	p = Popen("./checkout.py",stdout = PIPE)
	printedout = p.stdout.read()
	value = p.wait()
	wx.MessageBox("Output from checkout.py\n**************************\n"+printedout+"**************************\nexit value: "+str(value),"Message")

	frame.Close()

def init_values():
	global TC,RA,datTC,datRA
	for i in range(len(TC)):
		TC[i].SetValue(datTC[i][2])
		process_values(2,i)
	for i in range(len(RA)):
		RA[i].SetSelection(0)
		process_values(3,i)

def enablevalues(typeID,partID):
	arrID = [typeID,partID]
	for i in range(len(RA)):
		for j in range(RA[i].GetCount()):
			for k in range(len(RA[i].disabledBy[j])-1,-1,-1):
				if RA[i].disabledBy[j][k] == arrID:
					del RA[i].disabledBy[j][k]
					if len(RA[i].disabledBy[j]) == 0 :
						RA[i].EnableItem(j,True)
						if(RA[i].GetSelection() == j):
							RA[i].isValidate = True

def checkvalues(typeID,partID,arr):
	for i in range(len(arr[0])):
		if (typeID == 2) and (partID == i):
			continue

		BadValue = arr[0][i]
		if TC[i].isValidate and (BadValue != -1) and (BadValue != int(TC[i].GetValue())):
			TC[i].SetValue('')
			TC[i].isValidate = False

	for i in range(len(arr[1])):
		if (typeID == 3) and (partID == i):
			continue

		if arr[1][i] >= 0:
			RA[i].EnableItem(arr[1][i],False)
			RA[i].disabledBy[arr[1][i]].append([typeID,partID])
			if RA[i].GetSelection() == arr[1][i]:
				RA[i].isValidate = False

def event_handler(e):
	# process actions
	ID = e.GetId()
	typeID = ID%10
	partID = ID/10
	if(typeID == 9):   # close button
		if   partID == 98:
			output_and_run()
		elif partID == 99:
			frame.Close()
	else:
		process_values(typeID,partID)

def process_values(typeID,partID):
	if(  typeID == 1): # FileOpen
		OnOpen()
	elif(typeID == 2): # TextCtrl
		enablevalues(typeID,partID)
		value = TC[partID].GetValue()
		if(value.isdigit() and int(value)>=0):
			TC[partID].isValidate = True
			for arr in NG:
				BadValue = arr[0][partID]
				if (BadValue != -1) and (BadValue != int(value)):
					checkvalues(typeID,partID,arr)
		else:
			TC[partID].isValidate = False
	elif(typeID == 3): # radiobox
		enablevalues(typeID,partID)
		RA[partID].isValidate = True
		for arr in NG:
			if arr[1][partID] == RA[partID].GetSelection():
				checkvalues(typeID,partID,arr)
		if(partID == 5):
			option.change_plane(RA[partID].GetSelection())
	elif(typeID == 4): # option panel
		option.check_validation()

	#enable generate button
	lgcsum = BtnFileName.logical
	for var in TC:
		lgcsum &= var.isValidate
	for var in RA:
		lgcsum &= var.isValidate
	lgcsum &= option.isValidate

	if(lgcsum):
		button_generate.Enable()
	else:
		button_generate.Disable()

def OnOpen():
	dirName = os.getcwd()
	dialog = wx.FileDialog(frame, "Choose the grid file", dirName, "", "*.x", wx.OPEN)
	if dialog.ShowModal() == wx.ID_OK:
		fileName = dialog.GetPath()
		LblFileName.fullpath = fileName
		LblFileName.SetLabel(os.path.basename(fileName))
		BtnFileName.logical = True
	dialog.Destroy()


if __name__ == "__main__":
	[RA,TC,NG] = [[],[],[]]
	height = 550

	app = wx.App()
	frame = wx.Frame(None,wx.ID_ANY,"Shimada Lab. CODE")
	frame.SetBackgroundColour('#FFFFFF')
	layout = wx.BoxSizer(wx.HORIZONTAL)

	################# for LeftPanel ############################
	LeftPanel = wx.Panel(frame,wx.ID_ANY)
	LeftPanel.SetBackgroundColour('#FFFFFF')
	layoutLeft = wx.BoxSizer(wx.VERTICAL)

	# for ImagePanel
	ImagePanel = wx.Panel(LeftPanel,wx.ID_ANY)
	pil = Image.open('store/gui/top.jpg')
	ratio = float(height-10*2)/float(pil.size[1])
	new_size = (int(pil.size[0]*ratio),int(pil.size[1]*ratio))
	pil.thumbnail(new_size,Image.ANTIALIAS)
	image = wx.EmptyImage(pil.size[0],pil.size[1])
	image.SetData(pil.convert('RGB').tostring())
	wx.StaticBitmap(ImagePanel, wx.ID_ANY, image.ConvertToBitmap())
	layoutLeft.Add(ImagePanel,flag=wx.ALL,border=10)

	# for LogoPanel
	LogoPanel = wx.Panel(LeftPanel,wx.ID_ANY)
	pil = Image.open('store/gui/logo.png')
	ratio = float(image.GetWidth()-10*2)/float(pil.size[0])
	new_size = (int(pil.size[0]*ratio),int(pil.size[1]*ratio))
	pil.thumbnail(new_size,Image.ANTIALIAS)
	image = wx.EmptyImage(pil.size[0],pil.size[1])
	image.SetData(pil.convert('RGB').tostring())
	wx.StaticBitmap(LogoPanel, wx.ID_ANY, image.ConvertToBitmap())
	LogoPanel.SetSize(new_size)
	layoutLeft.Add(LogoPanel,flag=wx.ALL,border=10)

	################# for RightPanel ############################
	RightPanel = wx.Panel(frame,wx.ID_ANY)
	RightPanel.SetBackgroundColour('#FFFFFF')
	layoutRight = wx.BoxSizer(wx.VERTICAL)

	#data set
	data_set()

	# Grid File Name
	tmpPanel = wx.Panel(RightPanel,wx.ID_ANY)
	tmpPanel.SetBackgroundColour('#FFFFFF')
	tmpbox    = wx.StaticBox(tmpPanel,wx.ID_ANY,"Grid File Name")
	tmplayout = wx.StaticBoxSizer(tmpbox,wx.HORIZONTAL)
	BtnFileName = wx.Button(tmpPanel,1,"Choose File.")
	BtnFileName.Bind(wx.EVT_BUTTON,event_handler)
	tmplayout.Add(BtnFileName,flag=wx.ALL|wx.ALIGN_CENTER,border=5)
	LblFileName =wx.StaticText(tmpPanel,wx.ID_ANY,"No File Selected yet.") 
	LblFileName.inphead = "grid file name    : "
	LblFileName.inptail = " #Grid file name (*.x)"
	tmplayout.Add(LblFileName,flag=wx.ALL|wx.ALIGN_CENTER,border=5)
	tmpPanel.SetSizer(tmplayout)
	layoutRight.Add(tmpPanel,flag=wx.GROW,border=10)
	BtnFileName.logical = False

	# TextCtrl
	for i in range(len(datTC)):
		tmpPanel = wx.Panel(RightPanel,wx.ID_ANY)
		tmpPanel.SetBackgroundColour('#FFFFFF')
		tmpbox    = wx.StaticBox(tmpPanel,wx.ID_ANY,datTC[i][0])
		tmplayout = wx.StaticBoxSizer(tmpbox,wx.HORIZONTAL)
		TC.append(wx.TextCtrl(tmpPanel,i*10+2,datTC[i][2]))
		TC[i].Bind(wx.EVT_TEXT,event_handler)
		TC[i].isValidate = True
		TC[i].inphead = datTC[i][3]
		TC[i].inptail = datTC[i][4]
		tmplayout.Add(TC[i],flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		tmplayout.Add(wx.StaticText(tmpPanel,wx.ID_ANY,datTC[i][1]),\
			flag=wx.ALL|wx.ALIGN_CENTER,border=5)
		tmpPanel.SetSizer(tmplayout)
		layoutRight.Add(tmpPanel,flag=wx.GROW,border=10)

	# Radio Buttons
	for i in range(len(datRA)):
		RA.append(wx.RadioBox(RightPanel,i*10+3,datRA[i][0],\
				choices=datRA[i][1],style = wx.RA_HORIZONTAL,\
				majorDimension=3))
		RA[i].Bind(wx.EVT_RADIOBOX,event_handler)
		RA[i].isValidate = True
		RA[i].inphead = datRA[i][2]
		RA[i].inptail = datRA[i][3]
		RA[i].disabledBy = [[] for ii in range(RA[i].GetCount())]
		layoutRight.Add(RA[i],flag=wx.GROW,border=10)

	# Options Plane
	option = OptionPanel(datOption,RightPanel)
	layoutRight.Add(option.panel,flag=wx.GROW,border=10)


	# close and generate Buttons
	tmpPanel = wx.Panel(RightPanel,wx.ID_ANY)
	tmpPanel.SetBackgroundColour('#FFFFFF')
	tmplayout = wx.BoxSizer(wx.HORIZONTAL)
	button_generate = wx.Button(tmpPanel,989,"Generate CODE!")
	button_close    = wx.Button(tmpPanel,999,"Close")
	button_generate.Bind(wx.EVT_BUTTON,event_handler)
	button_close.Bind(wx.EVT_BUTTON,event_handler)
	tmplayout.Add(button_generate,flag=wx.ALL|wx.ALIGN_RIGHT,border=5)
	tmplayout.Add(button_close   ,flag=wx.ALL|wx.ALIGN_RIGHT,border=5)
	tmpPanel.SetSizer(tmplayout)
	layoutRight.Add(tmpPanel,flag=wx.ALIGN_RIGHT)
	button_generate.Disable()

	# initialize input values
	init_values()
	
	################# for ALL ###################################
	LeftPanel.SetSizer(layoutLeft)
	RightPanel.SetSizer(layoutRight)
	layout.Add(LeftPanel)
	layout.Add(RightPanel,flag=wx.ALL,border=10)
	
	frame.SetSizer(layout)
	frame.SetAutoLayout(True)
	layout.Fit(frame)
	
	frame.Show()
	app.MainLoop()
