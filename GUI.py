#!/usr/bin/python
import wx
from subprocess import Popen,PIPE

import leftpanel
import gui_flow

def output_and_run():
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

def init_values():
	global TC,RA
	for i,v in enumerate(TC):
		v.SetValue(v.initial_value)
		process_values(2,i)
	for i,v in enumerate(RA):
		v.SetSelection(0)
		process_values(3,i)


def enablevalues(typeID,partID):
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

def checkvalues(typeID,partID):
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

def event_handler(e):
	ID = e.GetId()
	typeID = ID%10;partID = ID/10
	process_values(typeID,partID)

def process_values(typeID,partID):
	if(  typeID == 1): # FileOpen
		fn.OnOpen(frame)
	elif(typeID == 2 or typeID ==3): # TextCtrl
		enablevalues(typeID,partID)
		checkvalues(typeID,partID)
	elif(typeID == 9): # close button
		if   partID == 98:
			output_and_run()
		frame.Close()

	cb.enable_generate_button(fn,TC,RA)

if __name__ == "__main__":
	app = wx.App()
	frame = wx.Frame(None,wx.ID_ANY,"Shimada Lab. CODE")
	frame.SetBackgroundColour('#FFFFFF')
	layout = wx.BoxSizer(wx.HORIZONTAL)

	LeftPanel = leftpanel.leftpanel(frame)
	layout.Add(LeftPanel)

	[RightPanel,NG,fn,TC,RA,cb]=gui_flow.set_flowpanel(frame,event_handler)
	init_values()
	layout.Add(RightPanel,flag=wx.ALL,border=10)

	frame.SetSizer(layout)
	frame.SetAutoLayout(True)
	layout.Fit(frame)
	frame.Show()
	app.MainLoop()
