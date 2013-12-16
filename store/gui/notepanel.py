#!/usr/bin/python
import wx
import gui_flow
import gui_model
import gui_chem

class PageOne(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
	self.SetBackgroundColour('#FFFFFF')
        t = wx.StaticText(self, -1, "This is a PageOne object", (20,20))
class notepanel(wx.Panel):
	def __init__(self,frame):
		wx.Panel.__init__(self,frame)
		self.SetBackgroundColour('#FFFFFF')

		notebook = wx.Notebook(self,wx.ID_ANY,style=wx.NB_TOP)
		notebook.SetBackgroundColour('#FFFFFF')

		notebook.InsertPage(0,gui_chem. set_panel(frame,notebook),"chemical calculation")
		notebook.InsertPage(1,gui_flow. set_panel(frame,notebook),"flow")
		notebook.InsertPage(2,gui_model.set_panel(frame,notebook),"chemical model")

        	sizer = wx.BoxSizer()
        	sizer.Add(notebook, 2, wx.EXPAND)
        	self.SetSizer(sizer)
