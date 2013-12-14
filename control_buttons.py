#!/usr/bin/python
import wx
class CONTROL_BUTTONS:
	def __init__(self,RightPanel,event_handler):
		tmpPanel = wx.Panel(RightPanel,wx.ID_ANY)
		tmpPanel.SetBackgroundColour('#FFFFFF')
		tmplayout = wx.BoxSizer(wx.HORIZONTAL)
		button_generate = wx.Button(tmpPanel,989,"Generate CODE!")
		button_close    = wx.Button(tmpPanel,999,"Close")
		button_generate.Bind(wx.EVT_BUTTON,event_handler)
		button_close   .Bind(wx.EVT_BUTTON,event_handler)
		tmplayout.Add(button_generate,flag=wx.ALL|wx.ALIGN_RIGHT,border=5)
		tmplayout.Add(button_close   ,flag=wx.ALL|wx.ALIGN_RIGHT,border=5)
		tmpPanel.SetSizer(tmplayout)
		button_generate.Disable()
		self.panel=tmpPanel
		self.button_generate=button_generate
