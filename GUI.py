#!/usr/bin/python
import os
import wx
import store.gui.leftpanel
import store.gui.notepanel
if __name__ == "__main__":
	app = wx.App()
	frame = wx.Frame(None,wx.ID_ANY,"Shimada Lab. CODE")
	frame.SetBackgroundColour('#FFFFFF')
	layout = wx.BoxSizer(wx.HORIZONTAL)

	#leftpanel
	layout.Add(store.gui.leftpanel.leftpanel(frame))

	#rightpanel
	layout.Add(store.gui.notepanel.notepanel(frame),flag=wx.ALL,border=5)

	# show panel
	frame.SetSizer(layout)
	frame.SetAutoLayout(True)
	layout.Fit(frame)
	frame.Show()
	app.MainLoop()

	os.system("rm -f store/*.pyc store/checkout/*.pyc store/gui/*.pyc")
