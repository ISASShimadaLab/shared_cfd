#!/usr/bin/python
import wx
import leftpanel
import notepanel
if __name__ == "__main__":
	app = wx.App()
	frame = wx.Frame(None,wx.ID_ANY,"Shimada Lab. CODE")
	frame.SetBackgroundColour('#FFFFFF')
	layout = wx.BoxSizer(wx.HORIZONTAL)

	#leftpanel
	layout.Add(leftpanel.leftpanel(frame))

	#rightpanel
	layout.Add(notepanel.notepanel(frame),flag=wx.ALL,border=5)

	# show panel
	frame.SetSizer(layout)
	frame.SetAutoLayout(True)
	layout.Fit(frame)
	frame.Show()
	app.MainLoop()
