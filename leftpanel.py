#!/usr/bin/pytyoh
import wx
from PIL import Image
def leftpanel(frame):
	height = 530 #image height

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

	LeftPanel.SetSizer(layoutLeft)
	return LeftPanel
