#!/usr/bin/python
class BD_elm:
   def __init__(self,grid,plane,cind,lu,index):
      self.grid  = grid
      self.plane = plane
      self.cind  = cind		# index i or j
      self.lu    = lu		# lower and upper
      self.index = index	# index from which start
   def show(self):
      print "plane=%3i const=%2s lu=%6s index=%4i coordinate=(%12.5e, %12.5e)" % (self.plane,self.ccind(),self.clu(),self.index,self.grid[0],self.grid[1])
   def ccind(self):
      if self.cind==0:
         return 'i'
      else:
         return 'j'
   def clu(self):
      if self.lu==0:
         return 'lower'
      else:
         return 'upper'
   def isSameLine(self,other):
      return  (self.plane == other.plane) \
           and (self.cind  == other.cind ) \
           and (self.lu    == other.lu   )

class TOUCH_elm:
	def __init__(self,*argv):
		if len(argv)==3:
			[s1,s2,pm]=argv
			self.geo   = s1.cind*2+s2.cind+1
			self.order = pm
			self.width = 1
			self.p1    = s1.plane
			self.p2    = s2.plane
			self.cind1 = s1.cind
			self.cind2 = s2.cind
			self.lu1   = s1.lu
			self.lu2   = s2.lu
			self.s1    = s1.index+1
			if self.order>0:
			   self.s2    = s2.index+1
			else:
			   self.s2    = s2.index
			self.pm1   = -1+2*s1.lu
			self.pm2   = -1+2*s2.lu
		else:
			fromCopy = argv[0]
			self.geo   =  fromCopy.geo  
			self.order =  fromCopy.order
			self.width =  fromCopy.width
			self.p1    =  fromCopy.p1   
			self.p2    =  fromCopy.p2   
			self.cind1 =  fromCopy.cind1
			self.cind2 =  fromCopy.cind2
			self.lu1   =  fromCopy.lu1  
			self.lu2   =  fromCopy.lu2  
			self.s1    =  fromCopy.s1   
			self.s2    =  fromCopy.s2   
			self.pm1   =  fromCopy.pm1  
			self.pm2   =  fromCopy.pm2  

	def fetchc12(self,nijk):
		if self.pm1 <0:
			c1 =1
		else:
			c1 =nijk[self.p1][self.cind1]-1
                                                 
		if self.pm2 <0:                  
			c2 =1                    
		else:                            
			c2 =nijk[self.p2][self.cind2]-1
		return [c1,c2]

	def lineMPI(self,nijk):
		[c1,c2]=self.fetchc12(nijk)
		line1 = "%5i %5i %5i %5i %5i %5i %5i %5i\n" \
		            % (self.cind1,self.width,1,\
				self.p1+1,c1,self.s1,self.pm1,\
				self.myid2)
		line2 = "%5i %5i %5i %5i %5i %5i %5i %5i\n" \
		            % (self.cind2,self.width,self.order,\
				self.p2+1,c2,self.s2,self.pm2,\
				self.myid1)
		return [line1,line2]

	def lineCoPro(self,nijk):
		[c1,c2]=self.fetchc12(nijk)
		return"%5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n" \
	            % (self.geo,self.width,self.order,\
			self.p1+1,c1,self.s1,self.pm1,\
			self.p2+1,c2,self.s2,self.pm2)
	
	def line(self,nijk):
		[c1,c2]=self.fetchc12(nijk)

		if 'myid1' in self.__dict__.keys():
			return"%5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n" \
	        	    % (self.geo,self.width,self.order,\
				self.p1+1,c1,self.s1,self.pm1,\
				self.p2+1,c2,self.s2,self.pm2,self.myid1,self.myid2)
		elif 'myid' in self.__dict__.keys():
			return"%5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n" \
	        	    % (self.geo,self.width,self.order,\
				self.p1+1,c1,self.s1,self.pm1,\
				self.p2+1,c2,self.s2,self.pm2,self.myid)
		else:
			return"%5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n" \
	        	    % (self.geo,self.width,self.order,\
				self.p1+1,c1,self.s1,self.pm1,\
				self.p2+1,c2,self.s2,self.pm2)
	def show(self,nijk):
		if 'myid1' in self.__dict__.keys():
			print "# geo width order    p1    c1    s1   pm1    p2    c2    s2   pm2 myid1 myid2"
			print self.line(nijk),
		elif 'myid' in self.__dict__.keys():
			print "# geo width order    p1    c1    s1   pm1    p2    c2    s2   pm2  myid"
			print self.line(nijk),
		else:
			print "# geo width order    p1    c1    s1   pm1    p2    c2    s2   pm2"
			print self.line(nijk),
