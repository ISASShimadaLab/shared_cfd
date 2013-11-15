#!/usr/bin/python
import copy
import os

#### set constants ##############################
ngx=17
ngy=7

#### sort data ##############################
os.system("sort balance.dat -uo balance.dat")

##### read all of the data #######################
#file open
fi = open("balance.dat") 
fo = open("trimed_balance.dat","w")

for i in range(1,ngx+1):
   #read iso-x grid line
   buf = (ngy+1)*[0]
   for j in range(1,ngy+1):
      line  = fi.readline()
      array = line.split()
      array[1] = int(array[1])
      array[2] = int(array[2])
      buf[array[2]] = array
   
   for j in range(1,ngy+1):
      fo.write("%s %3.3i %3.3i %s %s\n" % (buf[j][0],buf[j][1]-1,buf[j][2]-1,buf[j][3],buf[j][4]))
      fo.write("%s %3.3i %3.3i %s %s\n" % (buf[j][0],buf[j][1]-1,buf[j][2]  ,buf[j][3],buf[j][4]))
   fo.write("\n")                                                           
                                                                            
   for j in range(1,ngy+1):                                                 
      fo.write("%s %3.3i %3.3i %s %s\n" % (buf[j][0],buf[j][1]  ,buf[j][2]-1,buf[j][3],buf[j][4]))
      fo.write("%s %3.3i %3.3i %s %s\n" % (buf[j][0],buf[j][1]  ,buf[j][2]  ,buf[j][3],buf[j][4]))
   fo.write("\n")


#file close
fi.close()
fo.close()
