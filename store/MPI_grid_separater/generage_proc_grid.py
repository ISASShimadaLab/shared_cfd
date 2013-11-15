#!/usr/bin/python
def makelist(arr,ng):
   sum_r=0
   sum_a=0
   for i in range(len(arr)):
      if arr[i]<0:
         sum_r -= arr[i]
      else:
         sum_a += arr[i]
   
   unit = (ng-sum_a)/sum_r
   res  = (ng-sum_a)%sum_r
   
   now=1
   line=""
   for i in range(len(arr)):
      line+=("%3.3i " % now)
      if (i+1)%5 == 0:
         line+="   "
   
      if arr[i]<0:
         delta = -arr[i]*unit
      else:
         delta =  arr[i]
   
      if i<res:
         delta+=1
      now+=delta
   
   return line

if __name__ == "__main__":
   # set ratio
#   bw0=-2.1
#   bw1=-2.2
#   arr_x=[bw0,bw0,bw0,bw0,bw0,\
#          -1,-1,-1,-1,-1,\
#          -1,-1,-1,-1,bw1,\
#          bw1,bw1]
#   bw1=16
#   arr_y=[bw1,bw1,bw1,bw1,bw1,\
#          -1,-1]

   arr_x=6*[-1]
   arr_y=3*[-1]
   print 6*3
   
   #parameter
   ni = 99
   nj = 49
   
   
   f=open("grid_separation.dat","w")
   line = makelist(arr_x,ni)
   f.write(line+"\n")
   line = makelist(arr_y,nj)
   f.write(line+"\n")
   f.close()

   import os
   os.system("./plot_grid.py")
