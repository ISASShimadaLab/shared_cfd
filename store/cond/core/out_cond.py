#!/usr/bin/python
import os
from mod_cond import *

def out_cond(toSetManually,MaxMPIcomm,ABOUTNV):
	os.system("cat checkout/condition.head.f90 >  checkout/condition.raw.f90")
	fp = open("checkout/condition.raw.f90",'a')
	
	########## Declaration
	fp.write("   integer,parameter::MaxMPIcomm = %i\n\n" % MaxMPIcomm)

	########## Cut Lines
	fp.write("   call section_exchange\n")
	fp.write("   call MPI_COMMUNICATIONS_CUT_LINE\n\n")

	########## I-DIRECTION
	fp.write("""
   !boundary right and left
   do plane=nps,npe
      select case(plane)""")
	for i,tsmp in enumerate(toSetManually):
	   plane = i+1
	   tsmp=tsmp.pop(0)
	   fp.write("""
      case(%3i)""" % plane)
	
	   # i=1/2
	   b_list=tsmp.pop(0)
	   if len(b_list)!=0:
	      fp.write("""
         if(gx(plane) .eq. 1) then
            !i=1/2""")
	      for arr in b_list:
	         fp.write("""
            do j=max(%5i,nys(plane)),min(nye(plane),%5i)
               !w(:,   0,j,plane) =
               !vhi(:, 0,j,plane) =""" % tuple(arr))
		 if ABOUTNV == "with-nV": fp.write("""
               !Yv(:,  0,j,plane) =""")
		 fp.write("""
               !w(:,  -1,j,plane) =
               !vhi(:,-1,j,plane) =""")
		 if ABOUTNV == "with-nV": fp.write("""
               !Yv(:, -1,j,plane) =""")
	         fp.write("""
            end do
         end if\n""")
	   # i=ni/2
	   b_list=tsmp.pop(0)
	   if len(b_list)!=0:
	      fp.write("""
         if(gx(plane) .eq. ngx(plane)) then
            !i=ni+1/2""")
	      for arr in b_list:
	         fp.write("""
            do j=max(%5i,nys(plane)),min(nye(plane),%5i)
               !w(:,  ni(plane)+1,j,plane) =
               !vhi(:,ni(plane)+1,j,plane) =""" % tuple(arr))
		 if ABOUTNV == "with-nV": fp.write("""
               !Yv(:, ni(plane)+1,j,plane) =""")
		 fp.write("""
               !w(:,  ni(plane)+2,j,plane) =
               !vhi(:,ni(plane)+2,j,plane) =""")
		 if ABOUTNV == "with-nV": fp.write("""
               !Yv(:, ni(plane)+2,j,plane) =""")
	         fp.write("""
            end do
         end if""")
	fp.write("""
      end select
   end do

   call MPI_COMMUNICATIONS_I_DIRECTION
""")
	########## END OF I-DIRECTION
	
	########## J-DIRECTION
	fp.write("""
   !boundary upper and lower
   do plane=nps,npe
      select case(plane)""")
	for i,tsmp in enumerate(toSetManually):
	   plane = i+1
	   tsmp=tsmp.pop(0)
	   fp.write("""
      case(%3i)""" % plane)
	
	   # i=1/2
	   b_list=tsmp.pop(0)
	   if len(b_list)!=0:
	      fp.write("""
         if(gy(plane) .eq. 1) then
            !j=1/2""")
	      for arr in b_list:
	         fp.write("""
            do i=max(%5i,nxs(plane)),min(nxe(plane),%5i)
               !w(:,  i, 0,plane) =
               !vhi(:,i, 0,plane) =""" % tuple(arr))
		 if ABOUTNV == "with-nV": fp.write("""
               !Yv(:, i, 0,plane) =""")
	         fp.write("""
               !w(:,  i,-1,plane) =
               !vhi(:,i,-1,plane) =""")
		 if ABOUTNV == "with-nV": fp.write("""
               !Yv(:, i,-1,plane) =""")
	         fp.write("""
            end do
         end if\n""")
	   # i=ni/2
	   b_list=tsmp.pop(0)
	   if len(b_list)!=0:
	      fp.write("""
         if(gy(plane) .eq.  ngy(plane)) then
            !j=nj+1/2""")
	      for arr in b_list:
	         fp.write("""
            do i=max(%5i,nxs(plane)),min(nxe(plane),%5i)
               !w(:,  i,nj(plane)+1,plane) =
               !vhi(:,i,nj(plane)+1,plane) =""" % tuple(arr))
		 if ABOUTNV == "with-nV": fp.write("""
               !Yv(:, i,nj(plane)+1,plane) =""")
	         fp.write("""
               !w(:,  i,nj(plane)+2,plane) =
               !vhi(:,i,nj(plane)+2,plane) =""")
		 if ABOUTNV == "with-nV": fp.write("""
               !Yv(:, i,nj(plane)+2,plane) =""")
	         fp.write("""
            end do
         end if""")
	fp.write("""
      end select
   end do

   call MPI_COMMUNICATIONS_J_DIRECTION

   call set_corners
""")
	########## END OF J-DIRECTION
	fp.close()
	os.system("cat checkout/condition.tail.f90 >> checkout/condition.raw.f90")


def out_cut(touch,nijk,Nproc):
	file_cut = Nproc*[""]
	for var in touch:
		file_cut[var.myid]+=var.lineCoPro(nijk)

	for i,line in enumerate(file_cut):
		if line != "":
			fo = open("checkout/cut_copro%03i.inp" % i,'w')
			fo.write("# myid = %03i\n" % i)
			fo.write("%i\n" % line.count('\n'))
			fo.write("# geo width order    p1    c1    s1   pm1    p2    c2    s2   pm2\n")
			fo.write(line)
			fo.close()

def out_MPI(MPIcomm,nijk,Nproc):
	MaxMPIcomm = 0
	file_MPI = Nproc*[""]
	for var in MPIcomm:
		[line1,line2] = var.lineMPI(nijk)
		file_MPI[var.myid1] +=line1
		file_MPI[var.myid2] +=line2
	for i,line in enumerate(file_MPI):
		if line != "":
			fp = open("checkout/MPIcomm%03i.inp" % i,"w")
			fp.write("# myid = %03i\n" % i)
			fp.write("%i\n" % line.count('\n'))
			fp.write("#cind width order     p     c     s    pm 2send\n")
			fp.write(line)
			fp.close()
			MaxMPIcomm = max(MaxMPIcomm, line.count('\n'))
	return MaxMPIcomm
