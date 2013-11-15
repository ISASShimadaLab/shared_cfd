program main
   use grbl_prmtr
   use prmtr
   implicit none
   double precision,dimension((dimq+3)*ni*nj)::buf
   integer step,step_new
   integer i,j,k
   character*50 filename
   character*50 str_buf
   character*50,parameter::outfile="out_change.bin"

   if(iargc() .ne. 2) then
      print *,"Wrong command-line arguments."
      call exit(0)
   end if

   call getarg(1,filename)
   call getarg(2,str_buf)
   read(str_buf,*) step_new
   print *,filename,step_new

   print '(a,a)',"inputfile :", trim(filename)
   print '(a,a)',"outputfile:", trim(outfile)

   !read q binary
   open(55,file=trim(filename),form="unformatted")
   open(66,file=trim(outfile), form="unformatted")

   read(55) step
   print '(i12,a,i12)',step,"->",step_new
   write(66) step_new

   read (55) buf
   write(66) buf

   close(55)
   close(66)
   print '(a)',"done."
end program main
