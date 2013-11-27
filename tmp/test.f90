program main
   implicit none
   character*200 line
   character*200 bufnext
   bufnext=""
   open(55,file="chem.inp")

   do
      line=next_word()
      print '(a)',trim(line)
   end do

   close(55)
contains
character*200 function next_word()
   implicit none
   integer ind
   if(bufnext .eq. "") bufnext=next_line()
   ind = index(bufnext,' ')
   if(ind == 0) then
      next_word=bufnext
      bufnext=""
   else
      next_word=bufnext(:ind-1)
      bufnext=trim(adjustl(bufnext(ind:)))
   end if
end function next_word
character*200 function next_line()
   implicit none
   integer ind
   do while(bufnext .eq. "")
      read(55,'(a)') bufnext
      bufnext=trim(adjustl(bufnext))
      ind = index(bufnext,'!')
      if(ind == 1) then
         bufnext=""
      else if(ind > 1) then
         bufnext=bufnext(:ind-1)
      end if
   end do
   next_line=bufnext
   bufnext=""
end function next_line
end program main

