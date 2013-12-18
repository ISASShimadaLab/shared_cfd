program driver
   use conditions
   implicit none
   logical flag

   !read files
   call init_pack_cea('hp')
   call read_conditions
   call calcTeq(flag)
   if(.not.flag) print *,"Calculation failed to converge."
end program driver
