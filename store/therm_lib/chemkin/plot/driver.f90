program driver
   use const_chem
   use conditions
   implicit none

   !read files
   call init_therm
   call read_conditions
   call init_vrho
   call set_Tref
   call out_Tref
   call out_plt
end program driver
