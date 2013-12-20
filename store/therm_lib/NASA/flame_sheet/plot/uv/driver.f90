program driver
   use conditions
   implicit none

   !read files
   call init_pack_flame_sheet
   call read_conditions
   call calcTeq
   call out_plt
   print *,"done"
end program driver
