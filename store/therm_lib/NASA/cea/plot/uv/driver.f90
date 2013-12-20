program driver
   use conditions
   implicit none

   !read files
   call init_pack_cea('uv')
   call read_conditions
   call calcTeq
   call out_plt
   print *,"done"
end program driver
