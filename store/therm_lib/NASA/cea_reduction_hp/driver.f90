program driver
   use conditions
   implicit none

   !read files
   call init_pack_cea('hp')
   call read_conditions
   call set_list_Teq_org
   call count_reactants
   call sort_species(.true.)
   call reduction
   call out_chem
   deallocate(list_p,list_T,list_Yf)
   deallocate(list_Teq_org)
   stop
end program driver
