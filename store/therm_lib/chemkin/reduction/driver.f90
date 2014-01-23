program driver
   use const_chem
   use conditions
   implicit none

   !read files
   call init_therm
   call read_conditions
   call init_vrho
   call set_Tref

   print *,"original                    (ns,nr)=",ns,nr
   call reduction_species(.true.)
   print *,"after 1st species reduction (ns,nr)=",ns,nr
   call reduction_reactions
   print *,"after reaction reduction    (ns,nr)=",ns,nr
   call remove_remaining_reactions
   call reduction_species(.false.)
   print *,"after 2nd species reduction (ns,nr)=",ns,nr
   call out_cheminp
end program driver
