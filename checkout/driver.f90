program driver
   use const_chem
   implicit none
   double precision xi,E,T, MWave,kappa,mu,MWini
   double precision,dimension(1:ne)::n,b0

   call read_fo_composition
   call set_therm_data
   call set_trans_data

   xi=0d0
   T=300d0
   call calc_initial_prop(xi,T, n,E,MWini,b0)
   print '(10es9.1)',n
   print '(10es9.1)',E
   print '(10es9.1)',MWini
   print '(10es9.1)',b0
   call calc_therm_flame_sheet(xi,E, T, MWave,kappa,mu)
   print '(10es9.1)',T
   print '(10es9.1)',MWave
   print '(10es9.1)',kappa
   print '(10es9.1)',mu
end program driver
