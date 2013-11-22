program driver
   use const_chem
   implicit none
   double precision xi,E,T, MWave,kappa,mu,MWini
   double precision,dimension(1:ne)::n,b0

   call set_therm_data
   call set_trans_data

   xi=0d0
   T=300d0
   call calc_initial_prop(xi,T, n,E,MWini,b0)
   print *,"Hello",n,E,MWini,b0
   call calc_therm_flame_sheet(xi,E, T, MWave,kappa,mu)
   print *,"Hello",T,MWave,kappa,mu
end program driver
