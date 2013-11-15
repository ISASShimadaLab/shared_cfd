module var_lusgs
   use grbl_prmtr
   double precision    Ap(1:dimq, 1:dimq, 0:ni+1, 0:nj+1)
   double precision    Am(1:dimq, 1:dimq, 0:ni+1, 0:nj+1)
   double precision    Bp(1:dimq, 1:dimq, 0:ni+1, 0:nj+1)
   double precision    Bm(1:dimq, 1:dimq, 0:ni+1, 0:nj+1)

   double precision   alpha(  1:ni,1:nj)
   double precision     phi(  1:ni,1:nj)
   double precision pre_sca(4,1:ni,1:nj)

   double precision dsci(1:ni,1:nj)
   double precision dscj(1:ni,1:nj)

   double precision vnci(1:2,1:ni,1:nj)
   double precision vncj(1:2,1:ni,1:nj)

   double precision pre1(dimq,1:ni,1:nj)
   double precision pre2(dimq,1:ni,1:nj)
   double precision dqdt(dimq,1:ni,1:nj)

   double precision cfl_tau
   double precision omega_max
   double precision omega_min
   double precision Dqmax
   double precision res_rate_warning
   double precision res_rate_error
   logical          ILwrite
   integer          NumInternalLoop

   logical          first_dual
   double precision res1
   double precision res
end module var_lusgs

