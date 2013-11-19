module var_lusgs
   use grbl_prmtr
   double precision    Ap(1:dimq, 1:dimq, 0:nimax+1, 0:njmax+1,Nplane)
   double precision    Am(1:dimq, 1:dimq, 0:nimax+1, 0:njmax+1,Nplane)
   double precision    Bp(1:dimq, 1:dimq, 0:nimax+1, 0:njmax+1,Nplane)
   double precision    Bm(1:dimq, 1:dimq, 0:nimax+1, 0:njmax+1,Nplane)

   double precision alpha(1:nimax, 1:njmax,Nplane)
   double precision   phi(1:nimax, 1:njmax,Nplane)
   double precision pre_sca(4,1:nimax,1:njmax,Nplane)

   double precision dsci(0:nimax+1,0:njmax+1,Nplane)
   double precision dscj(0:nimax+1,0:njmax+1,Nplane)

   double precision vnci(1:2,1:nimax,1:njmax,Nplane)
   double precision vncj(1:2,1:nimax,1:njmax,Nplane)

   double precision pre1(dimq,1:nimax,1:njmax,Nplane)
   double precision dpdq(dimq,1:nimax,1:njmax,Nplane)
   double precision dqdt(dimq,1:nimax,1:njmax,Nplane)

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

