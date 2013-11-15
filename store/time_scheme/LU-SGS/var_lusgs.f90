module var_lusgs
   use grbl_prmtr
   double precision    Ap(1:dimq, 1:dimq, 1:ni, 1:nj)
   double precision    Am(1:dimq, 1:dimq, 1:ni, 1:nj)
   double precision    Bp(1:dimq, 1:dimq, 1:ni, 1:nj)
   double precision    Bm(1:dimq, 1:dimq, 1:ni, 1:nj)

   double precision alpha(1:ni, 1:nj)

   double precision dsci(1:ni,1:nj)
   double precision dscj(1:ni,1:nj)

   double precision vnci(1:2,1:ni,1:nj)
   double precision vncj(1:2,1:ni,1:nj)
end module var_lusgs

