module variable
   use grbl_prmtr
   double precision    q(1:dimq,-1:ni+2,-1:nj+2)
   double precision   qp(1:dimq,-1:ni+2,-1:nj+2)
   double precision  qpp(1:dimq,-1:ni+2,-1:nj+2)
   double precision    w(1:dimw,-1:ni+2,-1:nj+2)

   double precision      DHi(nY, 1:ni,   1:nj  )
   double precision      vhi(nY,-1:ni+2,-1:nj+2)

   double precision wHli(1:dimw,   0:ni,   1:nj)
   double precision wHri(1:dimw,   0:ni,   1:nj)
   double precision wHlj(1:dimw,   1:ni,   0:nj)
   double precision wHrj(1:dimw,   1:ni,   0:nj)

   double precision   TGi(1:dimq,   0:ni,   1:nj)
   double precision   TGj(1:dimq,   1:ni,   0:nj)
   double precision  TGvi(1:dimq,   0:ni,   1:nj)
   double precision  TGvj(1:dimq,   1:ni,   0:nj)

   double precision   Vol(1:ni,1:nj)
   double precision   dsi(0:ni,1:nj)
   double precision   dsj(1:ni,0:nj)
   double precision geojaci(1:2,1:2,0:ni,1:nj) !for viscous flow
   double precision geojacj(1:2,1:2,1:ni,0:nj) !for viscous flow

   double precision  vni(1:2,0:ni,1:nj)
   double precision  vnj(1:2,1:ni,0:nj)

   double precision    xh(0:ni  ,0:nj  )
   double precision    rh(0:ni  ,0:nj  )
   double precision     x(0:ni+1,0:nj+1)
   double precision     r(0:ni+1,0:nj+1)

   double precision dt_mat(1:ni,1:nj)
   double precision dt_grbl

   double precision cfl
   integer          OutPeriod
   double precision tt

   double precision max_step
   double precision min_RMS
end module variable

