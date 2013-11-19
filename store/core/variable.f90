module variable
   use grbl_prmtr
   double precision    q(1:dimq,-1:nimax+2,-1:njmax+2,Nplane)
   double precision   qp(1:dimq,-1:nimax+2,-1:njmax+2,Nplane)
   double precision  qpp(1:dimq,-1:nimax+2,-1:njmax+2,Nplane)
   double precision    w(1:dimw,-1:nimax+2,-1:njmax+2,Nplane)

   double precision      DHi(nY, 1:nimax,   1:njmax  ,Nplane)
   double precision      vhi(nY,-1:nimax+2,-1:njmax+2,Nplane)

   double precision wHli(1:dimw,   0:nimax,   1:njmax,Nplane)
   double precision wHri(1:dimw,   0:nimax,   1:njmax,Nplane)
   double precision wHlj(1:dimw,   1:nimax,   0:njmax,Nplane)
   double precision wHrj(1:dimw,   1:nimax,   0:njmax,Nplane)

   double precision   TGi(1:dimq,   0:nimax,   1:njmax,Nplane)
   double precision   TGj(1:dimq,   1:nimax,   0:njmax,Nplane)
   double precision    Sq(1:dimq,   1:nimax,   1:njmax,Nplane)
   double precision  TGvi(1:dimq,   0:nimax,   1:njmax,Nplane)
   double precision  TGvj(1:dimq,   1:nimax,   0:njmax,Nplane)
   double precision   Svq(1:dimq,   1:nimax,   1:njmax,Nplane)

   double precision   Vol(1:nimax,1:njmax,Nplane)
   double precision  Area(1:nimax,1:njmax,Nplane)
   double precision   dsi(0:nimax,1:njmax,Nplane)
   double precision   dsj(1:nimax,0:njmax,Nplane)
   double precision geojacc(1:2,1:2,1:nimax,1:njmax,Nplane) !for viscous flow
   double precision geojaci(1:2,1:2,0:nimax,1:njmax,Nplane) !for viscous flow
   double precision geojacj(1:2,1:2,1:nimax,0:njmax,Nplane) !for viscous flow

   double precision  vni(1:2,0:nimax,1:njmax,Nplane)
   double precision  vnj(1:2,1:nimax,0:njmax,Nplane)

   double precision    xh(0:nimax  ,0:njmax  ,Nplane)
   double precision    rh(0:nimax  ,0:njmax  ,Nplane)
   double precision     x(0:nimax+1,0:njmax+1,Nplane)
   double precision     r(0:nimax+1,0:njmax+1,Nplane)

   double precision dt_mat(1:nimax,1:njmax,Nplane)
   double precision dt_grbl

   double precision cfl
   integer          OutPeriod
   double precision tt

   double precision max_step
   double precision min_RMS
end module variable

