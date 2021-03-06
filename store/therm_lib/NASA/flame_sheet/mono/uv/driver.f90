program driver
   use chem
   use chem_var
   implicit none
   character*100    buf
   double precision rho,E,T, MWave,kappa,mu
   double precision Y(2),Yv(3),DHi(2),vhi(3)

   !read files
   call init_pack_flame_sheet

   !calc stoichiometry
   Y(1)=1d0/(1d0+of)
   Y(2)=1d0-Y(1)
   E  =  Ef*Y(1)+  Eo*Y(2)
   rho=rhof*Y(1)+rhoo*Y(2)
   T=To
   call flame_sheet(Y,E, T, MWave,kappa,mu,DHi,Yv,vhi)
   !print '(a,  f15.7)',"MWave     ",MWave
   !print '(a,  f15.7)',"kappa     ",kappa
   !print '(a,  f15.7)',"T         ",T
   !print '(a, es15.7)',"mu(mPoise)",mu*1d3*1d1
   !print '(a,3es15.7)',"DHi       ",DHi
   !print '(a,3es15.7)',"Yv        ",Yv
   !print '(a, 3f15.7)',"vhi       ",vhi*1d-3

   print '(a)',"*** pure oxidizer *****************************"
   print '(a,  f15.7)',"        P(bar)      =",Po*1d-5
   print '(a,  f15.7)',"        rho(kg/m^3) =",rhoo
   print '(a,  f15.7)',"        T(K)        =",To
   print '(a, es15.7)',"        E(J/kg)     =",Eo
   print '(a,  f15.7)',"        average MW  =",MWo
   print '(a,  f15.7)',"        kappa       =",kappao
   print '(a, es15.7)',"        mu(mPoise)  =",muo*1d3*1d1
   print '(a)',"*** pure fuel *********************************"
   print '(a,  f15.7)',"        P(bar)      =",Pf*1d-5
   print '(a,  f15.7)',"        rho(kg/m^3) =",rhof
   print '(a,  f15.7)',"        T(K)        =",Tf
   print '(a, es15.7)',"        E(J/kg)     =",Ef
   print '(a,  f15.7)',"        average MW  =",MWf
   print '(a,  f15.7)',"        kappa       =",kappaf
   print '(a, es15.7)',"        mu(mPoise)  =",muf*1d3*1d1
   print '(a)',"*** stoichiometry *****************************"
   print '(a,es15.7)', "        o/f ratio   =",of
   print '(a,  f15.7)',"        P(bar)      =",rho*Ru*1d3/MWave*T*1d-5
   print '(a,  f15.7)',"        rho(kg/m^3) =",rho
   print '(a,  f15.7)',"        T(K)        =",T
   print '(a, es15.7)',"        E(J/kg)     =",E
   print '(a,  f15.7)',"        average MW  =",MWave
   print '(a,  f15.7)',"        kappa       =",kappa
   print '(a, es15.7)',"        mu(mPoise)  =",mu*1d3*1d1
end program driver
