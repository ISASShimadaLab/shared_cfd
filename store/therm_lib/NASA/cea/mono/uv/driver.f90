program driver
   use chem
   use chem_var
   implicit none
   character*100    buf
   double precision rho,E,T, MWave,kappa,mu
   double precision Y(2),DHi(2)
   double precision,dimension(max_ns)::n,vhi,Yv

   !read files
   call init_pack_cea('uv')

   !calc stoichiometry
   Y(1)=1d0/(1d0+of)
   Y(2)=1d0-Y(1)
   E  =  Ef*Y(1)+  Eo*Y(2)
   rho=1d0/(Y(1)/rhof+Y(2)/rhoo)
   n  = n_save(:,1,1,1)
   T=300d0
   call cea(rho,Y,E, T,n, MWave,kappa,mu,Yv,vhi,DHi)
   !print '(a,  f15.7)',"Y of f    ",Y(1)
   !print '(a,  f15.7)',"MWave     ",MWave
   !print '(a,  f15.7)',"kappa     ",kappa
   !print '(a,  f15.7)',"T         ",T
   !print '(a, es15.7)',"mu(mPoise)",mu*1d3*1d1
   !!print '(a,(10es9.1))',"Yv        ",Yv
   !!print '(a,(10es9.1))',"vhi       ",vhi
   !print '(a, es15.7)',"Enthalpy",  sum(Yv*vhi,1)
   !print '(a, es15.7)',"Enthalpy",  E+Ru*1d3/MWave*T

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
   print '(a, es15.7)',"        o/f ratio   =",of
   print '(a,  f15.7)',"        P(bar)      =",rho*Ru*1d3/MWave*T*1d-5
   print '(a,  f15.7)',"        rho(kg/m^3) =",rho
   print '(a,  f15.7)',"        T(K)        =",T
   print '(a, es15.7)',"        E(J/kg)     =",E
   print '(a,  f15.7)',"        average MW  =",MWave
   print '(a,  f15.7)',"        kappa       =",kappa
   print '(a, es15.7)',"        mu(mPoise)  =",mu*1d3*1d1
end program driver
