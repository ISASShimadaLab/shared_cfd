program driver
   use chem
   use chem_var
   implicit none
   character*100    buf
   double precision rho,E,T, MWave,kappa,mu
   double precision Y(2),DHi(2)
   double precision,dimension(max_ns)::n,vhi,Yv
   integer i,Ntic

   !read files
   call init_pack_cea

   !calc stoichiometry
   !Y(1)=1d0/(1d0+of)
   Y(1)=1.0d0
   Y(2)=1d0-Y(1)
   E  =  Ef*Y(1)+  Eo*Y(2)
   !rho=rhof*Y(1)+rhoo*Y(2)
   rho=1d1
   n  =nf*Y(1)+no*Y(2)+initial_eps
   T=300d0
   call cea(rho,Y,E, T,n, MWave,kappa,mu)
   print '(a,  f15.7)',"Y of f    ",Y(1)
   print '(a,  f15.7)',"MWave     ",MWave
   print '(a,  f15.7)',"kappa     ",kappa
   print '(a,  f15.7)',"T         ",T
   print '(a, es15.7)',"mu(mPoise)",mu*1d3*1d1
   !print '(a,3es15.7)',"DHi       ",DHi
   !print '(a,3es15.7)',"Yv        ",Yv
   !print '(a, 3f15.7)',"vhi       ",vhi*1d-3

   !print '(a)',"*** pure oxidizer *****************************"
   !print '(a,  f15.7)',"        P(bar)      =",Po*1d-5
   !print '(a,  f15.7)',"        rho(kg/m^3) =",rhoo
   !print '(a,  f15.7)',"        T(K)        =",To
   !print '(a, es15.7)',"        E(J/kg)     =",Eo
   !print '(a,  f15.7)',"        average MW  =",MWo
   !print '(a,  f15.7)',"        kappa       =",kappao
   !print '(a, es15.7)',"        mu(mPoise)  =",muo*1d3*1d1
   !print '(a)',"*** pure fuel *********************************"
   !print '(a,  f15.7)',"        P(bar)      =",Pf*1d-5
   !print '(a,  f15.7)',"        rho(kg/m^3) =",rhof
   !print '(a,  f15.7)',"        T(K)        =",Tf
   !print '(a, es15.7)',"        E(J/kg)     =",Ef
   !print '(a,  f15.7)',"        average MW  =",MWf
   !print '(a,  f15.7)',"        kappa       =",kappaf
   !print '(a, es15.7)',"        mu(mPoise)  =",muf*1d3*1d1
   !print '(a)',"*** stoichiometry *****************************"
   !print '(a,es15.7)', "        o/f ratio   =",of
   !print '(a,  f15.7)',"        P(bar)      =",rho*Ru*1d3/MWave*T*1d-5
   !print '(a,  f15.7)',"        rho(kg/m^3) =",rho
   !print '(a,  f15.7)',"        T(K)        =",T
   !print '(a, es15.7)',"        E(J/kg)     =",E
   !print '(a,  f15.7)',"        average MW  =",MWave
   !print '(a,  f15.7)',"        kappa       =",kappa
   !print '(a, es15.7)',"        mu(mPoise)  =",mu*1d3*1d1

   !open(55,file="control.inp")
   !read(55,'(a)') buf
   !read(buf(26:),*) Ntic
   !close(55)

   !open(55,file="out.dat")

   !write(55,'(a)') "# MF = mass fraction, oxid = oxidizer, prod = product, &
   !                &br = before reaction, ar = after reaction, &
   !                &MW = Molecular Weight, kappa = specific heat ratio, &
   !                &vis = viscosity coefficient"
   !write(55,'(a1,a14,100a15)') "#","br fuel MF","br oxid MF",&
   !                     "ar fuel MF","ar oxid MF","ar prod MF",&
   !                     "Temperature(K)","Energy(J/kg)","MW (g/mol)",&
   !                     "kappa","vis(Pa*s)"
   !do i=0,Ntic
   !   Y(1)=dble(i)/dble(Ntic)
   !   Y(2)=1d0-Y(1)
   !   E  =  Ef*Y(1)+  Eo*Y(2)
   !   rho=rhof*Y(1)+rhoo*Y(2)
   !   T=To
   !   call flame_sheet(Y,E, T, MWave,kappa,mu,DHi,Yv,vhi)
   !   write(55,'(100es15.7)') Y,Yv,T,E,MWave,kappa,mu
   !end do
   !close(55)
end program driver
