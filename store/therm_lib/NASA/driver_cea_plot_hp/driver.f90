program driver
   use chem
   use chem_var
   implicit none
   character*100    buf
   double precision p,H,T, MWave,kappa,mu
   double precision Y(2)
   double precision,dimension(max_ns)::n,vhi,Yv
   integer i,Ntic
   logical flag_cea

   !read files
   call init_pack_cea('hp')

   !calc stoichiometry
   Y(1)=1d0/(1d0+of)
   Y(2)=1d0-Y(1)
   H  = Hf*Y(1)+Ho*Y(2)
   p  = pf*Y(1)+po*Y(2)
   n  = n_save(:,1,1,1)
   T=300d0
   call cea_hp(p,Y,H, T,n, MWave,kappa,mu,Yv,vhi,flag_cea)
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
   print '(a, es15.7)',"        H(kJ/kg)    =",Ho*1d-3
   print '(a,  f15.7)',"        average MW  =",MWo
   print '(a,  f15.7)',"        kappa       =",kappao
   print '(a, es15.7)',"        mu(mPoise)  =",muo*1d3*1d1
   print '(a)',"*** pure fuel *********************************"
   print '(a,  f15.7)',"        P(bar)      =",Pf*1d-5
   print '(a,  f15.7)',"        rho(kg/m^3) =",rhof
   print '(a,  f15.7)',"        T(K)        =",Tf
   print '(a, es15.7)',"        H(kJ/kg)    =",Hf*1d-3
   print '(a,  f15.7)',"        average MW  =",MWf
   print '(a,  f15.7)',"        kappa       =",kappaf
   print '(a, es15.7)',"        mu(mPoise)  =",muf*1d3*1d1
   print '(a)',"*** stoichiometry *****************************"
   print '(a, es15.7)',"        o/f ratio   =",of
   print '(a,  f15.7)',"        P(bar)      =",p*1d-5
   print '(a,  f15.7)',"        rho(kg/m^3) =",p/(Ru*1d3/MWave*T)
   print '(a,  f15.7)',"        T(K)        =",T
   print '(a, es15.7)',"        H(kJ/kg)    =",H*1d-3
   print '(a,  f15.7)',"        average MW  =",MWave
   print '(a,  f15.7)',"        kappa       =",kappa
   print '(a, es15.7)',"        mu(mPoise)  =",mu*1d3*1d1

   open(55,file="control.inp")
   read(55,'(a)') buf
   read(buf(26:),*) Ntic
   close(55)

   open(55,file="out.dat")

   write(55,'(a)') "# MF = mass fraction, oxid = oxidizer, prod = product, &
                   &br = before reaction, &
                   &MW = Molecular Weight, kappa = specific heat ratio, &
                   &vis = viscosity coefficient"
   write(55,'(a1,a14,100a15)') "#","br fuel MF","br oxid MF",&
                        "Temperature(K)","Enthalpy(J/kg)","MW (g/mol)",&
                        "kappa","vis(Pa*s)"

   n=no+initial_eps
   do i=0,Ntic
      Y(1)=dble(i)/dble(Ntic)
      Y(2)=1d0-Y(1)
      H = Hf*Y(1)+Ho*Y(2)
      p = pf*Y(1)+po*Y(2)
      call cea_hp(p,Y,H, T,n, MWave,kappa,mu,Yv,vhi,flag_cea)
      write(55,'(100es15.7)') Y,T,H,MWave,kappa,mu
   end do
   close(55)
end program driver
