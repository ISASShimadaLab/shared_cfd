program driver
   use chem
   use chem_var
   implicit none
   double precision E,T, MWave,kappa,mu
   double precision Y(2),Yv(3),DHi(2),vhi(3)

   !read files
   call read_elmspc_to_use
   call set_therm_data
   call set_trans_data
   call read_fo_composition
   call read_p_composition

   !initialization
   call initialize_flame_sheet
   !Y(1)=0d0;E=Eo;T=To
   !Y(1)=1d0;E=Ef;T=Tf
   Y(1)=5d-1;E=0.5*(Eo+Ef);T=To
   Y(2)=1d0-Y(1)
   call flame_sheet(Y,E, T, MWave,kappa,mu,DHi,Yv,vhi)
   print '(a,  f15.7)',"MWave     ",MWave
   print '(a,  f15.7)',"kappa     ",kappa
   print '(a,  f15.7)',"T         ",T
   print '(a, es15.7)',"mu(mPoise)",mu*1d3*1d1
   print '(a,3es15.7)',"DHi       ",DHi
   print '(a,3es15.7)',"Yv        ",Yv
   print '(a, 3f15.7)',"vhi       ",vhi*1d-3
end program driver
