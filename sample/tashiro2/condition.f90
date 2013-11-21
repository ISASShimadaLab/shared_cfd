subroutine set_IC
   use grbl_prmtr
   use prmtr
   use variable
   use mod_mpi
   use gas
   implicit none
   integer i,j,plane

   do plane=nps,npe
      do i=nxs(plane),nxe(plane)
         do j=nys(plane),nye(plane)
            q(1,    i,j,plane) = 1.013d5/(R_gas*300d0)
            q(2,    i,j,plane) = 0d0
            q(3,    i,j,plane) = 0d0
            q(4,    i,j,plane) = 1.013d5/(kappa_gas-1d0)
            w(4,    i,j,plane) = 1.013d5
            w(indxg,i,j,plane) = kappa_gas
            w(indxR,i,j,plane) = R_gas
         end do
      end do
   end do
end subroutine set_IC

subroutine set_BC(step)
   use grbl_prmtr
   use prmtr
   use variable
   use mod_mpi
   use gas
   implicit none
   real(8) rho,u,v,p,T,a
   integer,intent(in)::step
   integer i,j,plane

   integer,parameter::DLength=dimw+nY !for MPI Communication
   integer,parameter::MaxMPIcomm = 3

   call section_exchange
   call MPI_COMMUNICATIONS_CUT_LINE


   !boundary right and left
   do plane=nps,npe
      select case(plane)
      case(  1)
         if(gx(plane) .eq. 1) then
            !i=1/2
            do j=max(    1,nys(plane)),min(nye(plane),   40)
               if (w(2,1,j,plane)<=0.d0) then
                  rho=w(1,1,j,plane)
                  u  =w(2,1,j,plane)
                  v  =w(3,1,j,plane)
                  P  =101300.d0
                  T  =P/(rho*R_gas)
                  w(1,     -1:0,j,plane) = rho
                  w(2,     -1:0,j,plane) = u
                  w(3,     -1:0,j,plane) = v
                  w(4,     -1:0,j,plane) = P
                  w(5,     -1:0,j,plane) = 1.d0
                  w(indxg ,-1:0,j,plane) = kappa_gas
                  w(indxht,-1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                  w(indxR ,-1:0,j,plane) = R_gas
                  w(indxMu,-1:0,j,plane) = nu_gas*rho
                  vhi(:,   -1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               else
                  rho=1.17d0
                  u=0.d0
                  v=0.d0
                  P=w(4,1,j,plane)
                  T=P/(rho*R_gas)
                  w(1,     -1:0,j,plane) = rho
                  w(2,     -1:0,j,plane) = u
                  w(3,     -1:0,j,plane) = v
                  w(4,     -1:0,j,plane) = P
                  w(5,     -1:0,j,plane) = 1.d0
                  w(indxg ,-1:0,j,plane) = kappa_gas
                  w(indxht,-1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                  w(indxR ,-1:0,j,plane) = R_gas
                  w(indxMu,-1:0,j,plane) = nu_gas*rho
                  vhi(:,   -1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               end if
            end do
         end if

      case(  2)
         if(gx(plane) .eq. 1) then
            !i=1/2
            do j=max(    1,nys(plane)),min(nye(plane),   48)
               if (j<=29) then
                  select case(j)
                  case(1)
                     rho=3.167383062d0     ;u=995.246420859d0 ;v=2.104905015d0     ;P=279711.27439433d0
                  case(2)
                     rho=3.168145503d0     ;u=995.312408515d0   ;v=6.22470697d0      ;P=279721.074085528d0
                  case(3)
                     rho=3.168211591d0     ;u=995.3054576455d0  ;v=8.2655706405d0    ;P=279727.874991708d0
                  case(4)
                     rho=3.168277679d0     ;u=995.298506776d0   ;v=10.306434311d0    ;P=279734.675897888d0
                  case(5)
                     rho=3.1689936995d0    ;u=995.292219636d0   ;v=16.2576495005d0   ;P=279760.624990298d0
                  case(6)
                     rho=3.169497614d0     ;u=995.23333168d0    ;v=22.030979307d0    ;P=279781.631660681d0
                  case(7)
                     rho=3.169897687d0     ;u=995.151532104d0   ;v=25.715833684d0    ;P=279811.37740765d0
                  case(8)
                     rho=3.170420326d0     ;u=995.0226565615d0  ;v=31.04058071d0     ;P=279847.641012173d0
                  case(9)
                     rho=3.171405707d0     ;u=994.799489601d0   ;v=37.6801195735d0   ;P=279954.075724197d0
                  case(10)
                     rho=3.172087204d0     ;u=994.596339691d0   ;v=42.308989928d0    ;P=280036.728946251d0
                  case(11)
                     rho=3.1723267955d0    ;u=994.383528276d0   ;v=46.6190467475d0   ;P=280082.644065858d0
                  case(12)
                     rho=3.1724257415d0    ;u=994.066543096d0   ;v=51.9049782355d0   ;P=280131.593504522d0
                  case(13)
                     rho=3.172157651d0     ;u=993.7456617005d0  ;v=56.6492932665d0   ;P=280153.425768816d0
                  case(14)
                     rho=3.1715049205d0    ;u=993.430614839d0   ;v=60.8698496025d0   ;P=280141.946565951d0
                  case(15)
                     rho=3.170509638d0     ;u=993.1343573335d0  ;v=64.5926269745d0   ;P=280102.546507717d0
                  case(16)
                     rho=3.1692652575d0    ;u=992.8607101285d0  ;v=67.8465367605d0   ;P=280048.376660701d0
                  case(17)
                     rho=3.1679055d0       ;u=992.612845238d0   ;v=70.6697366185d0   ;P=279995.670290441d0
                  case(18)
                     rho=3.16634615966667d0;u=992.360926137333d0;v=73.7051589543333d0;P=279960.956428361d0
                  case(19)
                     rho=3.1626102055d0    ;u=992.095260013d0   ;v=76.329610336d0    ;P=279668.908797877d0
                  case(20)
                     rho=3.1584860415d0    ;u=991.9684951685d0  ;v=78.271882902d0    ;P=279319.031066d0
                  case(21)
                     rho=3.1520137275d0    ;u=991.9258894275d0  ;v=80.1874876325d0   ;P=278707.49302713d0
                  case(22)
                     rho=3.1417842165d0    ;u=991.932726448d0   ;v=82.365244924d0    ;P=277700.408719343d0
                  case(23)
                     rho=3.124435417d0     ;u=992.2623310925d0  ;v=85.347157343d0    ;P=276188.17532897d0
                  case(24)
                     rho=3.0789424835d0    ;u=992.710413077d0   ;v=89.5851639675d0   ;P=272243.031591964d0
                  case(25)
                     rho=2.876491363d0     ;u=983.0899833915d0  ;v=96.047096339d0    ;P=263850.099636762d0
                  case(26)
                     rho=2.504247941d0     ;u=955.522111511d0   ;v=104.335599984d0   ;P=251938.389405018d0
                  case(27)
                     rho=1.933276308d0     ;u=862.9919822615d0  ;v=111.907831096d0   ;P=234951.773684619d0
                  case(28)
                     rho=1.2880461265d0    ;u=663.542857364d0   ;v=121.1072499615d0  ;P=196628.940012597d0
                  case(29)
                     rho=1.036782427d0     ;u=467.764126376d0   ;v=103.470942947d0   ;P=181811.987230449d0
                  end select
                  T=P/(rho*R_gas)
                  w(1,     -1:0,j,plane) = rho
                  w(2,     -1:0,j,plane) = u
                  w(3,     -1:0,j,plane) = v
                  w(4,     -1:0,j,plane) = P
                  w(5,     -1:0,j,plane) = 1.d0
                  w(indxg ,-1:0,j,plane) = kappa_gas
                  w(indxht,-1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                  w(indxR ,-1:0,j,plane) = R_gas
                  w(indxMu,-1:0,j,plane) = nu_gas*rho
                  vhi(:,   -1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               else
                  w(:,   0,j,plane) = w(:,   1,j,plane)
                  w(2:3, 0,j,plane) =-w(2:3, 1,j,plane)
                  vhi(:, 0,j,plane) = vhi(:, 1,j,plane)

                  w(:,  -1,j,plane) = w(:,   0,j,plane)
                  vhi(:,-1,j,plane) = vhi(:, 0,j,plane)
               end if
            end do
         end if

         if(gx(plane) .eq. ngx(plane)) then
            !i=ni+1/2
            do j=max(    1,nys(plane)),min(nye(plane),   88)
               T=w(4,ni(plane),j,plane)/(w(1,ni(plane),j,plane)*R_gas) 
               a=sqrt(kappa_gas*R_gas*T)
               if (w(2,ni(plane),j,plane)>=0.d0) then
                  if (w(2,ni(plane),j,plane)>=a) then
                     w(:,  ni(plane)+1,j,plane) = w(:,  ni(plane),j,plane)
                     vhi(:,ni(plane)+1,j,plane) = vhi(:,ni(plane),j,plane)

                     w(:,  ni(plane)+2,j,plane) = w(:,  ni(plane)+1,j,plane)
                     vhi(:,ni(plane)+2,j,plane) = vhi(:,ni(plane)+1,j,plane)
                  else
                     rho= w(1,ni(plane),j,plane)
                     u  = w(2,ni(plane),j,plane)
                     v  = w(3,ni(plane),j,plane)
                     P  = 101300.d0
                     T  = P/(rho*R_gas)
                     w(1     ,ni(plane)+1:ni(plane)+2,j,plane) = rho
                     w(2     ,ni(plane)+1:ni(plane)+2,j,plane) = u
                     w(3     ,ni(plane)+1:ni(plane)+2,j,plane) = v
                     w(4     ,ni(plane)+1:ni(plane)+2,j,plane) = P
                     w(5     ,ni(plane)+1:ni(plane)+2,j,plane) = 1.d0
                     w(indxg ,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas
                     w(indxht,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                     w(indxR ,ni(plane)+1:ni(plane)+2,j,plane) = R_gas
                     w(indxMu,ni(plane)+1:ni(plane)+2,j,plane) = nu_gas*rho
                     vhi(:   ,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
                  end if
               else
                  rho= 1.17d0
                  u  = 0.d0
                  v  = 0.d0
                  P  = w(4,ni(plane),j,plane)
                  T  = P/(rho*R_gas)
                  w(1     ,ni(plane)+1:ni(plane)+2,j,plane) = rho
                  w(2     ,ni(plane)+1:ni(plane)+2,j,plane) = u
                  w(3     ,ni(plane)+1:ni(plane)+2,j,plane) = v
                  w(4     ,ni(plane)+1:ni(plane)+2,j,plane) = P
                  w(5     ,ni(plane)+1:ni(plane)+2,j,plane) = 1.d0
                  w(indxg ,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas
                  w(indxht,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                  w(indxR ,ni(plane)+1:ni(plane)+2,j,plane) = R_gas
                  w(indxMu,ni(plane)+1:ni(plane)+2,j,plane) = nu_gas*rho
                  vhi(:   ,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               end if
            end do
         end if
      case(  3)
         if(gx(plane) .eq. 1) then
            !i=1/2
            do j=max(    1,nys(plane)),min(nye(plane),   39)
               if (w(2,1,j,plane)<=0.d0) then
                  rho=w(1,1,j,plane)
                  u  =w(2,1,j,plane)
                  v  =w(3,1,j,plane)
                  P  =101300.d0
                  T  =P/(rho*R_gas)
                  w(1     ,-1:0,j,plane) = rho
                  w(2     ,-1:0,j,plane) = u
                  w(3     ,-1:0,j,plane) = v
                  w(4     ,-1:0,j,plane) = P
                  w(5     ,-1:0,j,plane) = 1.d0
                  w(indxg ,-1:0,j,plane) = kappa_gas
                  w(indxht,-1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                  w(indxR ,-1:0,j,plane) = R_gas
                  w(indxMu,-1:0,j,plane) = nu_gas*rho
                  vhi(:   ,-1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               else
                  rho=1.17d0
                  u=0.d0
                  v=0.d0
                  P=w(4,1,j,plane)
                  T=P/(rho*R_gas)
                  w(1     ,-1:0,j,plane) = rho
                  w(2     ,-1:0,j,plane) = u
                  w(3     ,-1:0,j,plane) = v
                  w(4     ,-1:0,j,plane) = P
                  w(5     ,-1:0,j,plane) = 1.d0
                  w(indxg ,-1:0,j,plane) = kappa_gas
                  w(indxht,-1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                  w(indxR ,-1:0,j,plane) = R_gas
                  w(indxMu,-1:0,j,plane) = nu_gas*rho
                  vhi(:   ,-1:0,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               end if
            end do
         end if

         if(gx(plane) .eq. ngx(plane)) then
            !i=ni+1/2
            do j=max(    1,nys(plane)),min(nye(plane),   39)
               if (w(2,ni(plane),j,plane)>=0.d0) then
                  rho= w(1,ni(plane),j,plane)
                  u  = w(2,ni(plane),j,plane)
                  v  = w(3,ni(plane),j,plane)
                  P  = 101300.d0
                  T  = P/(rho*R_gas)
                  w(1     ,ni(plane)+1:ni(plane)+2,j,plane) = rho
                  w(2     ,ni(plane)+1:ni(plane)+2,j,plane) = u
                  w(3     ,ni(plane)+1:ni(plane)+2,j,plane) = v
                  w(4     ,ni(plane)+1:ni(plane)+2,j,plane) = P
                  w(5     ,ni(plane)+1:ni(plane)+2,j,plane) = 1.d0
                  w(indxg ,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas
                  w(indxht,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                  w(indxR ,ni(plane)+1:ni(plane)+2,j,plane) = R_gas
                  w(indxMu,ni(plane)+1:ni(plane)+2,j,plane) = nu_gas*rho
                  vhi(:   ,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               else
                  rho= 1.17d0
                  u  = 0.d0
                  v  = 0.d0
                  P  = w(4,ni(plane),j,plane)
                  T  = P/(rho*R_gas)
                  w(1     ,ni(plane)+1:ni(plane)+2,j,plane) = rho
                  w(2     ,ni(plane)+1:ni(plane)+2,j,plane) = u
                  w(3     ,ni(plane)+1:ni(plane)+2,j,plane) = v
                  w(4     ,ni(plane)+1:ni(plane)+2,j,plane) = P
                  w(5     ,ni(plane)+1:ni(plane)+2,j,plane) = 1.d0
                  w(indxg ,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas
                  w(indxht,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                  w(indxR ,ni(plane)+1:ni(plane)+2,j,plane) = R_gas
                  w(indxMu,ni(plane)+1:ni(plane)+2,j,plane) = nu_gas*rho
                  vhi(:   ,ni(plane)+1:ni(plane)+2,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               end if
            end do
         end if
      end select
   end do

   call MPI_COMMUNICATIONS_I_DIRECTION

   !boundary upper and lower
   do plane=nps,npe
      select case(plane)
      case(  1)
         if(gy(plane) .eq. 1) then
            !j=1/2
            do i=max(    1,nxs(plane)),min(nxe(plane),   59)
               w(:,  i, 0,plane) = w(:,  i,1,plane)
               w(2:3,i, 0,plane) =-w(2:3,i,1,plane)
               vhi(:,i, 0,plane) = vhi(:,i,1,plane)

               w(:,  i,-1,plane) = w(:,  i,0,plane)
               vhi(:,i,-1,plane) = vhi(:,i,0,plane)
            end do
         end if

      case(  2)
         if(gy(plane) .eq. 1) then
            !j=1/2
            do i=max(    1,nxs(plane)),min(nxe(plane),  309)
               w(:,  i, 0,plane) = w(:,  i,1,plane)
               w(3,  i, 0,plane) =-w(3,  i,1,plane)
               vhi(:,i, 0,plane) = vhi(:,i,1,plane)

               w(:,  i,-1,plane) = w(:,  i,0,plane)
               vhi(:,i,-1,plane) = vhi(:,i,0,plane)
            end do
         end if

         if(gy(plane) .eq.  ngy(plane)) then
            !j=nj+1/2
            do i=max(   64,nxs(plane)),min(nxe(plane),  260)
               w(:,  i,nj(plane)+1,plane) = w(:,  i,nj(plane)  ,plane)
               w(2:3,i,nj(plane)+1,plane) =-w(2:3,i,nj(plane)  ,plane)
               vhi(:,i,nj(plane)+1,plane) = vhi(:,i,nj(plane)  ,plane)

               w(:,  i,nj(plane)+2,plane) = w(:,  i,nj(plane)+1,plane)
               vhi(:,i,nj(plane)+2,plane) = vhi(:,i,nj(plane)+1,plane)
            end do
         end if
      case(  3)
         if(gy(plane) .eq. 1) then
            !j=1/2
            do i=max(  123,nxs(plane)),min(nxe(plane),  319)
               w(:,  i, 0,plane) = w(:,  i,1,plane) 
               w(2:3,i, 0,plane) =-w(2:3,i,1,plane) 
               vhi(:,i, 0,plane) = vhi(:,i,1,plane) 

               w(:,  i,-1,plane) = w(:,  i,0,plane) 
               vhi(:,i,-1,plane) = vhi(:,i,0,plane) 
            end do
         end if

         if(gy(plane) .eq.  ngy(plane)) then
            !j=nj+1/2
            do i=max(    1,nxs(plane)),min(nxe(plane),  368)
               if (w(3,  i,nj(plane),plane)>=0.d0) then
                     rho= w(1,i,nj(plane),plane)
                     u  = w(2,i,nj(plane),plane)
                     v  = w(3,i,nj(plane),plane)
                     P  = 101300.d0
                     T  = P/(rho*R_gas)
                     w(1     ,i,nj(plane)+1:nj(plane)+2,plane) = rho
                     w(2     ,i,nj(plane)+1:nj(plane)+2,plane) = u
                     w(3     ,i,nj(plane)+1:nj(plane)+2,plane) = v
                     w(4     ,i,nj(plane)+1:nj(plane)+2,plane) = P
                     w(5     ,i,nj(plane)+1:nj(plane)+2,plane) = 1.d0
                     w(indxg ,i,nj(plane)+1:nj(plane)+2,plane) = kappa_gas
                     w(indxht,i,nj(plane)+1:nj(plane)+2,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                     w(indxR ,i,nj(plane)+1:nj(plane)+2,plane) = R_gas
                     w(indxMu,i,nj(plane)+1:nj(plane)+2,plane) = nu_gas*rho
                     vhi(:   ,i,nj(plane)+1:nj(plane)+2,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               else
                  rho= 1.17d0
                  u  = 0.d0
                  v  = 0.d0
                  P  = w(4,i,nj(plane),plane)
                  T  = P/(rho*R_gas)
                  w(1     ,i,nj(plane)+1:nj(plane)+2,plane) = rho
                  w(2     ,i,nj(plane)+1:nj(plane)+2,plane) = u
                  w(3     ,i,nj(plane)+1:nj(plane)+2,plane) = v
                  w(4     ,i,nj(plane)+1:nj(plane)+2,plane) = P
                  w(5     ,i,nj(plane)+1:nj(plane)+2,plane) = 1.d0
                  w(indxg ,i,nj(plane)+1:nj(plane)+2,plane) = kappa_gas
                  w(indxht,i,nj(plane)+1:nj(plane)+2,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T+0.5d0*(u**2+v**2)
                  w(indxR ,i,nj(plane)+1:nj(plane)+2,plane) = R_gas
                  w(indxMu,i,nj(plane)+1:nj(plane)+2,plane) = nu_gas*rho
                  vhi(:   ,i,nj(plane)+1:nj(plane)+2,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
               end if
            end do
         end if
      end select
   end do

   call MPI_COMMUNICATIONS_J_DIRECTION

   call set_corners
contains
subroutine section_exchange!{{{
   implicit none
   integer inc,i
   integer geo,width,order,p1,c1,s1,pm1,p2,c2,s2,pm2

   do i=1,num_cut_copro
      geo  = cut_copro( 1,i)
      width= cut_copro( 2,i)
      order= cut_copro( 3,i)
      p1   = cut_copro( 4,i)
      c1   = cut_copro( 5,i)
      s1   = cut_copro( 6,i)
      pm1  = cut_copro( 7,i)
      p2   = cut_copro( 8,i)
      c2   = cut_copro( 9,i)
      s2   = cut_copro(10,i)
      pm2  = cut_copro(11,i)
      select case(geo)
      case(1) !i-i
         do inc=0,width-1
            w(:,  c1+pm1,  s1+inc,      p1) = w(:,  c2,    s2+inc*order,p2)
            vhi(:,c1+pm1,  s1+inc,      p1) = vhi(:,c2,    s2+inc*order,p2)
            w(:,  c1+pm1*2,s1+inc,      p1) = w(:,  c2-pm2,s2+inc*order,p2)
            vhi(:,c1+pm1*2,s1+inc,      p1) = vhi(:,c2-pm2,s2+inc*order,p2)

            w(:,  c2+pm2,  s2+inc*order,p2) = w(:,  c1,    s1+inc,      p1)
            vhi(:,c2+pm2,  s2+inc*order,p2) = vhi(:,c1,    s1+inc,      p1)
            w(:,  c2+pm2*2,s2+inc*order,p2) = w(:,  c1-pm1,s1+inc,      p1)
            vhi(:,c2+pm2*2,s2+inc*order,p2) = vhi(:,c1-pm1,s1+inc,      p1)
         end do
      case(2) !i-j
         do inc=0,width-1
            w(:,  c1+pm1,  s1+inc,      p1) = w(:,  s2+inc*order,c2,    p2)
            vhi(:,c1+pm1,  s1+inc,      p1) = vhi(:,s2+inc*order,c2,    p2)
            w(:,  c1+pm1*2,s1+inc,      p1) = w(:,  s2+inc*order,c2-pm2,p2)
            vhi(:,c1+pm1*2,s1+inc,      p1) = vhi(:,s2+inc*order,c2-pm2,p2)

            w(:,  s2+inc*order,c2+pm2,  p2) = w(:,  c1,    s1+inc,      p1)
            vhi(:,s2+inc*order,c2+pm2,  p2) = vhi(:,c1,    s1+inc,      p1)
            w(:,  s2+inc*order,c2+pm2*2,p2) = w(:,  c1-pm1,s1+inc,      p1)
            vhi(:,s2+inc*order,c2+pm2*2,p2) = vhi(:,c1-pm1,s1+inc,      p1)
         end do
      case(3) !j-i
         do inc=0,width-1
            w(:,  s1+inc,      c1+pm1,  p1) = w(:,  c2,    s2+inc*order,p2)
            vhi(:,s1+inc,      c1+pm1,  p1) = vhi(:,c2,    s2+inc*order,p2)
            w(:,  s1+inc,      c1+pm1*2,p1) = w(:,  c2-pm2,s2+inc*order,p2)
            vhi(:,s1+inc,      c1+pm1*2,p1) = vhi(:,c2-pm2,s2+inc*order,p2)

            w(:,  c2+pm2,  s2+inc*order,p2) = w(:,  s1+inc,      c1,    p1)
            vhi(:,c2+pm2,  s2+inc*order,p2) = vhi(:,s1+inc,      c1,    p1)
            w(:,  c2+pm2*2,s2+inc*order,p2) = w(:,  s1+inc,      c1-pm1,p1)
            vhi(:,c2+pm2*2,s2+inc*order,p2) = vhi(:,s1+inc,      c1-pm1,p1)
         end do
      case(4) !j-j
         do inc=0,width-1
            w(:,  s1+inc,      c1+pm1,  p1) = w(:,  s2+inc*order,c2,    p2)
            vhi(:,s1+inc,      c1+pm1,  p1) = vhi(:,s2+inc*order,c2,    p2)
            w(:,  s1+inc,      c1+pm1*2,p1) = w(:,  s2+inc*order,c2-pm2,p2)
            vhi(:,s1+inc,      c1+pm1*2,p1) = vhi(:,s2+inc*order,c2-pm2,p2)
                                                                        
            w(:,  s2+inc*order,c2+pm2,  p2) = w(:,  s1+inc,      c1,    p1)
            vhi(:,s2+inc*order,c2+pm2,  p2) = vhi(:,s1+inc,      c1,    p1)
            w(:,  s2+inc*order,c2+pm2*2,p2) = w(:,  s1+inc,      c1-pm1,p1)
            vhi(:,s2+inc*order,c2+pm2*2,p2) = vhi(:,s1+inc,      c1-pm1,p1)
         end do
      end select
   end do
end subroutine section_exchange!}}}
subroutine MPI_COMMUNICATIONS_CUT_LINE!{{{
   implicit none
   double precision  tmp(DLength*2*bwmax)
   double precision tmpc(DLength*2*bwmax,MaxMPIcomm)
   integer cnt
   integer,dimension(MaxMPIcomm)::cntc,toProc,ierrc,ireqc
   integer inc,i
   integer cind,width,order,p,c,s,pm

   !send
   do i = 1,NumMPIComm
      cind      = MPIComm(1,i)
      width     = MPIComm(2,i)
      order     = MPIComm(3,i)
      p         = MPIComm(4,i)
      c         = MPIComm(5,i)
      s         = MPIComm(6,i)
      pm        = MPIComm(7,i)
      toProc(i) = MPIComm(8,i)

      cntc(i)=0
      if(cind .eq. 0) then
         do inc=0,width-1
            !print '(a7,4i3,9es9.1)',"send",myid,p,c,s+inc*order,w(:,  c,   s+inc*order,p)
            tmpc(cntc(i)+1:cntc(i)+nY,  i) = vhi(:,c,   s+inc*order,p); cntc(i)=cntc(i)+nY
            tmpc(cntc(i)+1:cntc(i)+dimw,i) = w(:,  c,   s+inc*order,p); cntc(i)=cntc(i)+dimw
            !print '(a7,4i3,9es9.1)',"send",myid,p,c-pm,s+inc*order,w(:,  c-pm,   s+inc*order,p)
            tmpc(cntc(i)+1:cntc(i)+nY,  i) = vhi(:,c-pm,s+inc*order,p); cntc(i)=cntc(i)+nY
            tmpc(cntc(i)+1:cntc(i)+dimw,i) = w(:,  c-pm,s+inc*order,p); cntc(i)=cntc(i)+dimw
         end do
      else
         do inc=0,width-1
            tmpc(cntc(i)+1:cntc(i)+nY,  i) = vhi(:,s+inc*order,c,   p); cntc(i)=cntc(i)+nY
            tmpc(cntc(i)+1:cntc(i)+dimw,i) = w(:,  s+inc*order,c,   p); cntc(i)=cntc(i)+dimw
            tmpc(cntc(i)+1:cntc(i)+nY,  i) = vhi(:,s+inc*order,c-pm,p); cntc(i)=cntc(i)+nY
            tmpc(cntc(i)+1:cntc(i)+dimw,i) = w(:,  s+inc*order,c-pm,p); cntc(i)=cntc(i)+dimw
         end do
      end if
      call MPI_Isend(tmpc(:,i),cntc(i),MPI_DOUBLE_PRECISION,toProc(i),0,MPI_COMM_WORLD,ireqc(i),ierrc(i))
   end do

   !receive
   do i = 1,NumMPIComm
      cind      = MPIComm(1,i)
      width     = MPIComm(2,i)
      order     = MPIComm(3,i)
      p         = MPIComm(4,i)
      c         = MPIComm(5,i)
      s         = MPIComm(6,i)
      pm        = MPIComm(7,i)

      cnt=DLength*2*width
      call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,toProc(i),0,MPI_COMM_WORLD, istatus, ierr)
      cnt=0
      if(cind .eq. 0) then
         do inc=0,width-1
            vhi(:,c+pm,  s+inc*order,p) = tmp(cnt+1:cnt+nY);   cnt=cnt+nY
            w(:,  c+pm,  s+inc*order,p) = tmp(cnt+1:cnt+dimw); cnt=cnt+dimw
            !print '(a7,4i3,9es9.1)',"receive",myid,p,c+pm,s+inc*order,w(:,  c+pm,  s+inc*order,p)
            vhi(:,c+pm*2,s+inc*order,p) = tmp(cnt+1:cnt+nY);   cnt=cnt+nY
            w(:,  c+pm*2,s+inc*order,p) = tmp(cnt+1:cnt+dimw); cnt=cnt+dimw
            !print '(a7,4i3,9es9.1)',"receive",myid,p,c+pm*2,s+inc*order,w(:,c+pm*2,s+inc*order,p)
         end do
      else
         do inc=0,width-1
            vhi(:,s+inc*order,c+pm,  p) = tmp(cnt+1:cnt+nY);   cnt=cnt+nY
            w(:,  s+inc*order,c+pm,  p) = tmp(cnt+1:cnt+dimw); cnt=cnt+dimw
            vhi(:,s+inc*order,c+pm*2,p) = tmp(cnt+1:cnt+nY);   cnt=cnt+nY
            w(:,  s+inc*order,c+pm*2,p) = tmp(cnt+1:cnt+dimw); cnt=cnt+dimw
         end do
      end if
   end do

   !wait
   do i = 1,NumMPIComm
      call MPI_Wait(ireqc(i),istatus,ierrc(i))
   end do
end subroutine MPI_COMMUNICATIONS_CUT_LINE!}}}
subroutine MPI_COMMUNICATIONS_I_DIRECTION!{{{
   implicit none
   double precision  tmp(DLength*2*(bwmax+2))
   double precision tmps(DLength*2*(bwmax+2),Nplane)
   integer ii,jj
   integer cnt,cnts(Nplane),flag(Nplane)

   !grid west to east
   !send
   do plane = nps,npe
      if(gx(plane) .ne. ngx(plane)) then
         cnts(plane)=0
         do jj=nys(plane),nye(plane)
            do ii=nxe(plane)-1,nxe(plane)
               tmps(cnts(plane)+1:cnts(plane)+nY,  plane) =vhi(:,ii,jj,plane)
               cnts(plane)=cnts(plane)+nY
               tmps(cnts(plane)+1:cnts(plane)+dimw,plane) =w(  :,ii,jj,plane)
               cnts(plane)=cnts(plane)+dimw
            end do
         end do
         flag(plane)=plane
         call MPI_Isend(tmps(:,plane),cnts(plane),MPI_DOUBLE_PRECISION,ge(plane),flag(plane),&
                         MPI_COMM_WORLD,ireq(plane),ierrs(plane))
      end if
   end do

   !receive
   do plane = nps,npe
      if(gx(plane) .ne. 1) then
         cnt=DLength*2*bwy(plane)
         call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gw(plane),plane,&
                         MPI_COMM_WORLD, istatus, ierr)

         cnt=0
         do jj=nys(plane),nye(plane)
            do ii=nxs(plane)-2,nxs(plane)-1
               vhi(:,ii,jj,plane)=tmp(cnt+1:cnt+nY)
               cnt=cnt+nY
               w(  :,ii,jj,plane)=tmp(cnt+1:cnt+dimw)
               cnt=cnt+dimw
            end do
         end do
      end if
   end do

   do plane = nps,npe
      if(gx(plane) .ne. ngx(plane)) call MPI_Wait(ireq(plane),istatus,ierrs(plane))
   end do

   !grid east to west
   !send
   do plane = nps,npe
      if(gx(plane) .ne. 1) then
         cnts(plane)=0
         do jj=nys(plane),nye(plane)
            do ii=nxs(plane),nxs(plane)+1
               tmps(cnts(plane)+1:cnts(plane)+nY  ,plane) =vhi(:,ii,jj,plane)
               cnts(plane)=cnts(plane)+nY                                   
               tmps(cnts(plane)+1:cnts(plane)+dimw,plane) =w(  :,ii,jj,plane)
               cnts(plane)=cnts(plane)+dimw
            end do
         end do
         flag(plane)=plane
         call MPI_Isend(tmps(:,plane),cnts(plane),MPI_DOUBLE_PRECISION,gw(plane),flag(plane),&
                         MPI_COMM_WORLD,ireq(plane),ierrs(plane))
      end if
   end do

   !receive
   do plane = nps,npe
      if(gx(plane) .ne. ngx(plane)) then
         cnt=DLength*2*bwy(plane)
         call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,ge(plane),plane,&
                         MPI_COMM_WORLD, istatus, ierr)

         cnt=0
         do jj=nys(plane),nye(plane)
            do ii=nxe(plane)+1,nxe(plane)+2
               vhi(:,ii,jj,plane)=tmp(cnt+1:cnt+nY)
               cnt=cnt+nY                           
               w(  :,ii,jj,plane)=tmp(cnt+1:cnt+dimw)
               cnt=cnt+dimw
            end do
         end do
      end if
   end do

   do plane = nps,npe
      if(gx(plane) .ne. 1) call MPI_Wait(ireq(plane),istatus,ierrs(plane))
   end do
end subroutine MPI_COMMUNICATIONS_I_DIRECTION!}}}
subroutine MPI_COMMUNICATIONS_J_DIRECTION!{{{
   implicit none
   double precision  tmp(DLength*2*(bwmax+2))
   double precision tmps(DLength*2*(bwmax+2),Nplane)
   integer ii,jj
   integer cnt,cnts(Nplane),flag(Nplane)

   !grid south to north
   !send
   do plane = nps,npe
      if(gy(plane) .ne. ngy(plane)) then
         cnts(plane)=0
         do jj=nye(plane)-1,nye(plane)
            do ii=nxs(plane)-1,nxe(plane)+1
               tmps(cnts(plane)+1:cnts(plane)+nY,  plane) =vhi(:,ii,jj,plane)
               cnts(plane)=cnts(plane)+nY
               tmps(cnts(plane)+1:cnts(plane)+dimw,plane) =w(  :,ii,jj,plane)
               cnts(plane)=cnts(plane)+dimw
            end do
         end do
         flag(plane)=plane
         call MPI_Isend(tmps(:,plane),cnts(plane),MPI_DOUBLE_PRECISION,gn(plane),flag(plane),&
                         MPI_COMM_WORLD,ireq(plane),ierrs(plane))
      end if
   end do

   !receive
   do plane = nps,npe
      if(gy(plane) .ne. 1) then
         cnt=DLength*2*(bwx(plane)+2)
         call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gs(plane),plane,&
                         MPI_COMM_WORLD, istatus, ierr)

         cnt=0
         do jj=nys(plane)-2,nys(plane)-1
            do ii=nxs(plane)-1,nxe(plane)+1
               vhi(:,ii,jj,plane)=tmp(cnt+1:cnt+nY)
               cnt=cnt+nY
               w(  :,ii,jj,plane)=tmp(cnt+1:cnt+dimw)
               cnt=cnt+dimw
            end do
         end do
      end if
   end do

   do plane = nps,npe
      if(gy(plane) .ne. ngy(plane)) call MPI_Wait(ireq(plane),istatus,ierrs(plane))
   end do

   !grid north to south
   !send
   do plane = nps,npe
      if(gy(plane) .ne. 1) then
         cnts(plane)=0
         do jj=nys(plane),nys(plane)+1
            do ii=nxs(plane)-1,nxe(plane)+1
               tmps(cnts(plane)+1:cnts(plane)+nY,  plane) =vhi(:,ii,jj,plane)
               cnts(plane)=cnts(plane)+nY
               tmps(cnts(plane)+1:cnts(plane)+dimw,plane) =w(  :,ii,jj,plane)
               cnts(plane)=cnts(plane)+dimw
            end do
         end do
         flag(plane)=plane
         call MPI_Isend(tmps(:,plane),cnts(plane),MPI_DOUBLE_PRECISION,gs(plane),flag(plane),&
                         MPI_COMM_WORLD,ireq(plane),ierrs(plane))
      end if
   end do

   !receive
   do plane = nps,npe
      if(gy(plane) .ne.  ngy(plane)) then
         cnt=DLength*2*(bwx(plane)+2)
         call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gn(plane),plane,&
                         MPI_COMM_WORLD, istatus, ierr)
      
         cnt=0
         do jj=nye(plane)+1,nye(plane)+2
            do ii=nxs(plane)-1,nxe(plane)+1
               vhi(:,ii,jj,plane)=tmp(cnt+1:cnt+nY)
               cnt=cnt+nY
               w(  :,ii,jj,plane)=tmp(cnt+1:cnt+dimw)
               cnt=cnt+dimw
            end do
         end do
      end if
   end do
   
   do plane = nps,npe
      if(gy(plane) .ne. 1) call MPI_Wait(ireq(plane),istatus,ierrs(plane))
   end do
end subroutine MPI_COMMUNICATIONS_J_DIRECTION!}}}
subroutine set_corners!{{{
   implicit none
   integer i,j,plane
   integer pi,pj

   do plane = nps,npe
      do pi = -1,1,2
         if(pi <0) then
            i = nxs(plane)
         else
            i = nxe(plane)
         end if
         do pj = -1,1,2
            if(pj <0) then
               j = nys(plane)
            else
               j = nye(plane)
            end if
            w(  :,i+pi,j+pj,plane) = &
                     1d0/3d0*(w(  :,i,j,plane)+w(  :,i,j+pj,plane)+w(  :,i+pi,j,plane))
            vhi(:,i+pi,j+pj,plane) = &
                     1d0/3d0*(vhi(:,i,j,plane)+vhi(:,i,j+pj,plane)+vhi(:,i+pi,j,plane))
         end do
      end do
   end do
end subroutine set_corners!}}}
end subroutine set_BC

