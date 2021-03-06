subroutine set_thermo_prop
   use mod_mpi
   use variable
   use prmtr
   implicit none
   double precision ei,T,MWave,kappa,mu
   integer i,j,plane

   do plane = nps,npe
      !$omp parallel do private(i,ei,T,MWave,kappa,mu)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            ei = q(nY+3,i,j,plane)/w(1,i,j,plane)-0.5d0*(w(2,i,j,plane)**2+w(3,i,j,plane)**2)
            T  = w(4,i,j,plane)/(w(1,i,j,plane)*w(indxR,i,j,plane))
            call flame_sheet(w(5:6,i,j,plane),ei, T, MWave,kappa,mu,DHi(:,i,j,plane),Yv(:,i,j,plane),vhi(:,i,j,plane))

            DHi(:,   i,j,plane) = 0d0
            w(4,     i,j,plane) = w(1,i,j,plane)*(R_uni/MWave)*T
            w(indxg, i,j,plane) = kappa
            w(indxht,i,j,plane) = (q(nY+3,i,j,plane)+w(4,i,j,plane))/w(1,i,j,plane)
            w(indxMu,i,j,plane) = mu
            w(indxR, i,j,plane) = R_uni/MWave
         end do
      end do
      !$omp end parallel do
   end do
end subroutine set_thermo_prop

subroutine YPT2w(Yf,p,T,wt,vhit)
   use grbl_prmtr
   use chem
   use prmtr
   implicit none
   double precision,intent(in) ::Yf
   double precision,intent(in) ::p
   double precision,intent(in) ::T
   double precision,intent(out)::wt(dimw)
   double precision,intent(out)::vhit(nV)

   double precision Y(2),E,Yv(3)
   double precision MWave,kappa,mu,rho


   Y(1)=Yf
   Y(2)=1d0-Yf

   call flame_sheet_PTgiven(Y,T, E,MWave,kappa,mu,Yv,vhit)
   rho = P/(R_uni/MWave*T)
   wt(1)=rho
   wt(2)=0d0
   wt(3)=0d0
   wt(4)=p
   wt(5)=Y(1)
   wt(6)=Y(2)
   wt(indxg ) = kappa
   wt(indxht) = E+R_uni/MWave*T
   wt(indxR ) =   R_uni/MWave
   wt(indxMu) = mu
end subroutine YPT2w

