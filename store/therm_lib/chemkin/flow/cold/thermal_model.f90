!data consistency check between flow and chemistry
subroutine check_chem_and_flow!{{{
   use grbl_prmtr
   use const_chem
   implicit none
   if(nY .ne. ns) then
      print *,"Error. nY and ns is different. (nY,ns)=",nY,ns
      stop
   end if
end subroutine check_chem_and_flow!}}}

! flame sheet model
subroutine set_thermo_prop!{{{
   use mod_mpi
   use prmtr
   use variable
   use chem_var
   implicit none
   double precision ei,T,MW,kappa,xi
   integer i,j,plane

   do plane = nps,npe
      !$omp parallel do private(i,ei,T,MW,kappa,xi)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            !values necessary to calculate thermo values
            ei=q(nY+3,i,j,plane)/w(1,i,j,plane)-0.5d0*(w(2,i,j,plane)**2+w(3,i,j,plane)**2)
            T=w(4,i,j,plane)/(w(1,i,j,plane)*w(indxR,i,j,plane))

            !calc therm
            call calc_T(q(1:nY,i,j,plane),w(1,i,j,plane),ei, T, kappa,MW,DHi(:,i,j,plane),vhi(:,i,j,plane),w(indxMu,i,j,plane))

            !calculate R_gas and w
            w(indxR, i,j,plane)=R_uni/MW
            w(     4,i,j,plane)=w(1,i,j,plane)*w(indxR,i,j,plane)*T
            w(indxg, i,j,plane)=kappa
            w(indxht,i,j,plane)=(q(nY+3,i,j,plane)+w(4,i,j,plane))/w(1,i,j,plane)
         end do
      end do
      !$omp end parallel do
   end do
end subroutine set_thermo_prop!}}}

!boundary condition
subroutine calc_boundary(p,T,deg, wt,vhi)!{{{
   use grbl_prmtr
   use prmtr
   use const_chem
   implicit none
   double precision,intent(in)::p
   double precision,intent(inout)::T
   double precision,intent(in)::deg

   double precision,intent(out)::wt(dimw)
   double precision,intent(out)::vhi(ns)

   double precision,dimension(ns)::vrho,DHi
   double precision rho,E,MWave
   integer j

   call calc_vrho(p,T,deg, rho,vrho,E)
   call calc_T(vrho,rho,E, T, wt(indxg),MWave,DHi,vhi,wt(indxMu))

   !set wt
   wt(1)=rho
   wt(2)=0d0
   wt(3)=0d0
   wt(4)=p
   do j=1,ns
      wt(4+j)=vrho(j)/rho
   end do
   wt(indxht)=E+p/rho
   wt(indxR) =R_uni/MWave
end subroutine calc_boundary!}}}
