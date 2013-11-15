module gas
   double precision,parameter::kappa_gas = KAPPA
   double precision,parameter::R_gas     = RGAS
   double precision,parameter::nu_gas    = NU
end module gas

subroutine set_thermo_prop
   use mod_mpi
   use variable
   use gas
   implicit none
   double precision ei,T
   integer i,j

   !$omp parallel do private(i,ei,T)
   do j=nys,nye
      do i=nxs,nxe
         ei = q(nY+3,i,j)/w(1,i,j)-0.5d0*(w(2,i,j)**2+w(3,i,j)**2)
         T  = ei*(kappa_gas-1d0)/R_gas

         w(4,     i,j) = w(1,i,j)*R_gas*T
         w(indxg, i,j) = kappa_gas
         w(indxht,i,j) = (q(nY+3,i,j)+w(4,i,j))/w(1,i,j)
         w(indxMu,i,j) = w(1,i,j)*nu_gas
         w(indxR, i,j) = R_gas
         DHi(1,   i,j) = 0d0
         vhi(1,   i,j) = kappa_gas/(kappa_gas-1d0)*R_gas*T
      end do
   end do
   !$omp end parallel do
end subroutine set_thermo_prop
