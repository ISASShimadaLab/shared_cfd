module gas
   double precision,parameter::kappa_gas = 1.4
   double precision,parameter::R_gas     = 287
   double precision,parameter::nu_gas    = 1.6e-5
end module gas

subroutine set_thermo_prop
   use mod_mpi
   use variable
   use gas
   implicit none
   double precision ei,T
   integer i,j,plane

   do plane = nps,npe
      !$omp parallel do private(i,ei,T)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            ei = q(nY+3,i,j,plane)/w(1,i,j,plane)-0.5d0*(w(2,i,j,plane)**2+w(3,i,j,plane)**2)
            T  = ei*(kappa_gas-1d0)/R_gas

            w(4,     i,j,plane) = w(1,i,j,plane)*R_gas*T
            w(indxg, i,j,plane) = kappa_gas
            w(indxht,i,j,plane) = (q(nY+3,i,j,plane)+w(4,i,j,plane))/w(1,i,j,plane)
            w(indxMu,i,j,plane) = w(1,i,j,plane)*nu_gas
            w(indxR, i,j,plane) = R_gas
            DHi(1,   i,j,plane) = 0d0
            vhi(1,   i,j,plane) = kappa_gas/(kappa_gas-1d0)*R_gas*T
         end do
      end do
      !$omp end parallel do
   end do
end subroutine set_thermo_prop
