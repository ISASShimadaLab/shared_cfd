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
            call cold_flow(w(5:6,i,j,plane),ei, T, MWave,kappa,mu,DHi(:,i,j,plane),vhi(:,i,j,plane))

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
