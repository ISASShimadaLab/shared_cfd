subroutine calc_next_step
   use mod_mpi
   use variable
   implicit none
   integer i,j,plane
   double precision vol_dqdt(1:dimq)

   do plane = nps,npe
      !$omp parallel do private(i,vol_dqdt)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            vol_dqdt=+dsi(i-1,  j,plane)*(TGi(1:dimq,i-1,j  ,plane)-TGvi(1:dimq,i-1,j  ,plane))&
                     -dsi(  i,  j,plane)*(TGi(1:dimq,i  ,j  ,plane)-TGvi(1:dimq,i  ,j  ,plane))&
                     +dsj(  i,j-1,plane)*(TGj(1:dimq,i  ,j-1,plane)-TGvj(1:dimq,i  ,j-1,plane))&
                     -dsj(  i,  j,plane)*(TGj(1:dimq,i  ,j  ,plane)-TGvj(1:dimq,i  ,j  ,plane))&
                     +Area( i,  j,plane)*(Sq( 1:dimq,i  ,j  ,plane)+Svq( 1:dimq,i  ,j  ,plane))

            q(:,i,j,plane)=qp(:,i,j,plane)+DT_LOCAL_GLOBAL/Vol(i,j,plane)*vol_dqdt
         end do
      end do
      !$omp end parallel do
   end do
end subroutine calc_next_step
