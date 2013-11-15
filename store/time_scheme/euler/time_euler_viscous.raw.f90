subroutine calc_next_step
   use mod_mpi
   use variable
   implicit none
   integer i,j
   double precision vol_dqdt(1:dimq)

   !$omp parallel do private(i,vol_dqdt)
   do j=nys,nye
      do i=nxs,nxe
         vol_dqdt=+dsi(i-1,  j)*(TGi(1:dimq,i-1,j  )-TGvi(1:dimq,i-1,j  ))&
                  -dsi(  i,  j)*(TGi(1:dimq,i  ,j  )-TGvi(1:dimq,i  ,j  ))&
                  +dsj(  i,j-1)*(TGj(1:dimq,i  ,j-1)-TGvj(1:dimq,i  ,j-1))&
                  -dsj(  i,  j)*(TGj(1:dimq,i  ,j  )-TGvj(1:dimq,i  ,j  ))

         q(:,i,j)=qp(:,i,j)+DT_LOCAL_GLOBAL/Vol(i,j)*vol_dqdt
      end do
   end do
   !$omp end parallel do
end subroutine calc_next_step
