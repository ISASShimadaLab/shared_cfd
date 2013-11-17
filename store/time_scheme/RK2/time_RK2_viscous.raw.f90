subroutine calc_next_step_RK2_first
   use grbl_prmtr
   use mod_mpi
   use variable
   implicit none
   integer i,j,plane
   double precision vol_dqdt(1:dimq)

   do plane=nps,npe
      !$omp parallel do private(i,vol_dqdt)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            vol_dqdt=+dsi(i-1,  j,plane)*(TGi(1:dimq,i-1,j  ,plane)-TGvi(1:dimq,i-1,j  ,plane))&
                     -dsi(  i,  j,plane)*(TGi(1:dimq,i  ,j  ,plane)-TGvi(1:dimq,i  ,j  ,plane))&
                     +dsj(  i,j-1,plane)*(TGj(1:dimq,i  ,j-1,plane)-TGvj(1:dimq,i  ,j-1,plane))&
                     -dsj(  i,  j,plane)*(TGj(1:dimq,i  ,j  ,plane)-TGvj(1:dimq,i  ,j  ,plane))
      
            q(:,i,j,plane)=qp(:,i,j,plane)+0.5d0*DT_LOCAL_GLOBAL/Vol(i,j,plane)*vol_dqdt
         end do
      end do
      !$omp end parallel do
   end do
end subroutine calc_next_step_RK2_first

subroutine calc_next_step_RK2_second
   use grbl_prmtr
   use mod_mpi
   use variable
   implicit none
   integer i,j,k,plane
   double precision vol_dqdt(1:dimq)

   do plane=nps,npe
      !$omp parallel do private(i,k,vol_dqdt)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            vol_dqdt=+dsi(i-1,  j,plane)*(TGi(1:dimq,i-1,j  ,plane)-TGvi(1:dimq,i-1,j  ,plane))&
                     -dsi(  i,  j,plane)*(TGi(1:dimq,i  ,j  ,plane)-TGvi(1:dimq,i  ,j  ,plane))&
                     +dsj(  i,j-1,plane)*(TGj(1:dimq,i  ,j-1,plane)-TGvj(1:dimq,i  ,j-1,plane))&
                     -dsj(  i,  j,plane)*(TGj(1:dimq,i  ,j  ,plane)-TGvj(1:dimq,i  ,j  ,plane))
      
            q(:,i,j,plane)=qp(:,i,j,plane)+DT_LOCAL_GLOBAL/Vol(i,j,plane)*vol_dqdt
            do k=1,nY
               if(q(k,i,j,plane)<0d0) then
                  print *,"Error:Negative partial density"
                  print '(a)',"(i,j)="
                  print '(2i3)',i,j
                  print '(a)',"q="
                  print '(100es11.3)',q(:,i,j,plane)
                  call exit(1)
               end if
            end do
         end do
      end do
      !$omp end parallel do
   end do
end subroutine calc_next_step_RK2_second

