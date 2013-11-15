subroutine calc_next_step_implicit(step_internal)
   use grbl_prmtr
   use mod_mpi
   use variable
   use var_lusgs
   implicit none
   integer,intent(in)::step_internal
   integer i,j,k,sm
   double precision Dq(1:dimq,nxs-1:nxe+1,nys-1:nye+1)
   double precision omega,tmp

   !set RHS
   !$omp parallel do private(i)
   do j=nys,nye
      do i=nxs,nxe
         Dq(:,i,j)= dsi(i-1,  j)*(TGi(1:dimq,i-1,j  )-TGvi(1:dimq,i-1,j  ))&
                   -dsi(  i,  j)*(TGi(1:dimq,i  ,j  )-TGvi(1:dimq,i  ,j  ))&
                   +dsj(  i,j-1)*(TGj(1:dimq,i  ,j-1)-TGvj(1:dimq,i  ,j-1))&
                   -dsj(  i,  j)*(TGj(1:dimq,i  ,j  )-TGvj(1:dimq,i  ,j  ))
         Dq(:,i,j)= Dq(:,i,j)-Vol(i,j)*dqdt(:,i,j)

         Dq(:,i,j)=Dq(:,i,j)+pre_sca(2,i,j)*dot_product(pre2(:,i,j),Dq(:,i,j))*pre1(:,i,j)
      end do
   end do
   !$omp end parallel do

   !!!!!!!!!!!!! FIRST STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !set BC Delta q
   Dq(:,nxs-1,:)=0d0
   Dq(:,:,nys-1)=0d0

   !calculate forward
   do sm=nxs+nys,nxe+nye
      !$omp parallel do private(j)
      do i=max(nxs,sm-nye),min(sm-nys,nxe)
         j=sm-i
         Dq(:,i,j)=alpha(i,j)*(Dq(:,i,j)&
                              +dsci(i-1,j  )*matmul(Ap(:,:,i-1,j  ),Dq(:,i-1,j  ))&
                              +dscj(i  ,j-1)*matmul(Bp(:,:,i  ,j-1),Dq(:,i  ,j-1)))

         Dq(:,i,j)=Dq(:,i,j)+pre_sca(4,i,j) *dot_product(pre2(:,i,j),Dq(:,i,j)) *pre1(:,i,j)
         Dq(:,i,j)=1d0/pre_sca(1,i,j)*Dq(:,i,j)
      end do
      !$omp end parallel do
   end do

   !!!!!!!!!!!!! SECOND STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !$omp parallel do private(i)
   do j=nys,nye
      do i=nxs,nxe
         Dq(:,i,j)=Dq(:,i,j)+pre_sca(3,i,j) *dot_product(pre2(:,i,j),Dq(:,i,j)) *pre1(:,i,j)
         Dq(:,i,j)=pre_sca(1,i,j)*Dq(:,i,j)
      end do
   end do
   !$omp end parallel do

   !!!!!!!!!!!!! THIRD STEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !set BC Delta q
   Dq(:,nxe+1,:)=0d0
   Dq(:,:,nye+1)=0d0

   omega=1d300
   res = 0d0
   !calculate backward
   do sm=nxe+nye,nxs+nys,-1
      !$omp parallel do private(j,k) reduction(min:omega) reduction(+:res)
      do i=max(nxs,sm-nye),min(sm-nys,nxe)
         j=sm-i
         Dq(:,i,j)=Dq(:,i,j)&
                  -alpha(i,j)*(dsci(i+1,j  )*matmul(Am(:,:,i+1,j  ),Dq(:,i+1,j  ))&
                              +dscj(i  ,j+1)*matmul(Bm(:,:,i  ,j+1),Dq(:,i  ,j+1)))

         Dq(:,i,j)=Dq(:,i,j)+pre_sca(4,i,j) *dot_product(pre2(:,i,j),Dq(:,i,j)) *pre1(:,i,j)
         Dq(:,i,j)=1d0/pre_sca(1,i,j)*Dq(:,i,j)

         res =res +Dq(nY+3,i,j)**2
         !do k=1,dimq
         !   res =res +Dq(k,i,j)**2
         !end do

         do k=1,nY
            if(Dq(k,i,j)<0d0) then
               omega = min(omega,-q(k,i,j)/Dq(k,i,j)*Dqmax)
               if(omega<omega_min) then
                  print *,"Error: Omega becomes too small. omega = ",omega
                  call exit(1)
               end if
            end if
         end do
      end do
      !$omp end parallel do
   end do

   !!!!!!!!! FOR PARAMETERS ADJUSTMENTS !!!!!!!!!!!!!!!!
   tmp = omega
   call MPI_Reduce(tmp,omega,1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
   tmp = res
   call MPI_Reduce(tmp,  res,1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   if(myid .eq. 0) then
      if(res .ne. res) then
         print *,"res is NaN at internal step=",step_internal
         call exit(1)
      end if
      if(step_internal .eq. 1) res1 = res
      write(66,'(i3,100e12.3)') step_internal,omega,res/res1,dt_grbl
   end if
   !!!!!!!!! END OF FOR PARAMETERS ADJUSTMENTS !!!!!!!!!

   omega      = min(omega_max,omega)

   !add Dq
   !$omp parallel do private(i)
   do j=nys,nye
      do i=nxs,nxe
         q(:,i,j)=q(:,i,j)+omega*Dq(:,i,j)
      end do
   end do
   !$omp end parallel do
end subroutine calc_next_step_implicit

