subroutine calc_next_step_implicit
   use mod_mpi
   use variable
   use var_lusgs
   implicit none
   integer i,j,sm
   double precision Dq(1:dimq,nxs-1:nxe+1,nys-1:nye+1)

   !set RHS
   !$omp parallel do private(i)
   do j=nys,nye
      do i=nxs,nxe
         Dq(:,i,j)= dsi(i-1,  j)*(TGi(1:dimq,i-1,j  )-TGvi(1:dimq,i-1,j  ))&
                   -dsi(  i,  j)*(TGi(1:dimq,i  ,j  )-TGvi(1:dimq,i  ,j  ))&
                   +dsj(  i,j-1)*(TGj(1:dimq,i  ,j-1)-TGvj(1:dimq,i  ,j-1))&
                   -dsj(  i,  j)*(TGj(1:dimq,i  ,j  )-TGvj(1:dimq,i  ,j  ))
         Dq(:,i,j)= Dq(:,i,j)+pre1(:,i,j)*dot_product(pre2(:,i,j),Dq(:,i,j))
      end do
   end do
   !$omp end parallel do

   !set BC Delta q
   Dq(:,nxs-1,:)=0d0
   Dq(:,:,nys-1)=0d0

   !calculate forward
   do sm=nxs+nys,nxe+nye
      !$omp parallel do private(j)
      do i=max(nxs,sm-nye),min(sm-nys,nxe)
         j=sm-i
         Dq(:,i,j)=alpha(i,j)*(Dq(:,i,j)&
                              +dsci(i,j)*matmul(Ap(:,:,i,j),Dq(:,i-1,j  ))&
                              +dscj(i,j)*matmul(Bp(:,:,i,j),Dq(:,i  ,j-1)))
      end do
      !$omp end parallel do
   end do

   !set BC Delta q
   Dq(:,nxe+1,:)=0d0
   Dq(:,:,nye+1)=0d0

   !calculate backward
   do sm=nxe+nye,nxs+nys,-1
      !$omp parallel do private(j)
      do i=max(nxs,sm-nye),min(sm-nys,nxe)
         j=sm-i
         Dq(:,i,j)=Dq(:,i,j)&
                  -alpha(i,j)*(dsci(i,j)*matmul(Am(:,:,i,j),Dq(:,i+1,j  ))&
                              +dscj(i,j)*matmul(Bm(:,:,i,j),Dq(:,i  ,j+1)))
      end do
      !$omp end parallel do
   end do

   !$omp parallel do private(i)
   !add Dq
   do j=nys,nye
      do i=nxs,nxe
         q(:,i,j)=q(:,i,j)+Dq(:,i,j)
      end do
   end do
   !$omp end parallel do
end subroutine calc_next_step_implicit

