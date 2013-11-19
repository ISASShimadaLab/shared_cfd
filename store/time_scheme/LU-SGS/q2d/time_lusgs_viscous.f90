subroutine calc_next_step_implicit
   use mod_mpi
   use variable
   use var_lusgs
   implicit none
   integer i,j,l,sm,plane
   double precision Dq(1:dimq,0:nimax+1,0:njmax+1,Nplane)
   double precision Dq_small(1:dimq)
   double precision tmp

   do plane = nps,npe
      !set RHS
      !$omp parallel do private(i)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            Dq(:,i,j,plane)= dsi(i-1,  j,plane)*(TGi(1:dimq,i-1,j  ,plane)-TGvi(1:dimq,i-1,j  ,plane))&
                            -dsi(  i,  j,plane)*(TGi(1:dimq,i  ,j  ,plane)-TGvi(1:dimq,i  ,j  ,plane))&
                            +dsj(  i,j-1,plane)*(TGj(1:dimq,i  ,j-1,plane)-TGvj(1:dimq,i  ,j-1,plane))&
                            -dsj(  i,  j,plane)*(TGj(1:dimq,i  ,j  ,plane)-TGvj(1:dimq,i  ,j  ,plane))&
                            +Area( i,  j,plane)*(Sq( 1:dimq,i  ,j  ,plane)+Svq( 1:dimq,i  ,j  ,plane))
         end do
      end do
      !$omp end parallel do

      !set BC Delta q
      Dq(:,nxs(plane)-1,nys(plane)-1:nye(plane),plane)=0d0
      Dq(:,nxs(plane)-1:nxe(plane),nys(plane)-1,plane)=0d0

      !calculate forward
      do sm=nxs(plane)+nys(plane),nxe(plane)+nye(plane)
         !$omp parallel do private(j,l,tmp,Dq_small)
         do i=max(nxs(plane),sm-nye(plane)),min(sm-nys(plane),nxe(plane))
            j=sm-i
            Dq_small = dsci(i-1,j  ,plane)*matmul(Ap(:,:,i-1,j  ,plane),Dq(:,i-1,j  ,plane))&
                      +dscj(i  ,j-1,plane)*matmul(Bp(:,:,i  ,j-1,plane),Dq(:,i  ,j-1,plane))
            Dq(:,i,j,plane)=alpha(i,j,plane)*(Dq(:,i,j,plane)+Dq_small)

            tmp=0d0
            do l=1,nY+1
               tmp = tmp + Dq(l,i,j,plane)*dpdq(l,i,j,plane)
            end do
            tmp = tmp + Dq(nY+3,i,j,plane)*dpdq(nY+3,i,j,plane)
            tmp = Dq(nY+2,i,j,plane) + alpha(i,j,plane)*Area(i,j,plane)*tmp
            Dq(nY+2,i,j,plane) = tmp/(1d0-alpha(i,j,plane)*Area(i,j,plane)*dpdq(nY+2,i,j,plane))
         end do
         !$omp end parallel do
      end do

      !set BC Delta q
      Dq(:,nxe(plane)+1,nys(plane):nye(plane)+1,plane)=0d0
      Dq(:,nxs(plane):nxe(plane)+1,nye(plane)+1,plane)=0d0

      !calculate backward
      do sm=nxe(plane)+nye(plane),nxs(plane)+nys(plane),-1
         !$omp parallel do private(j)
         do i=max(nxs(plane),sm-nye(plane)),min(sm-nys(plane),nxe(plane))
            j=sm-i
            Dq(:,i,j,plane)=Dq(:,i,j,plane)&
                           -alpha(i,j,plane)*(dsci(i+1,j  ,plane)*matmul(Am(:,:,i+1,j  ,plane),Dq(:,i+1,j  ,plane))&
                                             +dscj(i  ,j+1,plane)*matmul(Bm(:,:,i  ,j+1,plane),Dq(:,i  ,j+1,plane)))
         end do
         !$omp end parallel do
      end do

      !$omp parallel do private(i)
      !add Dq
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            q(:,i,j,plane)=q(:,i,j,plane)+Dq(:,i,j,plane)
         end do
      end do
      !$omp end parallel do
   end do
end subroutine calc_next_step_implicit

