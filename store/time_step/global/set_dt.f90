subroutine set_dt(step)
   use mod_mpi
   use grbl_prmtr
   use variable
   implicit none
   integer,intent(in)::step
   double precision rho,u,v,p,un,a
   double precision uds,tmp
   integer i,j

   tmp=1d300
   !$omp parallel do private(i,rho,u,v,p,un,a,uds) reduction(min:tmp)
   !!$omp parallel do private(i,rho,u,v,p,un,a,uds)
   do j=nys,nye
      do i=nxs,nxe
         rho=w(1,i,j)
         u  =w(2,i,j)
         v  =w(3,i,j)
         p  =w(4,i,j)
         a=sqrt(w(indxg,i,j)*p/rho)
         uds=0d0

         un=vni(1,i  ,j)*u+vni(2,i  ,j)*v
         uds=uds+(abs(un)+a)*dsi(i  ,j)
         un=vni(1,i-1,j)*u+vni(2,i-1,j)*v
         uds=uds+(abs(un)+a)*dsi(i-1,j)
         un=vnj(1,i,j  )*u+vnj(2,i,j  )*v
         uds=uds+(abs(un)+a)*dsj(i,j  )
         un=vnj(1,i,j-1)*u+vnj(2,i,j-1)*v
         uds=uds+(abs(un)+a)*dsj(i,j-1)

         dt_mat(i,j)=cfl*Vol(i,j)/max(uds,w(indxMu,i,j)/rho)
         if(tmp>dt_mat(i,j)) tmp=dt_mat(i,j)
      end do
   end do
   !$omp end parallel do

   call MPI_Allreduce(tmp, dt_grbl, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

   tt = tt+dt_grbl
end subroutine set_dt
