subroutine set_dt(step)
   use mod_mpi
   use grbl_prmtr
   use variable
   implicit none
   integer,intent(in)::step
   double precision rho,u,v,p,un,a
   double precision uds,tmp
   integer i,j,plane

   tmp=1d300
   do plane = nps,npe
      !$omp parallel do private(i,rho,u,v,p,un,a,uds) reduction(min:tmp)
      !!$omp parallel do private(i,rho,u,v,p,un,a,uds)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            rho=w(1,i,j,plane)
            u  =w(2,i,j,plane)
            v  =w(3,i,j,plane)
            p  =w(4,i,j,plane)
            a=sqrt(w(indxg,i,j,plane)*p/rho)
            uds=0d0

            un=vni(1,i  ,j,plane)*u+vni(2,i  ,j,plane)*v
            uds=uds+(abs(un)+a)      *dsi(i  ,j,plane)
            un=vni(1,i-1,j,plane)*u+vni(2,i-1,j,plane)*v
            uds=uds+(abs(un)+a)      *dsi(i-1,j,plane)
            un=vnj(1,i,j  ,plane)*u+vnj(2,i,j  ,plane)*v
            uds=uds+(abs(un)+a)      *dsj(i,j  ,plane)
            un=vnj(1,i,j-1,plane)*u+vnj(2,i,j-1,plane)*v
            uds=uds+(abs(un)+a)      *dsj(i,j-1,plane)

            dt_mat(i,j,plane)=cfl*Vol(i,j,plane)/max(uds,w(indxMu,i,j,plane)/rho)
            if(tmp>dt_mat(i,j,plane)) tmp=dt_mat(i,j,plane)
         end do
      end do
      !$omp end parallel do
   end do

   call MPI_Allreduce(tmp, dt_grbl, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

   tt = tt+dt_grbl
end subroutine set_dt
