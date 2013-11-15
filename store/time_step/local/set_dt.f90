subroutine set_dt(step)
   use mod_mpi
   use grbl_prmtr
   use variable
   implicit none
   integer,intent(in)::step
   double precision rho,u,v,p,un,a
   double precision uds
   integer i,j,plane

   do plane = nps,npe
      !$omp parallel do private(i,rho,u,v,p,un,a,uds)
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
         end do
      end do
      !$omp end parallel do
   end do
end subroutine set_dt
