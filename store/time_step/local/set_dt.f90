subroutine set_dt(step)
   use mod_mpi
   use grbl_prmtr
   use variable
   implicit none
   integer,intent(in)::step
   double precision rho,u,v,p,un,a
   double precision uds
   integer i,j

   !$omp parallel do private(i,rho,u,v,p,un,a,uds)
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
      end do
   end do
   !$omp end parallel do
end subroutine set_dt
