subroutine init_TGv
   use mod_mpi
   use variable
   implicit none
   integer i,j

   !set TGvi
   !$omp parallel do default(shared),private(i)
   do j=nys,nye
      do i=nxs-1,nxe
         TGvi(:,i,j)=0d0
      end do
   end do
   !$omp end parallel do

   !set TGvj
   !$omp parallel do default(shared),private(i)
   do j=nys-1,nye
      do i=nxs,nxe
         TGvj(:,i,j)=0d0
      end do
   end do
   !$omp end parallel do
end subroutine init_TGv

subroutine set_TGv
end subroutine set_TGv
