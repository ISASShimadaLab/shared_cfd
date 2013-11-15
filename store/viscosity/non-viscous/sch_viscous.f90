subroutine init_TGv
   use mod_mpi
   use variable
   implicit none
   integer i,j,plane

   do plane = nps,npe
      !set TGvi
      !$omp parallel do default(shared),private(i)
      do j=nys(plane),nye(plane)
         do i=nxs(plane)-1,nxe(plane)
            TGvi(:,i,j,plane)=0d0
         end do
      end do
      !$omp end parallel do

      !set TGvj
      !$omp parallel do default(shared),private(i)
      do j=nys(plane)-1,nye(plane)
         do i=nxs(plane),nxe(plane)
            TGvj(:,i,j,plane)=0d0
         end do
      end do
      !$omp end parallel do
   end do
end subroutine init_TGv

subroutine set_TGv
end subroutine set_TGv
