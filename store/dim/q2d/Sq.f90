subroutine set_Sq
   use mod_mpi
   use grbl_prmtr
   use variable
   implicit none
   integer i,j,plane

   do plane = nps,npe
      !$omp parallel do private(i)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            Sq(:,i,j,plane)=0d0
            Sq(nY+2,i,j,plane)=w(4,i,j,plane)
         end do
      end do
      !$omp end parallel do
   end do
end subroutine set_Sq

