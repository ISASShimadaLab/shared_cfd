subroutine set_IC
   use grbl_prmtr
   use prmtr
   use variable
   use mod_mpi
   implicit none
   integer i,j,plane

   do plane=nps,npe
      do i=nxs(plane),nxe(plane)
         do j=nys(plane),nye(plane)
            !q(:,    i,j,plane) =
            !w(4,    i,j,plane) =
            !w(indxg,i,j,plane) =
            !w(indxR,i,j,plane) =
         end do
      end do
   end do
end subroutine set_IC

subroutine set_BC(step)
   use grbl_prmtr
   use prmtr
   use variable
   use mod_mpi
   implicit none
   integer,intent(in)::step
   integer i,j,plane

   integer,parameter::DLength=dimw+2*nV !for MPI Communication
