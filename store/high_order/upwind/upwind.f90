module half_value
   use mod_mpi
   use grbl_prmtr
   implicit none
contains
subroutine set_HV(q,qHli,qHri,qHlj,qHrj)
   implicit none
   double precision,intent(in)::    q(1:dimw,-1:nimax+2,-1:njmax+2,Nplane)
   double precision,intent(out)::qHli(1:dimw, 0:nimax,   1:njmax,  Nplane)
   double precision,intent(out)::qHri(1:dimw, 0:nimax,   1:njmax,  Nplane)
   double precision,intent(out)::qHlj(1:dimw, 1:nimax,   0:njmax,  Nplane)
   double precision,intent(out)::qHrj(1:dimw, 1:nimax,   0:njmax,  Nplane)

   integer i,j,plane

  
   do plane = nps,npe
      !$omp parallel do private(i)
      do j=nys(plane),nye(plane)
         do i=nxs(plane)-1,nxe(plane)
            qHli(:,i,j,plane)=q(1:dimw,i  ,j,plane)
            qHri(:,i,j,plane)=q(1:dimw,i+1,j,plane)
         end do
      end do
      !$omp end parallel do

      !$omp parallel do private(i)
      do j=nys(plane)-1,nye(plane)
         do i=nxs(plane),nxe(plane)
            qHlj(:,i,j,plane)=q(1:dimw,i,j  ,plane)
            qHrj(:,i,j,plane)=q(1:dimw,i,j+1,plane)
         end do
      end do
      !$omp end parallel do
   end do
end subroutine set_HV
end module half_value

