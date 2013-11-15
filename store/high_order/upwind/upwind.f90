module half_value
   use mod_mpi
   use grbl_prmtr
   implicit none
contains
subroutine set_HV(q,qHli,qHri,qHlj,qHrj)
   implicit none
   double precision,intent(in)::    q(1:dimw,-1:ni+2,-1:nj+2)
   double precision,intent(out)::qHli(1:dimw, 0:ni,   1:nj)
   double precision,intent(out)::qHri(1:dimw, 0:ni,   1:nj)
   double precision,intent(out)::qHlj(1:dimw, 1:ni,   0:nj)
   double precision,intent(out)::qHrj(1:dimw, 1:ni,   0:nj)

   integer i,j

   !$omp parallel do default(none) private(i) shared(j,q,qHli,qHri,nxs,nxe,nys,nye)
   do j=nys,nye
      do i=nxs-1,nxe
         qHli(:,i,j)=q(1:dimw,i  ,j)
         qHri(:,i,j)=q(1:dimw,i+1,j)
      end do
   end do
   !$omp end parallel do

   !$omp parallel do default(none) private(i) shared(j,q,qHlj,qHrj,nxs,nxe,nys,nye)
   do j=nys-1,nye
      do i=nxs,nxe
         qHlj(:,i,j)=q(1:dimw,i,j  )
         qHrj(:,i,j)=q(1:dimw,i,j+1)
      end do
   end do
   !$omp end parallel do
end subroutine set_HV
end module half_value

