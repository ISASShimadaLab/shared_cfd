module half_value
   use mod_mpi
   use grbl_prmtr
   implicit none
   double precision phi_mus,k_mus
contains
subroutine set_HV(q,qHli,qHri,qHlj,qHrj)!van albada
   implicit none
   double precision,intent(in)::    q(1:dimw,-1:ni+2,-1:nj+2)
   double precision,intent(out)::qHli(1:dimw, 0:ni,   1:nj)
   double precision,intent(out)::qHri(1:dimw, 0:ni,   1:nj)
   double precision,intent(out)::qHlj(1:dimw, 1:ni,   0:nj)
   double precision,intent(out)::qHrj(1:dimw, 1:ni,   0:nj)

   double precision          delta_ui(1:dimw,0:ni+1,  1:nj)
   double precision          delta_li(1:dimw,0:ni+1,  1:nj)
   double precision          delta_uj(1:dimw,  1:ni,0:nj+1)
   double precision          delta_lj(1:dimw,  1:ni,0:nj+1)
   double precision,dimension(dimw)::Du,Dl,s
   double precision,parameter::epsln=1d-300
   !double precision,parameter::epsln=1d-6
   integer i,j
   
   !$omp parallel do default(none) private(i,Du,Dl,s) shared(j,q,delta_ui,delta_li,nxs,nxe,nys,nye,phi_mus,k_mus)
   do j=nys,nye
      do i=nxs-1,nxe+1
         Du(1:dimw)=q(:,i+1,j)-q(:,i  ,j)
         Dl(1:dimw)=q(:,i  ,j)-q(:,i-1,j)

         !s=(2d0*Du*Dl+epsln)/(Du**2+Dl**2+epsln)
         s=(2d0*Du*Dl)/(Du**2+Dl**2+epsln)

         delta_ui(:,i,j)=0.5d0*s*((1d0+k_mus*s)*Du+(1d0-k_mus*s)*Dl)
         delta_li(:,i,j)=0.5d0*s*((1d0+k_mus*s)*Dl+(1d0-k_mus*s)*Du)
      end do
   end do
   !$omp end parallel do

   !$omp parallel do default(none) private(i) shared(j,q,delta_ui,delta_li,qHli,qHri,nxs,nxe,nys,nye,phi_mus,k_mus)
   do j=nys,nye
      do i=nxs-1,nxe
         qHli(:,i,j)=q(1:dimw,i  ,j)+phi_mus*0.5d0*delta_ui(1:dimw,i  ,j)
         qHri(:,i,j)=q(1:dimw,i+1,j)-phi_mus*0.5d0*delta_li(1:dimw,i+1,j)
      end do
   end do
   !$omp end parallel do

   !$omp parallel do default(none) private(i,Du,Dl,s) shared(j,q,delta_uj,delta_lj,nxs,nxe,nys,nye,phi_mus,k_mus)
   do j=nys-1,nye+1
      do i=nxs,nxe
         Du(1:dimw)=q(:,i,j+1)-q(:,i,j  )
         Dl(1:dimw)=q(:,i,j  )-q(:,i,j-1)

         !s=(2d0*Du*Dl+epsln)/(Du**2+Dl**2+epsln)
         s=(2d0*Du*Dl)/(Du**2+Dl**2+epsln)

         delta_uj(:,i,j)=0.5d0*s*((1d0+k_mus*s)*Du+(1d0-k_mus*s)*Dl)
         delta_lj(:,i,j)=0.5d0*s*((1d0+k_mus*s)*Dl+(1d0-k_mus*s)*Du)
      end do
   end do
   !$omp end parallel do

   !$omp parallel do default(none) private(i) shared(j,q,delta_uj,delta_lj,qHlj,qHrj,nxs,nxe,nys,nye,phi_mus,k_mus)
   do j=nys-1,nye
      do i=nxs,nxe
         qHlj(:,i,j)=q(1:dimw,i,j  )+phi_mus*0.5d0*delta_uj(1:dimw,i,j  )
         qHrj(:,i,j)=q(1:dimw,i,j+1)-phi_mus*0.5d0*delta_lj(1:dimw,i,j+1)
      end do
   end do
   !$omp end parallel do
end subroutine set_HV
end module half_value

subroutine read_muscl_parameter
   use half_value
   use mod_mpi
   implicit none
   character*100 buf,buf2
   integer i

   if(myid .eq. 0) then
      open(8,file="control.inp")
      do
         read(8,'(a)',end=99) buf
         if(buf(1:25) .eq. 'MUSCL ON(1)/OFF(0)      :') then
            read(buf(26:),*) i
            phi_mus = dble(i)
         else if(buf(1:25) .eq. 'MUSCL Precision Order   :') then
            read(buf(26:),*) i
            exit
         end if
      end do

99    close(8)

      select case(i)
         case(2)
            k_mus = -1d0
         case(3)
            k_mus =  1d0/3d0
         case default
            print *,"Odd Muscl Precision Order. 2 or 3. : value = ",i
            stop
      end select

      if(phi_mus .ne. 0d0 .and. phi_mus .ne. 1d0) then
         print *,"Odd Muscl Parameter.:phi_mus = ",phi_mus
         stop
      end if
   end if

   call MPI_Bcast(phi_mus,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(k_mus  ,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end subroutine read_muscl_parameter


