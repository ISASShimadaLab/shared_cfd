subroutine set_ABpm
   use grbl_prmtr
   use mod_mpi
   use variable
   use var_lusgs

   implicit none
   double precision rho,u,v,H,c,p,Vs
   double precision Ur,Vv,alp,ud,cd,phh
   double precision tmp,uds,dtau
   double precision unA, unB, nuA, nuB
   double precision,dimension(dimq,dimq)::A,B,Atilde,Btilde
   double precision,dimension(dimq)::tA,tB
   double precision gmmad,gamm
   integer i,j,k,l

   !$omp parallel do private(i,k,l,&
   !$omp                     rho,u,v,H,Vs,&
   !$omp                     Ur,Vv,alp,ud,cd,phh,&
   !$omp                     p,c,tmp,&
   !$omp                     unA,unB,nuA,nuB,&
   !$omp                     A,B,Atilde,Btilde,&
   !$omp                     gmmad,gamm,uds,dtau,&
   !$omp                     tA,tB)
   do j=nys,nye
      do i=nxs,nxe
         !set w
         rho  = w(     1,i,j)
         u    = w(     2,i,j)
         v    = w(     3,i,j)
         p    = w(     4,i,j)
         gamm = w(indxg ,i,j)
         H    = w(indxht,i,j)
         c    = sqrt(gamm*p/rho)
         gmmad= gamm-1d0

         Vs   = u**2+v**2

         !unA,unB,nuA,nuB
         unA = vnci(1,i,j)*u+vnci(2,i,j)*v
         unB = vncj(1,i,j)*u+vncj(2,i,j)*v

         Vv  = w(indxMu,i,j)/rho/min(dsci(i,j),dscj(i,j))
         Ur  = min(max(Vv,sqrt(Vs)),c)
         phi(i,j) = 1d0/Ur**2 -1d0/c**2
         alp = 0.5d0*(1d0-(Ur/c)**2)

         ud  = unA*(1d0-alp)
         cd  = sqrt(alp**2 *unA**2 +Ur**2)
         nuA = abs(ud)+cd

         ud  = unB*(1d0-alp)
         cd  = sqrt(alp**2 *unB**2 +Ur**2)
         nuB = abs(ud)+cd

         !set new dt
         uds= nuA*(dsi(i,j)+dsi(i-1,j))+nuB*(dsj(i,j)+dsj(i,j-1))
         dtau = cfl_tau*Vol(i,j)/max(uds,w(indxMu,i,j)/rho)

         !set preconditioner
         phh = -1d0/(c**2+1d0/phi(i,j))

         do k=1,nY
            pre1(k,i,j) = w(4+k,i,j)
         end do
         pre1(nY+1,i,j) = u
         pre1(nY+2,i,j) = v
         pre1(nY+3,i,j) = H

         do k=1,nY
            pre2(k,i,j) = 0.5d0*gmmad*Vs+DHi(k,i,j)
         end do
         pre2(nY+1,i,j) =-gmmad*u
         pre2(nY+2,i,j) =-gmmad*v
         pre2(nY+3,i,j) = gmmad


         !set A
         A(nY+1,nY+1)= u*(3d0-gamm)
         A(nY+2,nY+1)= v
         A(nY+3,nY+1)= H-u**2*gmmad

         A(nY+1,nY+2)=-v*gmmad
         A(nY+2,nY+2)= u
         A(nY+3,nY+2)=-u*v*gmmad

         A(nY+1,nY+3)= gmmad
         A(nY+2,nY+3)= 0d0
         A(nY+3,nY+3)= u*gamm

         !set tA
         do k=1,nY
            tA(k)     =-w(4+k,i,j)*u
            A( k,nY+1)= w(4+k,i,j)
            A( k,nY+2)= 0d0
            A( k,nY+3)= 0d0
         end do
         tA(nY+1)=-u**2 +0.5d0*gmmad*Vs
         tA(nY+2)=-u*v
         tA(nY+3)= u*(-H+0.5d0*gmmad*Vs)

         !set B
         B(nY+1,nY+1)= 0d0
         B(nY+1,nY+2)= 1d0
         B(nY+1,nY+3)= 0d0

         B(nY+1,nY+1)= v
         B(nY+2,nY+1)=-u*gmmad
         B(nY+3,nY+1)=-u*v*gmmad

         B(nY+1,nY+2)= u
         B(nY+2,nY+2)= v*(3d0-gamm)
         B(nY+3,nY+2)= H-v**2*gmmad

         B(nY+1,nY+3)= 0d0
         B(nY+2,nY+3)= gmmad
         B(nY+3,nY+3)= v*gamm

         !set tB
         do k=1,nY
            tB(k)     =-w(4+k,i,j)*v
            B( k,nY+1)= 0d0
            B( k,nY+2)= w(4+k,i,j)
            B( k,nY+3)= 0d0
         end do
         tB(nY+1)=-v*u
         tB(nY+2)=-v**2 +0.5d0*gmmad*Vs
         tB(nY+3)= v*(-H+0.5d0*gmmad*Vs)

         do k=1,nY
            do l=1,nY
               A(l,k)=tA(l)
               B(l,k)=tB(l)
            end do

            A(k   ,k)=tA(k   )+u
            A(nY+1,k)=tA(nY+1)+DHi(k,i,j)
            A(nY+2,k)=tA(nY+2)
            A(nY+3,k)=tA(nY+3)+DHi(k,i,j)*u

            B(k   ,k)=tB(k   )+v
            B(nY+1,k)=tB(nY+1)
            B(nY+2,k)=tB(nY+2)+DHi(k,i,j)
            B(nY+3,k)=tB(nY+3)+DHi(k,i,j)*v
         end do


         !multiply preM
         do k=1,dimq
            tmp = 0d0
            do l=1,dimq
               tmp = tmp + pre2(l,i,j)*A(l,k)
            end do

            tmp = phh * tmp

            do l=1,dimq
               A(l,k) = A(l,k) + pre1(l,i,j)*tmp
            end do
         end do

         do k=1,dimq
            tmp = 0d0
            do l=1,dimq
               tmp = tmp + pre2(l,i,j)*B(l,k)
            end do

            tmp = phh * tmp

            do l=1,dimq
               B(l,k) = B(l,k) + pre1(l,i,j)*tmp
            end do
         end do


         !set Atilde, Apm
         Atilde = vnci(1,i,j)*A+vnci(2,i,j)*B

         Ap(:,:,i,j)= Atilde*0.5d0
         Am(:,:,i,j)= Atilde*0.5d0

         do k=1,dimq
            Ap(k,k,i,j)=Ap(k,k,i,j)+nuA*0.5d0
            Am(k,k,i,j)=Am(k,k,i,j)-nuA*0.5d0
         end do

         !set Btilde, Bpm
         Btilde= vncj(1,i,j)*A+vncj(2,i,j)*B

         Bp(:,:,i,j)= Btilde*0.5d0
         Bm(:,:,i,j)= Btilde*0.5d0

         do k=1,dimq
            Bp(k,k,i,j)=Bp(k,k,i,j)+nuB*0.5d0
            Bm(k,k,i,j)=Bm(k,k,i,j)-nuB*0.5d0
         end do

         !set alpha and beta
         tmp = Vol(i,j)/dtau + dsci(i,j)*nuA + dscj(i,j)*nuB
         alpha(i,j)= 1d0/tmp
         beta( i,j)= 3d0*Vol(i,j)/(2d0*DT_LOCAL_GLOBAL) *alpha(i,j)

         !set new pre1
         phh = -1d0/(c**2+(1d0+beta(i,j))/phi(i,j))
         pre1(:,i,j) = phh * pre1(:,i,j)
      end do
   end do
   !$omp end parallel do
end subroutine set_ABpm

subroutine set_dqdt
   use grbl_prmtr
   use mod_mpi
   use variable
   use var_lusgs
   implicit none
   integer i,j

   if(first_dual) then
      !$omp parallel do private(i)
      do j=nys,nye
         do i=nxs,nxe
            dqdt(:,i,j)=(q(:,i,j)-qp(:,i,j))/DT_LOCAL_GLOBAL
         end do
      end do
      !$omp end parallel do
   else
      !$omp parallel do private(i)
      do j=nys,nye
         do i=nxs,nxe
            dqdt(:,i,j)=(3d0*q(:,i,j)-4d0*qp(:,i,j)+qpp(:,i,j))/(2d0*DT_LOCAL_GLOBAL)
         end do
      end do
      !$omp end parallel do
   end if
end subroutine set_dqdt

subroutine set_vnc_dsc
   use grbl_prmtr
   use mod_mpi
   use variable
   use var_lusgs
   implicit none
   double precision,dimension(2)::a,b,c,d
   integer i,j

   do j=nys,nye
      do i=nxs,nxe
         a(1)=(xh(i-1,j-1)+xh(i,j-1))*0.5d0
         a(2)=(rh(i-1,j-1)+rh(i,j-1))*0.5d0
         b(1)=(xh(i,j-1)  +xh(i,j)  )*0.5d0
         b(2)=(rh(i,j-1)  +rh(i,j)  )*0.5d0
         c(1)=(xh(i-1,j)  +xh(i,j)  )*0.5d0
         c(2)=(rh(i-1,j)  +rh(i,j)  )*0.5d0
         d(1)=(xh(i-1,j-1)+xh(i-1,j))*0.5d0
         d(2)=(rh(i-1,j-1)+rh(i-1,j))*0.5d0
         dsci(i,j)=sqrt((c(1)-a(1))**2+(c(2)-a(2))**2)
         dscj(i,j)=sqrt((b(1)-d(1))**2+(b(2)-d(2))**2)
         vnci(1,i,j)=(b(1)-d(1))/dscj(i,j)
         vnci(2,i,j)=(b(2)-d(2))/dscj(i,j)
         vncj(1,i,j)=(c(1)-a(1))/dsci(i,j)
         vncj(2,i,j)=(c(2)-a(2))/dsci(i,j)

         !!only when cylindrical
         !dsci(i,j)=dsci(i,j)*(c(2)+a(2))*0.5d0
         !dscj(i,j)=dscj(i,j)*(b(2)+d(2))*0.5d0
      end do
   end do
end subroutine set_vnc_dsc



subroutine check_internal_loop_convergence(step,step_res)
   use mod_mpi
   use var_lusgs
   implicit none
   integer,intent(in)::step
   integer,intent(in)::step_res

   if(myid .eq. 0 .and. res/res1 > res_rate_warning) then
      print '(a,i12,es11.3)',"WARNING:NOT CONVERGED AT INTERNAL LOOP&
                            & (step,res_rate) = ", step,res/res1
      if(res/res1 > res_rate_error .and. step .ne. step_res + 1) &
         stop "ERROR. EXCEED MAXIMUM RESIDUAL."
   end if

   first_dual = .false.
end subroutine check_internal_loop_convergence

subroutine read_InternalLoop
   use mod_mpi
   use var_lusgs
   implicit none
   character*100 buf
   integer i

   if(myid .eq. 0) then
      open(8,file="control.inp")
      do
         read(8,'(a)',end=99) buf
         if(buf(1:25)      .eq. 'InternalLoop CFL tau    :') then
            read(buf(26:),*) cfl_tau
         else if(buf(1:25) .eq. 'InternalLoop Omega Max  :') then
            read(buf(26:),*) omega_max
         else if(buf(1:25) .eq. 'InternalLoop Omega Min  :') then
            read(buf(26:),*) omega_min
         else if(buf(1:25) .eq. 'InternalLoop Dqrate Max :') then
            read(buf(26:),*) Dqmax
         else if(buf(1:25) .eq. 'InternalLoop ResRateWarn:') then
            read(buf(26:),*) res_rate_warning
         else if(buf(1:25) .eq. 'InternalLoop ResRateErr :') then
            read(buf(26:),*) res_rate_error
         else if(buf(1:25) .eq. 'InternalLoop OutOmega   :') then
            read(buf(26:),*) ILwrite
         else if(buf(1:25) .eq. 'InternalLoop Max Number :') then
            read(buf(26:),*) NumInternalLoop
         end if
      end do
99    close(8)
   end if

   call MPI_Bcast(cfl_tau,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(omega_max,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(omega_min,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(Dqmax,           1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(NumInternalLoop, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)

   first_dual = .true.
end subroutine read_InternalLoop


