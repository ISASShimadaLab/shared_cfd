!init_package
subroutine init_therm!{{{
   call read_chemkin_parameter
   call set_reac_and_therm_data
   call set_trans_data
   call check_chem_and_flow
   call read_fo_composition
end subroutine init_therm!}}}

!data consistency check between flow and chemistry
subroutine check_chem_and_flow!{{{
   use grbl_prmtr
   use const_chem
   implicit none
   if(nY .ne. ns) then
      print *,"Error. nY and ns is different. (nY,ns)=",nY,ns
      stop
   end if
end subroutine check_chem_and_flow!}}}

! flame sheet model
subroutine set_thermo_prop!{{{
   use mod_mpi
   use prmtr
   use variable
   use chem_var
   implicit none
   double precision ei,T,MW,kappa,xi
   integer i,j,plane

   do plane = nps,npe
      !$omp parallel do private(i,ei,T,MW,kappa,xi)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            !values necessary to calculate thermo values
            ei=q(nY+3,i,j,plane)/w(1,i,j,plane)-0.5d0*(w(2,i,j,plane)**2+w(3,i,j,plane)**2)
            T=w(4,i,j,plane)/(w(1,i,j,plane)*w(indxR,i,j,plane))

            !calc therm
            call calc_T(q(1:nY,i,j,plane),w(1,i,j,plane),ei, T, kappa,MW,DHi(:,i,j,plane),vhi(:,i,j,plane),w(indxMu,i,j,plane))

            !calculate R_gas and w
            w(indxR, i,j,plane)=R_uni/MW
            w(     4,i,j,plane)=w(1,i,j,plane)*w(indxR,i,j,plane)*T
            w(indxg, i,j,plane)=kappa
            w(indxht,i,j,plane)=(q(nY+3,i,j,plane)+w(4,i,j,plane))/w(1,i,j,plane)
         end do
      end do
      !$omp end parallel do
   end do
end subroutine set_thermo_prop!}}}

! reaction
subroutine proceed_reaction!{{{
   use mod_mpi
   use variable
   implicit none
   integer i,j,plane
   double precision T

   do plane = nps,npe
      !$omp parallel do private(i,T)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            T=w(4,i,j,plane)/w(1,i,j,plane)/w(indxR,i,j,plane)
            if(T>1500d0 .or. (w(5,i,j,plane)<1d0-1d-11 .and. T>700d0)) then
               call reaction(T,q(1:nY,i,j,plane),dt_grbl)
            !   write(myid+100,'(2(i4.4,x),2(es15.7,x),i1,x,i3.3)') i,j,x(i,j),r(i,j),1,myid
            !else
            !   write(myid+100,'(2(i4.4,x),2(es15.7,x),i1,x,i3.3)') i,j,x(i,j),r(i,j),0,myid
            end if
         end do
      end do
      !$omp end parallel do
   end do

   !close(myid+100)

   !call MPI_Barrier(MPI_COMM_WORLD)
   !call exit(0)
end subroutine proceed_reaction!}}}
subroutine distribute_reaction!{{{
   use mod_mpi
   use variable
   implicit none
   integer i,j,k,kk,plane
   double precision T

   double precision tmp_send((nY+1)*bwmax*bwmax)
   integer data_send(3,bwmax*bwmax)
   integer tocalc

   integer num_send( 2,Nproc)
   integer num_receive(0:Nproc-1)
   integer list(2,1:Nproc)
   integer buf,buf2
   integer nsend,ndata
   integer bwrec,resrec
   integer res1,res2

   integer lists(2,1:Nproc)
   integer listr(2,1:Nproc)
   integer nums,numr

   double precision qrecv(nY+1,bwmax*bwmax)
   integer v_ireq(bwmax*bwmax)
   integer v_ierr(bwmax*bwmax)

   !!!! search S calculation
   nums=0
   tocalc=0
   do plane = nps,npe
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            T=w(4,i,j,plane)/w(1,i,j,plane)/w(indxR,i,j,plane)
            if(T>1500d0 .or. (w(5,i,j,plane)<1d0-1d-11 .and. T>700d0)) then
               tmp_send(tocalc*(nY+1)+1: tocalc*(nY+1)+nY)=q(1:nY,i,j,plane)
               tmp_send(tocalc*(nY+1)+nY+1)=T
               tocalc=tocalc + 1
               data_send(1,tocalc)=i
               data_send(2,tocalc)=j
               data_send(3,tocalc)=plane
            end if
         end do
      end do
   end do
   !!!! end search S calculation

   if(myid .eq. 0) then
      ! receive each data
      nsend=0
      ndata=0
      do i=0,Nproc-1
         if(i .eq. 0) then
            buf=tocalc
         else
            call MPI_Recv(buf,1,MPI_INTEGER,i,0,MPI_COMM_WORLD, istatus, ierr)
         end if

         if(buf>0) then
            nsend=nsend+1
            ndata=ndata+buf
            num_send(1,nsend)=i
            num_send(2,nsend)=buf
         end if
      end do
      ! end receive each data

      !calc num_receive
      bwrec=ndata/Nproc
      resrec=ndata-Nproc*bwrec
      do i=0,Nproc-1
         if(i<resrec) then
            num_receive(i)=bwrec+1
         else
            num_receive(i)=bwrec
         end if
      end do

      !make send list
      j=0
      res2=0
      do i=1,nsend
         res1=num_send(2,i)
         k=0
         do
            k=k+1
            if(num_receive(j)-res2>=res1) then
               list(1,k)=j
               list(2,k)=res1

               res2=res2+res1
               if(num_receive(j) .eq. res2) then
                  res2=0
                  j=j+1
               end if
               exit
            else
               list(1,k)=j
               list(2,k)=num_receive(j)-res2

               res1=res1-list(2,k)
               res2=0
               j=j+1
            end if
         end do

         if(num_send(1,i) .eq. 0) then
            nums        =k
            lists(:,1:k)=list(:,1:k)
         else
            call MPI_Send(   k,  1,MPI_INTEGER,num_send(1,i),0,MPI_COMM_WORLD,ierr)
            call MPI_Send(list,2*k,MPI_INTEGER,num_send(1,i),0,MPI_COMM_WORLD,ierr)
         end if
      end do
      !end make send list

      !make receive list
      j=1
      res2=0
      do i=0,Nproc-1
         res1=num_receive(i)
         k=0
         do
            k=k+1
            if(num_send(2,j)-res2>=res1) then
               list(1,k)=num_send(1,j)
               list(2,k)=res1

               res2=res2+res1
               if(num_send(2,j) .eq. res2) then
                  res2=0
                  j=j+1
               end if
               exit
            else
               list(1,k)=num_send(1,j)
               list(2,k)=num_send(2,j)-res2

               res1=res1-list(2,k)
               res2=0
               j=j+1
            end if
         end do

         if(i .eq. 0) then
            numr        =k
            listr(:,1:k)=list(:,1:k)
         else
            call MPI_Send(k,     1,MPI_INTEGER,i,0,MPI_COMM_WORLD,ierr)
            call MPI_Send(list,2*k,MPI_INTEGER,i,0,MPI_COMM_WORLD,ierr)
         end if
      end do
      !end make receive list
   else
      call MPI_Send(tocalc,  1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
      if(tocalc .ne. 0) then
         call MPI_Recv( nums,     1,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
         call MPI_Recv(lists,2*nums,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
      end if
      call MPI_Recv( numr,     1,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
      call MPI_Recv(listr,2*numr,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
   end if

   !!!!!!!!!!!!!!!!!!!! SEND, CALC, RECV !!!!!!!!!!!!!!!!!
   buf=1
   do i=1,nums
      if(lists(1,i) .ne. myid) then
         call MPI_Isend(tmp_send(buf),lists(2,i)*(nY+1), MPI_DOUBLE_PRECISION,&
                         lists(1,i),0,MPI_COMM_WORLD,v_ireq(i),v_ierr(i))
         buf=buf+lists(2,i)*(nY+1)
      else
         buf2=0
         do j=1,numr
            if(listr(1,j) .ne. myid) then
               buf2=buf2+listr(2,j)
            else
               do k=1,listr(2,j)
                  qrecv(:,buf2+k)=tmp_send(buf:buf+nY)
                  buf=buf+nY+1
               end do
               exit
            end if
         end do
      end if
   end do

   buf=1
   do i=1,numr
      if(listr(1,i) .ne. myid) then
         call MPI_Recv(qrecv(1,buf),listr(2,i)*(nY+1), MPI_DOUBLE_PRECISION,&
                         listr(1,i),0,MPI_COMM_WORLD, istatus, ierr)
      end if
      buf=buf+listr(2,i)
   end do

   !$omp parallel do private(i)
   do i=1,buf-1
      call reaction(qrecv(nY+1,i),qrecv(1:nY,i),dt_grbl)
   end do
   !$omp end parallel do

   do i=1,nums
      if(lists(1,i) .ne. myid) call MPI_Wait(v_ireq(i),istatus,ierr)
   end do

   buf=1
   do i=1,numr
      if(listr(1,i) .ne. myid) then
         call MPI_Isend(qrecv(1,buf),listr(2,i)*(nY+1),MPI_DOUBLE_PRECISION,&
                         listr(1,i),0,MPI_COMM_WORLD,v_ireq(i),v_ierr(i))
      else
         buf2=0
         do j=1,nums
            if(lists(1,j) .ne. myid) then
               buf2=buf2+lists(2,j)
            else
               do k=1,lists(2,j)
                  q(1:nY,data_send(1,buf2+k),&
                         data_send(2,buf2+k),&
                         data_send(3,buf2+k))&
                         =qrecv(1:nY,buf-1+k)
               end do
               exit
            end if
         end do
      end if
      buf=buf+listr(2,i)
   end do

   buf=0
   do i=1,nums
      if(lists(1,i) .ne. myid) then
         call MPI_Recv(tmp_send,lists(2,i)*(nY+1), MPI_DOUBLE_PRECISION,&
                         lists(1,i),0,MPI_COMM_WORLD, istatus,ierr)
         do j=1,lists(2,i)
            q(1:nY,data_send(1,buf+j),data_send(2,buf+j),data_send(3,buf+j))&
                    =tmp_send((j-1)*(nY+1)+1:(j-1)*(nY+1)+nY)
         end do
      end if
      buf=buf+lists(2,i)
   end do

   do i=1,numr
      if(listr(1,i) .ne. myid) call MPI_Wait(v_ireq(i),istatus,ierr)
   end do

   !print '(a,i3.2,a)',"ID:",myid," OK"
   !call MPI_Barrier(MPI_COMM_WORLD,ierr)
   !call exit(0)
end subroutine distribute_reaction!}}}

!boundary condition
subroutine calc_boundary(p,T,deg, wt,vhi)!{{{
   use grbl_prmtr
   use prmtr
   use const_chem
   implicit none
   double precision,intent(in)::p
   double precision,intent(inout)::T
   double precision,intent(in)::deg

   double precision,intent(out)::wt(dimw)
   double precision,intent(out)::vhi(ns)

   double precision,dimension(ns)::vrho,DHi
   double precision rho,E,MWave
   integer j

   call calc_vrho(p,T,deg, rho,vrho,E)
   call calc_T(vrho,rho,E, T, wt(indxg),MWave,DHi,vhi,wt(indxMu))

   !set wt
   wt(1)=rho
   wt(2)=0d0
   wt(3)=0d0
   wt(4)=p
   do j=1,ns
      wt(4+j)=vrho(j)/rho
   end do
   wt(indxht)=E+p/rho
   wt(indxR) =R_uni/MWave
end subroutine calc_boundary!}}}
