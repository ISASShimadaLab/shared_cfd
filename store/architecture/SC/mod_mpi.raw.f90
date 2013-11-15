module mod_mpi
   use grbl_prmtr
   implicit none
   include 'mpif.h'
   integer nps,npe
   integer myid,Nproc
   integer,dimension(Nplane)::ierrs,ireq,istatus(MPI_STATUS_SIZE)
   integer ierr
   integer,parameter::bwmax=BWMAX
   integer,parameter::ngxmax=NGXMAX
   integer,parameter::ngymax=NGYMAX
   integer,parameter::ngx(Nplane)=(/NGX/)
   integer,parameter::ngy(Nplane)=(/NGY/)
   integer,dimension(Nplane)::bwx,bwy
   integer,dimension(Nplane)::gx,gy
   integer,dimension(Nplane)::nxs,nys
   integer,dimension(Nplane)::nxe,nye
   integer,dimension(Nplane)::gn,gs,ge,gw
   integer nxs_mat(ngxmax,Nplane)
   integer nxe_mat(ngxmax,Nplane)
   integer nys_mat(ngymax,Nplane)
   integer nye_mat(ngymax,Nplane)
   integer proc_list(ngxmax,ngymax,Nplane)

   integer,allocatable::cut_copro(:,:)
   integer num_cut_copro

   integer,allocatable::MPIComm(:,:)
   integer NumMPIComm
end module mod_mpi

subroutine set_MPI
   use mod_mpi
   implicit none
   integer buf(7)
   integer i,j,myidt
   integer ng,plane

   CALL MPI_Init(ierr)
   CALL MPI_Comm_size(MPI_COMM_WORLD, Nproc, ierr)
   CALL MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)

   if(myid .eq. 0) then
      open(41,file="grid_separation.inp")
      do plane = 1,Nplane
         read(41,*) nxs_mat(1:ngx(plane),plane)
         read(41,*) nys_mat(1:ngy(plane),plane)
      end do
      close(41)

      ng=0
      do plane = 1,Nplane
         ng = ng+ngx(plane)*ngy(plane)
      end do

      if(Nproc .ne. ng) then
         print *,"Odd number of processors. Nproc:",Nproc," ng:",ng
         stop
      else
         print *,"Nproc=ng=",Nproc
      end if

      plane = 1; i = 1; j = 1
      do myidt=0,Nproc-1
         proc_list(i,j,plane) = myidt
         buf(1) = plane
         buf(2) = i
         buf(3) = j

         buf(4)=nxs_mat(i,plane)
         if(i .eq. ngx(plane)) then
            buf(5)=ni(plane)
         else
            buf(5)=nxs_mat(i+1,plane)-1
         end if
         nxe_mat(i,plane)=buf(5)

         buf(6)=nys_mat(j,plane)
         if(j .eq. ngy(plane)) then
            buf(7)=nj(plane)
         else
            buf(7)=nys_mat(j+1,plane)-1
         end if
         nye_mat(j,plane)=buf(7)

         if(myidt .eq. 0) then
            nps      = buf(1)
            gx(nps)  = buf(2)
            gy(nps)  = buf(3)
            nxs(nps) = buf(4)
            nxe(nps) = buf(5)
            nys(nps) = buf(6)
            nye(nps) = buf(7)
         else
            call MPI_Send(buf,7,MPI_INTEGER,myidt,0,MPI_COMM_WORLD,ierr)
         end if

         !calc next ngx,ngy,plane
         j=j+1
         if(j>ngy(plane)) then
            j=1
            i=i+1
            if(i>ngx(plane)) then
               i=1
               plane = plane +1
            end if
         end if
      end do
   else
      call MPI_Recv(buf,7,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
      nps      = buf(1)
      gx(nps)  = buf(2)
      gy(nps)  = buf(3)
      nxs(nps) = buf(4)
      nxe(nps) = buf(5)
      nys(nps) = buf(6)
      nye(nps) = buf(7)
   end if

   npe = nps
   bwx(nps)=nxe(nps)-nxs(nps)+1
   bwy(nps)=nye(nps)-nys(nps)+1

   gn(nps)=myid+1
   gs(nps)=myid-1
   ge(nps)=myid+ngy(nps)
   gw(nps)=myid-ngy(nps)

   if(bwx(nps) >bwmax .or. bwy(nps) > bwmax) then
      print *,"bwx or bwy exceed bwmax. (bwx, bwy, bwmax)=",bwx(nps),bwy(nps),bwmax
      stop
   end if

   call set_cut_copro

   call set_cut_mpi
end subroutine set_MPI

subroutine set_cut_copro
   use mod_mpi
   implicit none
   character*100 line
   integer*4,external::access
   integer i

   write(line,'("cut_copro",i3.3,".inp")') myid
   if(access(trim(line),"r") .ne. 0) then
        num_cut_copro=0
        return
   end if
   open(45,file=trim(line))
   read(45,'(a)') !for comment line
   read(45,*) num_cut_copro
   allocate(cut_copro(11,num_cut_copro))
   read(45,'(a)') !for comment line
   do i=1,num_cut_copro
      read(45,'(11i6)') cut_copro(:,i)
   end do
   close(45)
end subroutine set_cut_copro

subroutine set_cut_mpi
   use mod_mpi
   implicit none
   character*100 line
   integer*4,external::access
   integer i

   write(line,'("MPIcomm",i3.3,".inp")') myid
   if(access(trim(line),"r") .ne. 0) then
        NumMPIComm=0
        return
   end if
   open(45,file=trim(line))
   read(45,'(a)') !for comment line
   read(45,*) NumMPIComm
   allocate(MPIComm(8,NumMPIComm))
   read(45,'(a)') !for comment line
   do i=1,NumMPIComm
      read(45,'(8i6)') MPIComm(:,i)
   end do
   close(45)
end subroutine set_cut_mpi
