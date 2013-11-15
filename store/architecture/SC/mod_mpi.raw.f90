module mod_mpi
   use grbl_prmtr
   implicit none
   include 'mpif.h'
   integer myid,Nproc
   integer ierr,ireq,istatus(MPI_STATUS_SIZE)
   integer ierrs
   integer,parameter::bwmax=BWMAX
   integer,parameter::ngx=NGX
   integer,parameter::ngy=NGY
   integer,parameter::ng  =ngx*ngy
   integer bwx,bwy
   integer  gx, gy
   integer nxs,nys
   integer nxe,nye
   integer gn,gs,ge,gw
   integer nxs_mat(ngx)
   integer nxe_mat(ngx)
   integer nys_mat(ngy)
   integer nye_mat(ngy)
end module mod_mpi

subroutine set_MPI
   use mod_mpi
   implicit none
   integer buf(4)
   integer i,j,myidt

   CALL MPI_Init(ierr)
   CALL MPI_Comm_size(MPI_COMM_WORLD, Nproc, ierr)
   CALL MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)

   gx=myid/ngy+1
   gy=myid+1-(gx-1)*ngy

   gn=myid+1
   gs=myid-1
   ge=myid+ngy
   gw=myid-ngy

   if(myid .eq. 0) then
      open(41,file="grid_separation.inp")
      read(41,*) nxs_mat(:)
      read(41,*) nys_mat(:)
      close(41)

      do j=1,ngy
         do i=1,ngx
            myidt=j+ngy*(i-1)-1

            buf(1)=nxs_mat(i)
            if(i .eq. ngx) then
               buf(2)=ni
            else
               buf(2)=nxs_mat(i+1)-1
            end if
            nxe_mat(i)=buf(2)

            buf(3)=nys_mat(j)
            if(j .eq. ngy) then
               buf(4)=nj
            else
               buf(4)=nys_mat(j+1)-1
            end if
            nye_mat(j)=buf(4)

            if(myidt .eq. 0) then
               nxs=buf(1)
               nxe=buf(2)
               nys=buf(3)
               nye=buf(4)
               bwx=nxe-nxs+1
               bwy=nye-nys+1
            else
               call MPI_Send(buf,4,MPI_INTEGER,myidt,0,MPI_COMM_WORLD,ierr)
            end if
         end do
      end do
   else
      call MPI_Recv(buf,4,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
      nxs=buf(1)
      nxe=buf(2)
      nys=buf(3)
      nye=buf(4)
      bwx=nxe-nxs+1
      bwy=nye-nys+1
   end if

   if(bwx >bwmax .or. bwy > bwmax) then
      print *,"bwx or bwy exceed bwmax. (bwx, bwy, bwmax)=",bwx,bwy,bwmax
   end if
   if(Nproc .ne. ng) then
      if(myid .eq. 0) print *,"Odd number of processors. Nproc:",Nproc," ng:",ng
      stop
   else
      if(myid .eq. 0) print *,"Nproc=ng=",Nproc
   end if
end subroutine set_MPI
