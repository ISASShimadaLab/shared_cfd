module mod_mpi
   use grbl_prmtr
   implicit none
   integer nps,npe
   integer,dimension(Nplane)::nxs,nys,nxe,nye
   integer MPI_LOGICAL,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_COMM_WORLD,MPI_SUM,MPI_MIN,MPI_CHARACTER,MPI_MAX
   integer,dimension(Nplane)::ireq,istatus,ierrs
   integer ierr
   integer,parameter::myid=0
   integer,parameter::Nproc=1
   integer,parameter::bwmax =1 ! THIS VALUE IS NOT USED IN PERSONAL COMPUTER.
   integer,parameter::ngxmax=1
   integer,parameter::ngymax=1
   integer,parameter::ngx(Nplane)=  1
   integer,parameter::ngy(Nplane)=  1
   integer,parameter::gx(Nplane) =  1
   integer,parameter::gy(Nplane) =  1
   integer,parameter::gn(Nplane) = -1
   integer,parameter::gs(Nplane) = -1
   integer,parameter::ge(Nplane) = -1
   integer,parameter::gw(Nplane) = -1
   integer bwx(Nplane)
   integer bwy(Nplane)
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
   integer plane

   do plane=1,Nplane
      bwx(plane) = ni(plane)
      bwy(plane) = nj(plane)
      nxs(plane) = 1
      nys(plane) = 1
      nxe(plane) = ni(plane)
      nye(plane) = nj(plane)
   end do

   nps=1
   npe=Nplane

   proc_list(:,:,:)=0

   call set_cut_copro
end subroutine set_MPI

subroutine set_cut_copro
   use mod_mpi
   implicit none
   character*100 line
   integer*4,external::access
   integer i,j

   if(access('cut_copro000.inp',"r") .ne. 0) then
        num_cut_copro=0
        return
   end if
   open(45,file="cut_copro000.inp")
   read(45,'(a)') !for comment line
   read(45,*) num_cut_copro
   allocate(cut_copro(11,num_cut_copro))
   read(45,'(a)') !for comment line
   do i=1,num_cut_copro
      read(45,'(11i6)') cut_copro(:,i)
   end do
   close(45)
end subroutine set_cut_copro

subroutine     MPI_Finalize(ierr)
   integer ierr
end subroutine MPI_Finalize

subroutine     MPI_Allreduce(f, t, nod, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   double precision f, t
   integer  nod, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr
   t=f
end subroutine MPI_Allreduce

subroutine     MPI_Bcast(var, siz, MPI_DOUBLE_PRECISION, nod, MPI_COMM_WORLD, ierr)
   integer var
   integer siz, MPI_DOUBLE_PRECISION, nod, MPI_COMM_WORLD, ierr
end subroutine MPI_Bcast

subroutine     MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gw,tag,MPI_COMM_WORLD, istatus, ierr)
   double precision tmp
   integer cnt,MPI_DOUBLE_PRECISION,gw,tag,MPI_COMM_WORLD, istatus, ierr
end subroutine MPI_Recv
subroutine     MPI_Send(tmp,cnt,MPI_DOUBLE_PRECISION,ge,tag,MPI_COMM_WORLD,ierr)
   double precision tmp
   integer cnt,MPI_DOUBLE_PRECISION,ge,tag,MPI_COMM_WORLD,ierr
end subroutine MPI_Send

subroutine MPI_Isend(tmps,cnts,MPI_DOUBLE_PRECISION,ge,ge2,MPI_COMM_WORLD,ireq,ierrs)
    double precision tmps
    integer cnts,MPI_DOUBLE_PRECISION,ge,ge2,MPI_COMM_WORLD,ireq,ierrs
end subroutine MPI_Isend

subroutine MPI_Wait(ireq,istatus,ierrs)
    double precision ireq
    integer istatus,ierrs
end subroutine MPI_Wait

subroutine MPI_Reduce(tmpt, tmp, te, MPI_DOUBLE_PRECISION, MPI_SUM, te2, MPI_COMM_WORLD, ierr)
    double precision tmpt, tmp
    integer te1, MPI_DOUBLE_PRECISION, MPI_SUM, te2, MPI_COMM_WORLD, ierr
    tmp=tmpt
end subroutine MPI_Reduce

subroutine     chkelapse(resttime)
   integer*8 resttime
   resttime=1000000000
end subroutine chkelapse
