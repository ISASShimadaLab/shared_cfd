module mod_mpi
   use grbl_prmtr
   implicit none
   integer nxs,nys,nxe,nye
   integer MPI_LOGICAL,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_COMM_WORLD,MPI_SUM,MPI_MIN,MPI_CHARACTER
   integer ierr,ireq,istatus(1)
   integer ierrs
   integer,parameter::myid=0
   integer,parameter::Nproc=1
   integer,parameter::bw =ni+nj
   integer,parameter::ngx=1
   integer,parameter::ngy=1
   integer,parameter::ng =1
   integer,parameter::bwx=ni
   integer,parameter::bwy=nj
   integer,parameter::bwxs=ni/ngx
   integer,parameter::bwys=nj/ngy
   integer,parameter::gxl =ni-bwxs*ngx
   integer,parameter::gyl =nj-bwys*ngy
   integer,parameter::gx =1
   integer,parameter::gy =1
   integer,parameter::gn =-1
   integer,parameter::gs =-1
   integer,parameter::ge =-1
   integer,parameter::gw =-1
end module mod_mpi

subroutine set_MPI
   use mod_mpi
   implicit none
   nxs=1
   nys=1
   nxe=ni
   nye=nj
end subroutine set_MPI

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
   resttime=1000000000000000
end subroutine chkelapse
