subroutine time_start(timenow)
   implicit none
   real,intent(out)::timenow
   real vcpu(2)
   real,external::etime

   !timenow=etime(vcpu)
   timenow = -1d0
end subroutine time_start

subroutine time_barrier(timebefore,ediff)
   use mod_mpi
   implicit none
   real,intent(in)::timebefore
   real,intent(inout)::ediff
   real vcpu(2),e2
   real,external::etime

   !e2=etime(vcpu)
   !e2=e2-timebefore
   !ediff=ediff+e2
   !call MPI_Barrier(MPI_COMM_WORLD,ierr)
end subroutine time_barrier

subroutine time_barrier_mono(ediff)
   use mod_mpi
   implicit none
   real,intent(inout)::ediff
   real vcpu(2),e1,e2
   real,external::etime

   !e1=etime(vcpu)
   !call MPI_Barrier(MPI_COMM_WORLD,ierr)
   !e2=etime(vcpu)
   !e2=e2-e1
   !ediff=ediff+e2
end subroutine time_barrier_mono
