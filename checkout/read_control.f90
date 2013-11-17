subroutine read_control
   use mod_mpi
   use variable
   implicit none

   if(myid .eq. 0) then
      open(8,file="control.inp")

      read(8,'()')
      read(8,'(25x,es15.7)') max_step
      read(8,'(25x,es15.7)') min_RMS
      read(8,'(25x,es15.7)') cfl
      read(8,'(25x,i10)')    OutPeriod

      close(8)

      if(max_step <= 0 .or. min_RMS <= 0 .or. cfl <= 0 .or. OutPeriod <= 0) then
         print *,"Negative Values. max_step = ",max_step,"min_RMS = ", min_RMS, "cfl=",cfl, "OutPeriod=",OutPeriod
         stop
      end if
      if(cfl < 1d-3) then
         print *,"cfl is too small. Don't you forget '.'? cfl=",cfl
         stop
      end if
   end if


   call MPI_Bcast(max_step ,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(min_RMS  ,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(cfl      ,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(OutPeriod,1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end subroutine read_control
