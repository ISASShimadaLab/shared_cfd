      i=i+1
   end do

   if(myid .eq. 0) print *,"normal end"
   CALL MPI_Finalize(ierr)
stop
end program main
