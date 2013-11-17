logical function check_convergence(step,step_res)
   use mod_mpi
   use grbl_prmtr
   use variable
   implicit none
   integer,intent(in)::step
   integer,intent(in)::step_res
   double precision tmp,tmpt,sm
   integer i,j,plane
   integer*8 rest_time

   if(mod(step,OutPeriod) .ne. 0) then
      check_convergence = .true.
      return
   end if

   tmp =0d0
   sm = 0d0
   do plane = nps,npe
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            tmp =tmp +(q(nY+3,i,j,plane)-qp(nY+3,i,j,plane))**2
         end do
      end do
      sm = sm + dble((nxe(plane)-nxs(plane)+1)*(nye(plane)-nys(plane)+1))
   end do

   tmpt=sqrt(tmp)/sm
   if(tmpt .ne. tmpt) then
      print *,"NaN at procid=",myid
      call exit(1)
   end if
   call MPI_Reduce(tmpt, tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

   if(myid .eq. 0) then
      !RMS check file
      open(82,file="RMS.dat",position='append')
      write(82,*) step,tmp
      close(82)
      print '(a,i12,x,2(a,es15.7,x))',"step=",step,"log(RMS)=",log10(tmp),"dt_grbl=",dt_grbl

      call chkelapse(rest_time)

      check_convergence = .true.
      if(min_RMS>=tmp .and. step > 5000) then
         print *,"RMS accomplished at step=",step," at tmp=",tmp
         check_convergence = .false.
      else if(max_step+step_res<=step) then
         print *,"Error:Exceed Max Step. NOT CONVERTED."
         check_convergence = .false.
      else if(rest_time < 120) then
         print *,"MPI system time over."
         check_convergence = .false.
      end if
   end if

   call MPI_Bcast(check_convergence, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   call out_bin(step)
end function check_convergence!}}}

