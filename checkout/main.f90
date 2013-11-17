program main
   use mod_mpi
   use geometry
   use variable
   use half_value
   use var_lusgs
   implicit none
   integer step_res
   integer i
   logical,external::check_convergence
   integer*4,external::access
   integer j
   call set_MPI
   call init_geometry
   call read_control          !for general parameters
   call set_vnc_dsc           !for Jacobian
   call read_InternalLoop     !for Internal Loop
   call read_muscl_parameter  !for muscl_parameters
   call set_geojac            !for viscosity terms


   ! set initial condition
   if(access("restart.bin","r") .eq. 0) then
      if(myid .eq. 0) print *,"restart mode"
      call restart_bin(step_res)
      if(myid .eq. 0) print *,"the number of step spent before is ",step_res
   else
      if(myid .eq. 0) print *,"initial mode"
      call set_IC
      step_res=0
      tt = 0d0
      if(myid .eq. 0) print *,"read IC"
   end if

      call set_w
      call set_thermo_prop
   call set_BC(i)
   ! main loop
   i=step_res+1
   do while(check_convergence(i,step_res))

      call set_w
      call set_thermo_prop

      if(myid .eq. 0 .and. ILwrite) open(66,file="internal_res.dat",status="replace")

      qpp=qp
      qp =q

      j=1
      do while(j<NumInternalLoop)

      call set_w
      call set_thermo_prop

      call set_BC(i)
      call set_HV(w,wHli,wHri,wHlj,wHrj)
      call set_TG
      call set_TGv

         call set_ABpm(j)
         call set_dqdt
         call calc_next_step_implicit(j)
         j=j+1
      end do

      if(myid .eq. 0 .and. ILwrite) close(66)
      call check_internal_loop_convergence(i,step_res)
      i=i+1
   end do

   if(myid .eq. 0) print *,"normal end"
   CALL MPI_Finalize(ierr)
end program main
