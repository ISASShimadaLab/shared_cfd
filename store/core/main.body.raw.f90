
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
   i=step_res+1
part_primitive
   call set_BC(i)
   ! main loop
   do while(check_convergence(i,step_res))
