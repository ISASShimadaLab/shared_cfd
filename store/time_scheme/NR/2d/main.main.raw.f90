part_primitive
      call set_dt(i)

      if(myid .eq. 0 .and. ILwrite) open(66,file="internal_res.dat",status="replace")

      qpp=qp
      qp =q

      j=1
      do while(j<NumInternalLoop)
part_primitive
part_main
         call set_ABpm(j)
         call set_dqdt
         call calc_next_step_implicit(j)
         j=j+1
      end do

      if(myid .eq. 0 .and. ILwrite) close(66)
      call check_internal_loop_convergence(i,step_res)
