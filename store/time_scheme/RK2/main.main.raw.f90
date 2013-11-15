      qp = q
part_primitive
      call set_dt(i)
part_main
      call calc_next_step_RK2_first
part_primitive
part_main
      call calc_next_step_RK2_second
