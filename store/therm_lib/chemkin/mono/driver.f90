subroutine plot_time_history(T,vrho,dt,tout,outspc,Noutspc)
   use const_chem
   implicit none
   double precision,intent(inout)::T
   double precision,intent(inout)::vrho(ns)
   double precision,intent(in)   ::dt
   double precision,intent(in)   ::tout
   integer         ,intent(in)   ::outspc(max_ns)
   integer         ,intent(in)   ::Noutspc

   double precision tt,to
   double precision n(ns+1)
   double precision RWORK(LRW)
   integer          IWORK(LIW)
   integer          istate,i,j
   double precision atol
   logical          flag

   double precision rho

   tt    =0d0
   istate=1

   open(22,file="plt.dat")
   write(22,'(2es15.7)') tt,T
   i=1
   do
      to=dble(i)*dt
      if(to>=tout) to=tout
      call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol,flag)
      if(.not.flag) stop "Calculation failed at plot_time_history."
      write(22,'(2es15.7)',advance='no') tt,T
      if(Noutspc>0) then
         rho=0d0
         do j=1,ns
            rho=rho+vrho(j)
         end do
         rho=1d0/rho
         do j=1,Noutspc
            write(22,'(es15.7)',advance='no') vrho(outspc(j))*rho
         end do
      end if
      write(22,'()')

      if(to>=tout) exit
      i=i+1
   end do
   close(22)
end subroutine plot_time_history

subroutine read_controlinp(T,vrho,dt,tout,outspc,Noutspc)
   use chem
   implicit none
   double precision,intent(in)::T,vrho(max_ns)
   double precision,intent(out)::dt,tout
   integer         ,intent(out)::outspc(max_ns)
   integer         ,intent(out)::Noutspc

   double precision Tref(2)
   character*100 buf
   logical flag
   integer i

   !read tout,dt
   open(22,file="control.inp")
   read(22,'(a100)') buf
   read(buf(26:),*) tout
   read(22,'(a100)') buf
   read(buf(26:),*) dt

   if(tout <0d0) then
      print '(a)',"Time range will be determined from ignition time."
      call reaction_one_steps(T,vrho,0,Tref,flag)
      if(.not. flag) stop "Calculation failed at ignition time determination."
      tout=Tref(1)*1.5d0
      dt  =tout*1d-3
      print '(a,es15.7)',"Ignition Time(s)=",Tref(1)
   end if
   if(dt<0d0) stop "Negative time step."

   !read species
   read(22,'(a100)') buf
   Noutspc=1
   do
      read(22,'(a100)') buf
      if(buf .eq. 'all') then
         do i=1,ns
            outspc(i)=i
         end do
         Noutspc=ns+1
         exit
      end if
      if(buf .eq. 'end') exit
      call search_SYM(buf,outspc(Noutspc))
      if(outspc(Noutspc)<0) exit
      Noutspc=Noutspc+1
   end do
   Noutspc=Noutspc-1
   print '(a,i4)',"The number of species to plot=",Noutspc

   close(22)
contains
subroutine search_SYM(str,i)!{{{ 
   implicit none
   character(*),intent(inout)          ::str
   integer,intent(out)                 ::i
   do i=1,ns
      if(SYM_SPC(i) .eq. buf) return
   end do
   i=-1
end subroutine search_SYM!}}}
end subroutine read_controlinp

subroutine out_plotplt(outspc,Noutspc)
   use chem
   implicit none
   integer         ,intent(in)::outspc(max_ns)
   integer         ,intent(in)::Noutspc

   integer i

   open(22,file="plot.plt")
   write(22,'(a)') "#!/usr/bin/gnuplot"
   write(22,'(a)') "plot 'plt.dat' u 1:2 w l title 'temperature'"
   write(22,'(a)') "pause -1"
   if(Noutspc>0) then
      write(22,'(a)',advance='no') "plot 'plt.dat' u 1:3 w l title '"//trim(SYM_SPC(outspc(1)))//"'"
      do i=2,Noutspc
         write(22,'(a,i4,a)',advance='no') ",'' u 1:",i+2," w l title '"//trim(SYM_SPC(outspc(i)))//"'"
      end do
      write(22,'()')
      write(22,'(a)') "pause -1"
   end if
   close(22)
end subroutine out_plotplt

program driver
   use chem
   use chem_var
   implicit none
   double precision rho,T
   double precision dt,tout
   double precision Y(2)
   double precision vrho(max_ns)
   integer          outspc(max_ns)
   integer          Noutspc

   !read files
   call init_therm

   !calc stoichiometry
   Y(1)=1d0/(1d0+of)
   Y(2)=1d0-Y(1)
   rho =1d0/(Y(1)/rhof+Y(2)/rhoo)
   vrho=rho*(vwf*Y(1)+vwo*Y(2))
   T   =      Tf*Y(1)+ To*Y(2)

   call read_controlinp(T,vrho,dt,tout,outspc,Noutspc)
   print *,"(tout,dt)=",tout,dt

   !plot
   call plot_time_history(T,vrho,dt,tout,outspc,Noutspc)
   call out_plotplt(outspc,Noutspc)
   print *,"done"
end program driver

