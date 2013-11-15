program main
   use grbl_prmtr
   use prmtr
   implicit none
   double precision,allocatable,dimension(:,:)::x,y,rho_mat,T_mat,u_mat,v_mat,M_mat,p_mat,R_gas,gmma
   double precision,allocatable,dimension(:,:,:)::q,w_mat
   double precision rho,u,v,T,tmp
   integer nx,ny
   integer io,jo
   integer step,step_bin
   integer i,j,k

   double precision x_obj,y_obj
   double precision current_min

   !integer,parameter::Dstep=1000*1000
   integer,parameter::Dstep=10*1000
   character*50 filename
   !character*50,parameter::dirname="../CH4Air/"
   character*50,parameter::dirname=""
   integer*4,external::access

   !read geometry binary
   open(55,file=trim(dirname)//"geometry.bin",form="unformatted")
   read(55) nx
   read(55) ny
   allocate(x(nx,ny))
   allocate(y(nx,ny))
   read(55) ((x(i,j),y(i,j),i=1,nx),j=1,ny)
   close(55)

   !allocate other variables
   allocate(    q(dimq,nx,ny))
   allocate(w_mat(dimq,nx,ny))
   allocate(rho_mat(nx,ny))
   allocate(  T_mat(nx,ny))
   allocate(  u_mat(nx,ny))
   allocate(  v_mat(nx,ny))
   allocate(  M_mat(nx,ny))
   allocate(  p_mat(nx,ny))
   allocate(  R_gas(nx,ny))
   allocate(  gmma(nx,ny))

   call getarg(1,filename)
   if(trim(filename) .eq. "") then
      step=0
   else
      print *,filename
      read(filename,'(i12.0)') step
      !step=step-Dstep
   end if

   !set objective nx and ny
   x_obj=3d-2
   y_obj=0.5d-2
   current_min = 1d300
   do i=1,nx
      do j=1,ny
         tmp=(x(i,j)-x_obj)**2+(y(i,j)-y_obj)**2
         if(current_min>tmp) then
            current_min = tmp
            io = i
            jo = j
         end if
      end do
   end do

   !print infomation of graphed point.
   print '(a,es9.2,a,es9.2,a)',"The nearist point to (",x_obj,",",y_obj,") is"
   print '(a,es9.2,a,es9.2,a,i4,a,i4)',"(",x(io,jo),",",y(io,jo),") at i=",io," j=",jo

   !open output file
   open(50,file="flow_history.dat")
   
   do
      step=step+Dstep
      write(filename,'("result/result",i12.12)') step
      write(*       ,'("result/result",i12.12)') step
      if(access(trim(dirname)//trim(filename)//".bin","r") .ne. 0) exit

      !read q binary
      open(55,file=trim(dirname)//trim(filename)//".bin",form="unformatted")
      read(55) step_bin
      read(55) (((q(i,j,k),i=1,dimq),R_gas(j,k),gmma(j,k),p_mat(j,k),j=1,nx),k=1,ny)
      close(55)

      !trim data
      do j=1,ny
         do i=1,nx

            rho=q(1,i,j)
            u  =q(2,i,j)/rho
            v  =q(3,i,j)/rho
            T  =p_mat(i,j)/(rho*R_gas(i,j))

            rho_mat(i,j)=rho
              T_mat(i,j)=T
              u_mat(i,j)=u
              v_mat(i,j)=v
              M_mat(i,j)=sqrt((u**2+v**2)/(gmma(i,j)*R_gas(i,j)*T))

            do k=5,dimq
               w_mat(k,i,j)=q(k,i,j)/rho
            end do
         end do
      end do

      write(50,*) step,u_mat(io,jo),T_mat(io,jo)
   end do
   close(50)
end program main
