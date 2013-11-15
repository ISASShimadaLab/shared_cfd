module geometry
   use grbl_prmtr
   use variable
   implicit none
contains
subroutine init_geometry!{{{
   use mod_mpi
   implicit none

   if(myid .eq. 0) print *,"Initializing geometry matrix..."

   call init_xrh

   call set_xr

   call set_Area_Vol

   call set_dsij

   call set_vn

   if(myid .eq. 0) call out_geo_bin

   if(myid .eq. 0) print *,"Initialized geometry matrix."
end subroutine init_geometry!}}}

subroutine init_xrh!{{{
   use mod_mpi
   implicit none
   double precision t1,t2,t3
   double precision temp_mat((ni+1)*(nj+1))
   integer i,nr,nall

   if(myid .eq. 0) then
      open(77,file=trim(file_coordinate))

      !read the number of plane and the number of nodes
      read(77,*) t1
      !if(t1 .ne. 1) stop "Not one geometry plane. t1=",t1
      if(t1 .ne. 1) then
         print *,"t1=",t1
         stop "Not one geometry plane."
      end if
      read(77,*) t1,t2,t3
      if(t1 .ne. ni+1 .or. t2 .ne. nj+1 .or. t3 .ne. 1) stop "The number of node is odd."

      !calc the number of nodes
      nall=(ni+1)*(nj+1)
      nr=nall/4

      !read x_{i+1/2}
      do i=1,nr
         read(77,*) temp_mat((i-1)*4+1:i*4)
      end do
      if(nr*4 .ne. nall) read(77,*) temp_mat(nr*4+1:nall)
      xh(0:ni,0:nj)=reshape(temp_mat,(/ni+1,nj+1/))
      xh(:,:)=xh(:,:)

      !read r_{i+1/2}
      do i=1,nr
         read(77,*) temp_mat((i-1)*4+1:i*4)
      end do
      if(nr*4 .ne. nall) read(77,*) temp_mat(nr*4+1:nall)
      rh(0:ni,0:nj)=reshape(temp_mat,(/ni+1,nj+1/))

      close(77)
   end if

   call MPI_Bcast(xh(0,0), 2*(ni+1)*(nj+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(rh(0,0), 2*(ni+1)*(nj+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end subroutine init_xrh!}}}
subroutine set_xr!{{{
   implicit none
   double precision nx,nr
   integer i,j

   do j=1,nj
      do i=1,ni
         x(i,j)=1d0/4d0*(xh(i-1,j-1)+xh(i,j-1)+xh(i-1,j)+xh(i,j))
         r(i,j)=1d0/4d0*(rh(i-1,j-1)+rh(i,j-1)+rh(i-1,j)+rh(i,j))
      end do
   end do

   do j=1,nj
      x(   0,j)=2d0*x( 1,j)-x(   2,j)
      r(   0,j)=2d0*r( 1,j)-r(   2,j)
      x(ni+1,j)=2d0*x(ni,j)-x(ni-1,j)
      r(ni+1,j)=2d0*r(ni,j)-r(ni-1,j)
   end do

   do i=0,ni+1
      x(i,   0)=2d0*x(i, 1)-x(i,   2)
      r(i,   0)=2d0*r(i, 1)-r(i,   2)
      x(i,nj+1)=2d0*x(i,nj)-x(i,nj-1)
      r(i,nj+1)=2d0*r(i,nj)-r(i,nj-1)
   end do
end subroutine set_xr!}}}
subroutine set_Area_Vol!{{{
   implicit none
   integer i,j
   double precision Aabd,Abcd

   do j=1,nj
      do i=1,ni
         Aabd=Area_Triangle(xh(i-1,j-1),rh(i-1,j-1),xh(i,j-1),rh(i,j-1),xh(i-1,j),rh(i-1,j))
         Abcd=Area_Triangle(xh(  i,j-1),rh(  i,j-1),xh(i,  j),rh(i,  j),xh(i-1,j),rh(i-1,j))
         Vol(i,j)=Aabd+Abcd
      end do
   end do
end subroutine set_Area_Vol!}}}
subroutine set_dsij!{{{
   implicit none
   integer i,j

   do j=0,nj
      do i=1,ni
         dsj(i,j)=dis_point(xh(i-1,j),rh(i-1,j),xh(i,j),rh(i,j))
      end do
   end do

   do j=1,nj
      do i=0,ni
         dsi(i,j)=dis_point(xh(i,j-1),rh(i,j-1),xh(i,j),rh(i,j))
      end do
   end do
end subroutine set_dsij!}}}
subroutine set_vn!{{{
   implicit none
   double precision nx,nr,std
   integer i,j

   do j=1,nj
      do i=0,ni
        nx=  rh(i,j)-rh(i,j-1)
        nr=-(xh(i,j)-xh(i,j-1))
        std=sqrt(nx**2+nr**2)

        vni(1,i,j)=nx/std
        vni(2,i,j)=nr/std
      end do
   end do

   do j=0,nj
      do i=1,ni
        nx=-(rh(i,j)-rh(i-1,j))
        nr=  xh(i,j)-xh(i-1,j)
        std=sqrt(nx**2+nr**2)

        vnj(1,i,j)=nx/std
        vnj(2,i,j)=nr/std
      end do
   end do
end subroutine set_vn!}}}
double precision function Area_Triangle(a1,a2,b1,b2,c1,c2)!{{{
   implicit none
   double precision,intent(in)::a1,a2,b1,b2,c1,c2

   Area_Triangle=1d0/2d0*(a1*(b2-c2)+b1*(c2-a2)+c1*(a2-b2))
end function Area_Triangle!}}}
double precision function dis_point(x1,y1,x2,y2)!{{{
   implicit none
   double precision,intent(in)::x1,y1,x2,y2

   dis_point=sqrt((x1-x2)**2+(y1-y2)**2)
end function dis_point!}}}
subroutine out_geo_bin!{{{
   implicit none
   integer i,j

   open(45,file="geometry.bin",form="unformatted")
   write(45) ni
   write(45) nj
   write(45) ((x(i,j),r(i,j),i=1,ni),j=1,nj)
   close(45)
end subroutine out_geo_bin!}}}
end module geometry
