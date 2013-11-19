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
   double precision temp_mat((nimax+1)*(njmax+1))
   integer i,nr,nall,plane

   if(myid .eq. 0) then
      open(77,file=trim(file_coordinate))

      !read the number of plane and the number of nodes
      read(77,*) t1
      if(t1 .ne. Nplane) stop "The number of plane is odd"


      do plane = 1,Nplane
         read(77,*) t1,t2,t3
         if(t1 .ne. ni(plane)+1 .or. t2 .ne. nj(plane)+1 .or. t3 .ne. 1) stop "The number of node is odd."
      end do

      do plane = 1,Nplane
         !calc the number of nodes
         nall=(ni(plane)+1)*(nj(plane)+1)
         nr=nall/4

         !read x_{i+1/2}
         do i=1,nr
            read(77,*) temp_mat((i-1)*4+1:i*4)
         end do
         if(nr*4 .ne. nall) read(77,*) temp_mat(nr*4+1:nall)
         xh(0:ni(plane),0:nj(plane),plane)=reshape(temp_mat(1:nall),(/ni(plane)+1,nj(plane)+1/))

         !read r_{i+1/2}
         do i=1,nr
            read(77,*) temp_mat((i-1)*4+1:i*4)
         end do
         if(nr*4 .ne. nall) read(77,*) temp_mat(nr*4+1:nall)
         rh(0:ni(plane),0:nj(plane),plane)=reshape(temp_mat(1:nall),(/ni(plane)+1,nj(plane)+1/))

         !read z_{i+1/2} (now throw away)
         do i=1,nr
            read(77,*) temp_mat((i-1)*4+1:i*4)
         end do
         if(nr*4 .ne. nall) read(77,*) temp_mat(nr*4+1:nall)
      end do

      close(77)
   end if

   call MPI_Bcast(xh(0,0,1), 2*(nimax+1)*(njmax+1)*Nplane, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(rh(0,0,1), 2*(nimax+1)*(njmax+1)*Nplane, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end subroutine init_xrh!}}}
subroutine set_xr!{{{
   implicit none
   double precision nx,nr
   integer i,j,plane

   do plane = 1,Nplane
      do j=1,nj(plane)
         do i=1,ni(plane)
            x(i,j,plane)=1d0/4d0*(xh(i-1,j-1,plane)+xh(i,j-1,plane)+xh(i-1,j,plane)+xh(i,j,plane))
            r(i,j,plane)=1d0/4d0*(rh(i-1,j-1,plane)+rh(i,j-1,plane)+rh(i-1,j,plane)+rh(i,j,plane))
         end do
      end do

      do j=1,nj(plane)
         x(   0,j,plane)=2d0*x( 1,j,plane)-x(   2,j,plane)
         r(   0,j,plane)=2d0*r( 1,j,plane)-r(   2,j,plane)
         x(ni+1,j,plane)=2d0*x(ni,j,plane)-x(ni-1,j,plane)
         r(ni+1,j,plane)=2d0*r(ni,j,plane)-r(ni-1,j,plane)
      end do

      do i=0,ni(plane)+1
         x(i,   0,plane)=2d0*x(i, 1,plane)-x(i,   2,plane)
         r(i,   0,plane)=2d0*r(i, 1,plane)-r(i,   2,plane)
         x(i,nj+1,plane)=2d0*x(i,nj,plane)-x(i,nj-1,plane)
         r(i,nj+1,plane)=2d0*r(i,nj,plane)-r(i,nj-1,plane)
      end do
   end do
end subroutine set_xr!}}}
subroutine set_Area_Vol!{{{
   implicit none
   integer i,j,plane
   double precision Aabd,Abcd

   do plane = 1,Nplane
      do j=1,nj(plane)
         do i=1,ni(plane)
            Aabd=Area_Triangle(xh(i-1,j-1,plane),rh(i-1,j-1,plane),xh(i,j-1,plane),rh(i,j-1,plane),xh(i-1,j,plane),rh(i-1,j,plane))
            Abcd=Area_Triangle(xh(  i,j-1,plane),rh(  i,j-1,plane),xh(i,  j,plane),rh(i,  j,plane),xh(i-1,j,plane),rh(i-1,j,plane))
            Area(i,j,plane)=Aabd+Abcd
            Vol( i,j,plane)=(rh(i-1,j-1,plane)+rh(i,j-1,plane)+rh(i-1,j,plane))*Aabd/3d0&
                           +(rh(  i,j-1,plane)+rh(i,  j,plane)+rh(i-1,j,plane))*Abcd/3d0
         end do
      end do
   end do
end subroutine set_Area_Vol!}}}
subroutine set_dsij!{{{
   implicit none
   integer i,j,plane

   do plane = 1,Nplane
      do j=0,nj(plane)
         do i=1,ni(plane)
            dsj(i,j,plane)=dis_point(xh(i-1,j,plane),rh(i-1,j,plane),xh(i,j,plane),rh(i,j,plane))*(rh(i-1,j,plane)+rh(i,j,plane))/2d0
         end do
      end do

      do j=1,nj(plane)
         do i=0,ni(plane)
            dsi(i,j,plane)=dis_point(xh(i,j-1,plane),rh(i,j-1,plane),xh(i,j,plane),rh(i,j,plane))*(rh(i,j-1,plane)+rh(i,j,plane))/2d0
         end do
      end do
   end do
end subroutine set_dsij!}}}
subroutine set_vn!{{{
   implicit none
   double precision nx,nr,std
   integer i,j,plane

   do plane=1,Nplane
      do j=1,nj(plane)
         do i=0,ni(plane)
           nx=  rh(i,j,plane)-rh(i,j-1,plane)
           nr=-(xh(i,j,plane)-xh(i,j-1,plane))
           std=sqrt(nx**2+nr**2)

           vni(1,i,j,plane)=nx/std
           vni(2,i,j,plane)=nr/std
         end do
      end do

      do j=0,nj(plane)
         do i=1,ni(plane)
           nx=-(rh(i,j,plane)-rh(i-1,j,plane))
           nr=  xh(i,j,plane)-xh(i-1,j,plane)
           std=sqrt(nx**2+nr**2)

           vnj(1,i,j,plane)=nx/std
           vnj(2,i,j,plane)=nr/std
         end do
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
   integer i,j,plane

   open(45,file="geometry.bin",form="unformatted")
   do plane = 1,Nplane
      write(45) ni(plane)
      write(45) nj(plane)
      write(45) ((xh(i,j,plane),rh(i,j,plane),i=0,ni(plane)),j=0,nj(plane))
   end do
   close(45)
end subroutine out_geo_bin!}}}
end module geometry
