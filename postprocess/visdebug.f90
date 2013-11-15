program main
   use prmtr
   implicit none
   double precision,allocatable,dimension(:,:)::x,y
   double precision,allocatable,dimension(:,:,:)::q
   integer nx,ny
   integer i,j,k

   !read geometry binary{{{
   open(55,file="geometry.bin",form="unformatted")
   read(55) nx
   read(55) ny
   allocate(x(nx,ny))
   allocate(y(nx,ny))
   read(55) ((x(i,j),y(i,j),i=1,nx),j=1,ny)
   close(55)
   !}}}

   !read binary
   !allocate(q(nx,ny))
   allocate(q(2,nx,ny))
   open(55,file="debug.bin",form="unformatted")
   !read(55) ((q(i,j),i=1,nx),j=1,ny)
   read(55) (((q(i,j,k),i=1,2),j=1,nx),k=1,ny)
   close(55)

   !output file open
   open(66,file="debug.vtk")

   !header{{{
   write(66,'(a26)') '# vtk DataFile Version 2.0'
   write(66,'(a7)') '2D Data'
   write(66,'(a5)') 'ASCII'
   write(66,'(a23)') 'DATASET STRUCTURED_GRID'
   write(66,'(a11,i3,a1,i3,a1,i1)') 'DIMENSIONS ', nx, ' ', ny, ' ', 1
   !}}}

   !grid{{{
   write(66,'(a7,i7,a6)') 'POINTS ', nx*ny*1, ' float'
   do j = 1, ny
     do i = 1, nx
       write(66,'(e14.6e3,2(1x,e14.6e3))') x(i,j), y(i,j), 0.d0
     end do
   end do
   write(66,*)
 
   write(66,'(a11,i7)') 'POINT_DATA ', nx*ny*1
   !}}}

   !!scalar
   !write(66,'(a23)') 'SCALARS DEBUG float 1'
   !write(66,'(a20)') 'LOOKUP_TABLE default'
   !do j = 1, ny
   !  do i = 1, nx
   !    write(66,'(e14.6e3)') q(i,j)
   !  end do
   !end do
   !write(66,*)

   !velocity
   write(66,'(a22)') 'VECTORS DEBUG float'
   do j = 1, ny
     do i = 1, nx
       write(66,'(e14.6e3,2(1x,e14.6e3))') q(1,i,j), q(2,i,j), 0.d0
     end do
   end do
   write(66,*)

   !output file close
   close(66)
end program main
