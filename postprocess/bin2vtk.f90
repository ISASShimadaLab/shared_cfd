program main
   use prmtr
   implicit none
   double precision,allocatable,dimension(:,:)::x,y,rho_mat,T_mat,u_mat,v_mat,M_mat,p_mat
   double precision,allocatable,dimension(:,:,:)::q
   double precision rho,u,v,T,tmp
   integer nx,ny
   integer step
   integer i,j,k
   character*500 filename
   integer*4,external::access

   !read geometry binary{{{
   open(55,file="geometry.bin",form="unformatted")
   read(55) nx
   read(55) ny
   allocate(x(nx,ny))
   allocate(y(nx,ny))
   read(55) ((x(i,j),y(i,j),i=1,nx),j=1,ny)
   close(55)
   !}}}

   !allocate other variables{{{
   allocate(q(4,nx,ny))
   allocate(rho_mat(nx,ny))
   allocate(  T_mat(nx,ny))
   allocate(  u_mat(nx,ny))
   allocate(  v_mat(nx,ny))
   allocate(  M_mat(nx,ny))
   allocate(  p_mat(nx,ny))
   !}}}

   !read q binary{{{
   call getarg(1,filename)
   print *,"read ",trim(filename)
   open(55,file=trim(filename),form="unformatted")
   read(55) (((q(i,j,k),i=1,4),j=1,nx),k=1,ny)
   close(55)
   !}}}

   !trim data{{{
   do j=1,ny
      do i=1,nx

         rho=q(1,i,j)
         u  =q(2,i,j)/rho
         v  =q(3,i,j)/rho
         tmp=q(4,i,j)/rho
         tmp=tmp-0.5d0*(u**2+v**2)
         T  =tmp/(1d0/(gmma-1d0)*R_gas)

         rho_mat(i,j)=rho
           T_mat(i,j)=T
           u_mat(i,j)=u
           v_mat(i,j)=v
           M_mat(i,j)=sqrt((u**2+v**2)/(gmma*R_gas*T))
           p_mat(i,j)=rho*R_gas*T
      end do
   end do
   !}}}

   !output file open
   open(66,file="bin2vtk.vtk")

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

   !thermodynamical{{{
   !density{{{ 
   write(66,'(a23)') 'SCALARS Density float 1'
   write(66,'(a20)') 'LOOKUP_TABLE default'
   do j = 1, ny
     do i = 1, nx
       write(66,'(e14.6e3)') rho_mat(i,j)
     end do
   end do
   write(66,*)
   !}}}

   !Temperature{{{ 
   write(66,'(a27)') 'SCALARS Temperature float 1'
   write(66,'(a20)') 'LOOKUP_TABLE default'
   do j = 1, ny
     do i = 1, nx
       write(66,'(e14.6e3)') T_mat(i,j)
     end do
   end do
   write(66,*)
   !}}}

   !Pressure{{{ 
   write(66,'(a)') 'SCALARS Pressure float 1'
   write(66,'(a20)') 'LOOKUP_TABLE default'
   do j = 1, ny
     do i = 1, nx
       write(66,'(e14.6e3)') P_mat(i,j)
     end do
   end do
   write(66,*)
   !}}}
   !}}}

   !velocity{{{
   !Mach number{{{ 
   write(66,'(a)') 'SCALARS Machnumber float 1'
   write(66,'(a20)') 'LOOKUP_TABLE default'
   do j = 1, ny
     do i = 1, nx
       write(66,'(e14.6e3)') M_mat(i,j)
     end do
   end do
   write(66,*)
   !}}}

   !velocity{{{
   write(66,'(a22)') 'VECTORS Velocity float'
   do j = 1, ny
     do i = 1, nx
       write(66,'(e14.6e3,2(1x,e14.6e3))') u_mat(i,j), v_mat(i,j), 0.d0
     end do
   end do
   write(66,*)
   !}}}
   !}}}

   !output file close
   close(66)
end program main
