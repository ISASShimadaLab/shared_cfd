program main
   implicit none
   double precision,allocatable,dimension(:,:)::x,y,p_mat,R_gas,gmma
   double precision,allocatable,dimension(:,:,:)::q
   double precision,allocatable,dimension(:)::qmin
   integer step_bin
   integer nx,ny
   integer i,j,k
   integer ii,jj
   double precision rho, Yf, ei, eil, eih, dis, part
   double precision,parameter::R_uni=8.314472d3
   integer,parameter::numY=22
   integer,parameter::dimq=numY+3
   double precision,parameter::pi = 2d0*acos(0d0)
   character*50 filename

   !read geometry binary
   open(55,file="geometry.bin",form="unformatted")
   read(55) nx
   read(55) ny
   allocate(x(nx,ny))
   allocate(y(nx,ny))
   read(55) ((x(i,j),y(i,j),i=1,nx),j=1,ny)
   close(55)

   !allocate other variables
   allocate(   q(dimq,nx,ny))
   allocate(    R_gas(nx,ny))
   allocate(    gmma( nx,ny))
   allocate(    p_mat(nx,ny))

   allocate(qmin(dimq))

   call getarg(1,filename)

   !read q binary
   open(55,file=trim(filename),form="unformatted")
   read(55) step_bin
   read(55) (((q(i,j,k),i=1,dimq),R_gas(j,k),gmma(j,k),p_mat(j,k),j=1,nx),k=1,ny)
   close(55)


   !!!!!!!!!!!!!!!!!!!!!!   PROCESSING   !!!!!!!!!!!!!!!!!!!!!!
   qmin(:)=0d0
   do j=1,ny
      do i=1,nx
         dis = sqrt(x(i,j)**2 + y(i,j)**2 )
         if(dis <5d-3) then
            rho = 0d0
            do k=1,3
               rho = rho + q(k,i,j)
            end do
            Yf=q(1,i,j)/rho

            ei =  q(dimq,i,j)/rho -0.5d0*(q(dimq-1,i,j)**2+q(dimq-2,i,j)**2)/rho**2

            eih = Yf*(-4.5024d5) +(1d0-Yf)*( 9.1413d5)
            part = 0.5d0 * ( cos(pi *dis/5d-3) + 1d0)
            if(eih <ei) then
               print *,"ERROR!!!!"
               call exit(1)
            end if
            ei  = (eih-ei)*part + eil

            q(dimq,i,j) =  rho*ei+ 0.5d0*(q(dimq-1,i,j)**2+q(dimq-2,i,j)**2)/rho
            p_mat(i,j) = p_mat(i,j) / 300d0 * ( 1200d0 * part + 300d0)
         end if
      end do
   end do
   !!!!!!!!!!!!!!!!!!!!!!   END PROCESSING   !!!!!!!!!!!!!!!!!!

   !output file open
   open(66,file="modified.bin",form="unformatted")
   write(66) step_bin
   write(66) (((q(i,j,k),i=1,dimq),R_gas(j,k),gmma(j,k),p_mat(j,k),j=1,nx),k=1,ny)

   !output file close
   close(66)
end program main
