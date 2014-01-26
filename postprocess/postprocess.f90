module to_vis
   use grbl_prmtr
   implicit none
   integer,dimension(Nplane)::nxs_vis,nxe_vis,nys_vis,nye_vis
   character*50 dirname,filename
   integer Dstep,step,OverWrite
   double precision tt

   integer offset,ng
end module to_vis

program main
   use to_vis
   implicit none
   double precision,dimension(     0:nimax+1,0:njmax+1,     Nplane)::xh,rh
   double precision,dimension(dimq,0:nimax+1,0:njmax+1,     Nplane)::q
   double precision,dimension(     0:nimax+1,0:njmax+1,     Nplane)::rho_mat,T_mat,u_mat,v_mat,p_mat
   double precision,dimension(     0:nimax+1,0:njmax+1,     Nplane)::M_mat,R_gas,gmma,w_res,ei_mat,MW_mat
   double precision,dimension(     0:nimax+1,0:njmax+1,dimq,Nplane)::w_mat
   !double precision,dimension(     0:nimax+1,0:njmax+1,   2,Nplane)::debug
   double precision rho,u,v,T
   integer nx_bin,ny_bin,step_bin
   integer i,j,k,plane
   double precision,parameter::R_uni=8.314472d3
   character*50 title,str_num
   integer*4,external::access

   double precision tmp

   !read control parameters
   call read_control

   !read geometry binary
   open(55,file=trim(dirname)//"geometry.bin",form="unformatted")
   do plane = 1,Nplane
      read(55) nx_bin
      read(55) ny_bin
      read(55) ((xh(i,j,plane),rh(i,j,plane),i=0,ni(plane)),j=0,nj(plane))
   end do
   close(55)

   !process time history
   do
      write(filename,'("result/result",i12.12)') step
      if(access(trim(dirname)//trim(filename)//".bin","r") .ne. 0) exit
      write(*       ,'(a)') filename
      write(filename,'("result/result001.",i12.12)') step
      if(OverWrite .eq. 0 .and. access(trim(filename)//".vts","r") .eq. 0) then
         step=step+Dstep
         cycle
      end if
      write(filename,'("result/result",i12.12)') step

      !read q binary
      open(55,file=trim(dirname)//trim(filename)//".bin",form="unformatted")
      read(55) step_bin
      read(55) tt
      read(55) ((((q(i,j,k,plane),i=1,dimq),R_gas(j,k,plane),gmma(j,k,plane),p_mat(j,k,plane),j=0,ni(plane)+1),k=0,nj(plane)+1),plane=1,Nplane)
      close(55)

      !!read q binary
      !open(55,file=trim(dirname)//trim(filename)//".debug.bin",form="unformatted")
      !read(55) ((((debug(j,k,i,plane),i=1,2),j=0,ni(plane)+1),k=0,nj(plane)+1),plane=1,Nplane)
      !close(55)

      !Process data
      do plane=1,Nplane
         do j=0,nj(plane)+1
            do i=0,ni(plane)+1
               rho=0d0
               do k=1,nY
                  rho=rho+q(k,i,j,plane)
               end do
               u  =q(nY+1,i,j,plane)/rho
               v  =q(nY+2,i,j,plane)/rho
               T  =p_mat( i,j,plane)/rho/R_gas(i,j,plane)

               rho_mat(i,j,plane)=rho
                 T_mat(i,j,plane)=T
                 u_mat(i,j,plane)=u
                 v_mat(i,j,plane)=v
                 M_mat(i,j,plane)=sqrt((u**2+v**2) / (gmma(i,j,plane)*p_mat(i,j,plane)/rho))
                ei_mat(i,j,plane)=q(nY+3,i,j,plane)/rho -0.5d0*(u**2+v**2)
                MW_mat(i,j,plane)=R_uni/R_gas(i,j,plane)

               w_res(i,j,plane)=0d0
               do k=1,nY
                  w_mat(i,j,k,plane)=q(k,i,j,plane)/rho
                  if(k>3) then
                     w_res(i,j,plane)=w_res(i,j,plane)+w_mat(i,j,k,plane)
                  end if
               end do
            end do
         end do
      end do

      !calc half values
      do plane=1,Nplane
         call center_to_grid(plane,rho_mat(:,:,plane))
         call center_to_grid(plane,u_mat(:,:,plane))
         call center_to_grid(plane,v_mat(:,:,plane))
         call center_to_grid(plane,p_mat(:,:,plane))
         call center_to_grid(plane,gmma(:,:,plane))
         call center_to_grid(plane,T_mat(:,:,plane))
         call center_to_grid(plane,MW_mat(:,:,plane))
         call center_to_grid(plane,ei_mat(:,:,plane))

         call center_to_grid(plane,M_mat(:,:,plane))
         do k=1,nY
            call center_to_grid(plane,w_mat(:,:,k,plane))
         end do
         call center_to_grid(plane,w_res(:,:,plane))
      end do

      !write(filename,"('result/result',i3.3,'.',i12.12,'.plt')") plane,step
      !!open(28,file=trim(filename))
      !open(28,file="test.dat")
      !tmp=0d0
      !j=0
      !i=0
      !write(28,*) tmp,p_mat(i,j,1)/p_mat(i,0,1)
      !do i=1,ni(1)
      !   tmp=tmp+sqrt((xh(i,j,1)-xh(i-1,j,1))**2+(rh(i,j,1)-rh(i-1,j,1))**2)
      !   write(28,*) tmp/2d0/rh(ni(1),0,1),p_mat(i,j,1)/p_mat(0,j,1)
      !end do
      !close(28)

      do plane = 1,Nplane
         !output file open
         write(filename,"('result/result',i3.3,'.',i12.12,'.vts')") plane,step

         !header
         call vtk_header(plane)

         call vtk_scalar("Density")
         call vtk_scalar("Velocity_x")
         call vtk_scalar("Velocity_y")
         call vtk_scalar("Pressure")
         call vtk_scalar("SpecificHeatRatio")
         call vtk_vector("Velocity")
         call vtk_scalar("Temperature")
         call vtk_scalar("Molecular_Weight")
         call vtk_scalar("InternalEnergy")

         call vtk_scalar("MachNumber")

         do k=1,min(3,nY)
            write(title,'(a,i3.3)') 'phi',k
            call vtk_scalar(title)
         end do

         call vtk_scalar("phi_res")
         !call vtk_scalar("debug1")
         !call vtk_scalar("debug2")

         !binary
         call vtk_header_bin(plane,xh(:,:,plane),rh(:,:,plane))

         call vtk_scalar_bin(plane,rho_mat(:,:,plane))
         call vtk_scalar_bin(plane,u_mat(:,:,plane))
         call vtk_scalar_bin(plane,v_mat(:,:,plane))
         call vtk_scalar_bin(plane,p_mat(:,:,plane))
         call vtk_scalar_bin(plane,gmma(:,:,plane))
         call vtk_vector_bin(plane,u_mat(:,:,plane),v_mat(:,:,plane))
         call vtk_scalar_bin(plane,T_mat(:,:,plane))
         call vtk_scalar_bin(plane,MW_mat(:,:,plane))
         call vtk_scalar_bin(plane,ei_mat(:,:,plane))

         call vtk_scalar_bin(plane,M_mat(:,:,plane))

         do k=1,min(3,nY)
            call vtk_scalar_bin(plane,w_mat(:,:,k,plane))
         end do

         call vtk_scalar_bin(plane,w_res(:,:,plane))

         !call vtk_scalar_bin(plane,debug(:,:,1,plane))
         !call vtk_scalar_bin(plane,debug(:,:,2,plane))

         !output file close
         call vtk_footer_bin
      end do

      step=step+Dstep
   end do
contains
   subroutine vtk_header(plane)!{{{
      implicit none
      integer,intent(in)::plane

      ng=(nxe_vis(plane)-nxs_vis(plane)+1)*(nye_vis(plane)-nys_vis(plane)+1)

      open(55,file=trim(filename))
      write(55,'(a)') '<?xml version="1.0"?>'
      write(55,'(a)') '<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian">'
      write(55,'(a,6i4,a)') '<StructuredGrid WholeExtent="',&
                             nxs_vis(plane),&
                             nxe_vis(plane),&
                             nys_vis(plane),&
                             nye_vis(plane),&
                             0,0,'">'

      write(55,'(a)') '<FieldData>'
      write(55,'(a,e12.5,a)') '<DataArray type="Float32" Name="TIME" &
                           &NumberOfTuples="1" format="ascii">',&
                           tt,&
                           '</DataArray>'
      write(55,'(a,i12,a)') '<DataArray type="Int32" Name="CYCLE" &
                           &NumberOfTuples="1" format="ascii">',&
                           step,&
                           '</DataArray>'
      write(55,'(a)') '</FieldData>'

      write(55,'(a,6i4,a)') '<Piece Extent="',&
                             nxs_vis(plane),&
                             nxe_vis(plane),&
                             nys_vis(plane),&
                             nye_vis(plane),&
                            0,0,'">'
   
      offset=0
      write(55,'(a)') "<Points>"
      write(55,'(a,i8,a,es12.5,a,es12.5,a)') &
                      '<DataArray type="Float32" &
                      &Name="xyz" &
                      &NumberOfComponents="3" &
                      &format="appended" &
                      &offset="',offset,'">'
      write(55,'(a)') "</DataArray>"
      offset=offset+ng*4*3+4
      write(55,'(a)') "</Points>"
   
      write(55,'(a)') "<PointData>"

      ! Indices Vector
      write(55,'(a,i8,a,es12.5,a,es12.5,a)') &
                      '<DataArray type="Float32" &
                      &Name="Indices" &
                      &NumberOfComponents="3" &
                      &format="appended" &
                      &offset="',offset,'">'
      write(55,'(a)') "</DataArray>"
      offset=offset+ng*4*3+4
   end subroutine vtk_header!}}}
   subroutine vtk_header_bin(plane,xh,rh)!{{{
      implicit none
      integer,intent(in)::plane
      double precision,dimension(0:nimax+1,0:njmax+1),intent(in)::xh,rh
      integer i,j

      write(55,'(a)') "</PointData>"
   
      write(55,'(a)') "</Piece>"
      write(55,'(a)') "</StructuredGrid>"
      write(55,'(a)') '<AppendedData encoding="raw">'
   
      close(55)
   
      open(55,file=trim(filename),form='binary',recl=4,access='append')
      write(55) '_'
   
      write(55) 4*ng*3
      do j = nys_vis(plane), nye_vis(plane)
        do i = nxs_vis(plane), nxe_vis(plane)
            write(55) real(xh(i,j))
            write(55) real(rh(i,j))
            write(55) 0e0
         end do
      end do
 
      write(55) 4*ng*3
      do j = nys_vis(plane), nye_vis(plane)
        do i = nxs_vis(plane), nxe_vis(plane)
            write(55) real(i)
            write(55) real(j)
            write(55) 0e0
         end do
      end do
   end subroutine vtk_header_bin!}}}
   subroutine vtk_footer_bin!{{{
      implicit none
      close(55)
   
      open(55,file=trim(filename),access='append')
      write(55,'()')
      write(55,'(a)') "</AppendedData>"
      write(55,'(a)') "</VTKFile>"
      close(55)
   end subroutine vtk_footer_bin!}}}
   subroutine vtk_scalar(title)!{{{
      implicit none
      character(*),intent(in)::title

      write(55,'(a,i8,a,es12.5,a,es12.5,a)') &
                      '<DataArray type="Float32" &
                      &Name="'//trim(title)//'" &
                      &format="appended" &
                      &offset="',offset,'">'
      write(55,'(a)') "</DataArray>"
      offset=offset+ng*4+4;
   end subroutine vtk_scalar!}}}
   subroutine vtk_vector(title)!{{{
      implicit none
      character(*),intent(in)::title

      write(55,'(a,i8,a,es12.5,a,es12.5,a)') &
                      '<DataArray type="Float32" &
                      &Name="'//trim(title)//'" &
                      &NumberOfComponents="3" &
                      &format="appended" &
                      &offset="',offset,'">'
      write(55,'(a)') "</DataArray>"
      offset=offset+ng*4*3+4
   end subroutine vtk_vector!}}}
   subroutine vtk_scalar_bin(plane,arr)!{{{
      implicit none
      integer,intent(in)::plane
      double precision,dimension(0:nimax+1,0:njmax+1),intent(in)::arr

      integer i,j

      write(55) 4*ng
      do j = nys_vis(plane), nye_vis(plane)
        do i = nxs_vis(plane), nxe_vis(plane)
          write(55) real(arr(i,j))
        end do
      end do
   end subroutine vtk_scalar_bin!}}}
   subroutine vtk_vector_bin(plane,arr1,arr2)!{{{
      implicit none
      integer,intent(in)::plane
      double precision,dimension(0:nimax+1,0:njmax+1),intent(in)::arr1,arr2

      integer i,j

      write(55) 4*ng*3
      do j = nys_vis(plane), nye_vis(plane)
        do i = nxs_vis(plane), nxe_vis(plane)
          write(55) real(arr1(i,j))
          write(55) real(arr2(i,j))
          write(55) 0e0
        end do
      end do
   end subroutine vtk_vector_bin!}}}

   subroutine read_control!{{{
      implicit none
      integer i

      open(8,file="control.inp")

      read(8,'()')
      read(8,'(25x,i10)') Dstep
      read(8,'(25x,i10)') step
      read(8,'(25x,a)')   dirname
      read(8,'(25x,i10)') OverWrite
      do i=1,Nplane
         read(8,'(25x,i10)') nxs_vis(i)
         read(8,'(25x,i10)') nxe_vis(i)
         read(8,'(25x,i10)') nys_vis(i)
         read(8,'(25x,i10)') nye_vis(i)
      end do
      close(8)

      dirname = adjustl(dirname)
      if(OverWrite .ne. 0 .and. OverWrite .ne. 1) stop "Odd OverWrite Value"
      if(Dstep < 1)   stop "Odd Dstep."
      if(step < 1)    step = Dstep
      do i=1,Nplane
         if(nxe_vis(i) < 1) nxe_vis(i) = ni(i)
         if(nye_vis(i) < 1) nye_vis(i) = nj(i)
      end do
   end subroutine read_control!}}}
subroutine center_to_grid(plane,arr)!{{{
      implicit none
      integer,intent(in)::plane
      double precision,dimension(0:nimax+1,0:njmax+1),intent(inout)::arr

      integer i,j

      do j=0,nj(plane)
         do i=0,ni(plane)
            arr(i,j)=1d0/4d0*(arr(i  ,j  )&
                             +arr(i  ,j+1)&
                             +arr(i+1,j  )&
                             +arr(i+1,j+1))
         end do
      end do
end subroutine center_to_grid!}}}
end program main
