module to_vis
   implicit none
   integer nxs_vis,nxe_vis,nys_vis,nye_vis
   character*50 dirname
   integer Dstep,step,OverWrite
   double precision tt
end module to_vis

program main
   use grbl_prmtr
   use to_vis
   implicit none
   double precision,dimension(ni,nj)::x,r
   double precision,dimension(ni,nj)::rho_mat,T_mat,u_mat,v_mat,p_mat
   double precision,dimension(ni,nj)::M_mat,R_gas,gmma,w_res,ei_mat,MW_mat
   double precision,dimension(dimq,ni,nj)::q
   double precision,dimension(ni,nj,dimq)::w_mat
   double precision rho,u,v,T
   integer nx_bin,ny_bin,step_bin
   integer i,j,k
   double precision,parameter::R_uni=8.314472d3
   character*50 filename,title
   integer*4,external::access

   !read control parameters
   call read_control

   !read geometry binary
   open(55,file=trim(dirname)//"geometry.bin",form="unformatted")
   read(55) nx_bin
   read(55) ny_bin
   read(55) ((x(i,j),r(i,j),i=1,ni),j=1,nj)
   close(55)

   !process time history
   do
      write(filename,'("result/result",i12.12)') step
      if(access(trim(dirname)//trim(filename)//".bin","r") .ne. 0) exit
      write(*       ,'("result/result",i12.12)') step
      if(OverWrite .eq. 0 .and. access(trim(filename)//".vtk","r") .eq. 0) then
         step=step+Dstep
         cycle
      end if

      !read q binary
      open(55,file=trim(dirname)//trim(filename)//".bin",form="unformatted")
      read(55) step_bin
      read(55) tt
      read(55) (((q(i,j,k),i=1,dimq),R_gas(j,k),gmma(j,k),p_mat(j,k),j=1,ni),k=1,nj)
      close(55)

      !Process data
      do j=1,nj
         do i=1,ni
            rho=0d0
            do k=1,nY
               rho=rho+q(k,i,j)
            end do
            u  =q(nY+1,i,j)/rho
            v  =q(nY+2,i,j)/rho
            T  =p_mat( i,j)/rho/R_gas(i,j)

            rho_mat(i,j)=rho
              T_mat(i,j)=T
              u_mat(i,j)=u
              v_mat(i,j)=v
              M_mat(i,j)=sqrt((u**2+v**2) / (gmma(i,j)*p_mat(i,j)/rho))
             ei_mat(i,j)=q(nY+3,i,j)/rho -0.5d0*(u**2+v**2)
             MW_mat(i,j)=R_uni/R_gas(i,j)

            w_res(i,j)=0d0
            do k=1,nY
               w_mat(i,j,k)=q(k,i,j)/rho
               if(k>3) then
                  w_res(i,j)=w_res(i,j)+w_mat(i,j,k)
               end if
            end do
         end do
      end do

      !output file open
      open(66,file=trim(filename)//".vtk")

      call vtk_header(x,r)

      !thermodynamical
      call vtk_scalar(rho_mat,"Density")
      call vtk_scalar(u_mat,"Velocity_x")
      call vtk_scalar(v_mat,"Velocity_y")
      call vtk_scalar(p_mat,"Pressure")
      call vtk_scalar(gmma,"SpecificHeatRatio")
      call vtk_vector(u_mat,v_mat,"Velocity")
      call vtk_scalar(T_mat,"Temperature")
      call vtk_scalar(MW_mat,"Molecular_Weight")
      call vtk_scalar(ei_mat,"InternalEnergy")

      call vtk_scalar(M_mat,"MachNumber")

      !other conservative quantities
      do k=1,min(3,nY)
         write(title,'(a,i3.3)') 'phi',k
         call vtk_scalar(w_mat(:,:,k),title)
      end do

      call vtk_scalar(w_res,"phi_res")

      !output file close
      close(66)

      step=step+Dstep
   end do
contains
   !begin of subroutines{{{
   subroutine vtk_header(x,r)
      implicit none
      double precision,dimension(ni,nj),intent(in)::x,r
      integer i,j

      !header
      write(66,'(a26)') '# vtk DataFile Version 2.0'
      write(66,'(a7)') '2D Data'
      write(66,'(a5)') 'ASCII'
      write(66,'(a23)') 'DATASET STRUCTURED_GRID'
      write(66,'(a11,i3,a1,i3,a1,i1)') 'DIMENSIONS ', nxe_vis-nxs_vis+1, ' ', nye_vis-nys_vis+1, ' ', 1


      write(66,'(a17)') 'FIELD FieldData 2'

      write(66,'(a15)') 'TIME 1 1 double'
      write(66,'(e14.6e3)')  tt
      write(66,'(a13)') 'CYCLE 1 1 int'
      write(66,'(i0)')  step

      !grid
      write(66,'(a7,i7,a6)') 'POINTS ', (nxe_vis-nxs_vis+1)*(nye_vis-nys_vis+1)*1, ' float'
      do j = nys_vis, nye_vis
        do i = nxs_vis, nxe_vis
          write(66,'(e14.6e3,2(1x,e14.6e3))') x(i,j), r(i,j), 0.d0
        end do
      end do
      write(66,*)
 
      write(66,'(a11,i7)') 'POINT_DATA ', (nxe_vis-nxs_vis+1)*(nye_vis-nys_vis+1)*1
      write(66,'(a22)') 'VECTORS Indices float'
      do j = nys_vis, nye_vis
        do i = nxs_vis, nxe_vis
          write(66,'(e14.6e3,2(1x,e14.6e3))') dble(i), dble(j), 0.d0
        end do
      end do
      write(66,*)
   end subroutine vtk_header
   subroutine vtk_scalar(arr,title)
      implicit none
      double precision,dimension(ni,nj),intent(in)::arr
      character(*),intent(in)::title

      integer i,j

      write(66,'(a)') 'SCALARS '//trim(title)//' float 1'
      write(66,'(a20)') 'LOOKUP_TABLE default'
      do j = nys_vis, nye_vis
        do i = nxs_vis, nxe_vis
          write(66,'(e14.6e3)') arr(i,j)
        end do
      end do
      write(66,*)
   end subroutine vtk_scalar
   subroutine vtk_vector(arr1,arr2,title)
      implicit none
      double precision,dimension(ni,nj),intent(in)::arr1,arr2
      character(*),intent(in)::title

      integer i,j

      write(66,'(a22)') 'VECTORS '//title//' float'
      do j = nys_vis, nye_vis
        do i = nxs_vis, nxe_vis
          write(66,'(e14.6e3,2(1x,e14.6e3))') arr1(i,j), arr2(i,j), 0.d0
        end do
      end do
      write(66,*)
   end subroutine vtk_vector
   subroutine read_control
      implicit none

      open(8,file="control.inp")

      read(8,'()')
      read(8,'(25x,i10)') Dstep
      read(8,'(25x,i10)') step
      read(8,'(25x,a)')   dirname
      read(8,'(25x,i10)') OverWrite
      read(8,'(25x,i10)') nxs_vis
      read(8,'(25x,i10)') nxe_vis
      read(8,'(25x,i10)') nys_vis
      read(8,'(25x,i10)') nye_vis
      close(8)

      dirname = adjustl(dirname)
      if(OverWrite .ne. 0 .and. OverWrite .ne. 1) stop "Odd OverWrite Value"
      if(Dstep < 1)   stop "Odd Dstep."
      if(step < 1)    step = Dstep
      if(nxe_vis < 1) nxe_vis = ni
      if(nye_vis < 1) nye_vis = nj
   end subroutine read_control
   !end of subroutines}}}
end program main
