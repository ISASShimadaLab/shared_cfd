subroutine set_IC
   use grbl_prmtr
   use prmtr
   use variable
   use mod_mpi
   use chem_var
   implicit none
   double precision rho,u,v,T,p,E,gamm,xi
   double precision q_temp(dimq)
   integer i,j

   T=300d0

   u  =0.603d0
   v  =0d0

   q_temp(1:nY)=vrhof(1:nY)
   q_temp(nY+1)=rhof*u
   q_temp(nY+2)=rhof*v
   q_temp(nY+3)=rhof*(Ef+0.5d0*(u**2+v**2))

   do i=nxs,ni/2
     do j=nys,nye
        q(:,    i,j)=q_temp
        w(4,    i,j)=1d5
        w(indxg,i,j)=kappaf
        w(indxR,i,j)=R_uni/MWf
     end do
   end do

   q_temp(1:nY)=vrhoo(1:nY)
   q_temp(nY+1)=rhoo*u
   q_temp(nY+2)=rhoo*v
   q_temp(nY+3)=rhoo*(Eo+0.5d0*(u**2+v**2))

   do i=ni/2+1,nxe
     do j=nys,nye
        q(:,    i,j)=q_temp
        w(4,    i,j)=1d5
        w(indxg,i,j)=kappao
        w(indxR,i,j)=R_uni/MWo
     end do
   end do
end subroutine set_IC

subroutine set_BC(step)
   use mod_mpi
   use grbl_prmtr
   use prmtr
   use variable
   use chem_var
   implicit none
   integer,intent(in)::step
   double precision u,v,E,MW,gamm,xi,T
   integer i,j
   double precision,parameter::StepDuration=5d3 !combustion

   integer,parameter::DLength=dimw+nY !for MPI Communication


   !boundary right and left
   if(gx .eq. 1) then
      !i=1/2
      !subsonic inflow
      do j=nys,nye
         w(     1,0,j)=w(4,1,j)*MWf/(R_uni*Tf)
         w(     2,0,j)=0.603d0
         w(     3,0,j)=0d0
         w(     4,0,j)=w(4,1,j)

         w(4+1:4+nY,0,j)=vwf(1:nY)

         w(indxg, 0,j)=kappaf
         w(indxht,0,j)=Ef+R_uni/MWf*Tf+0.5d0*(w(2,0,j)**2+w(3,0,j)**2)
         w(indxR, 0,j)=R_uni/MWf
         w(indxMu,0,j)=muf
         vhi(:,   0,j)=vhif

         w(:,  -1,j) = w(:,  0,j)
         vhi(:,-1,j) = vhi(:,0,j)
      end do
   end if

   if(gx .eq. ngx) then
      !i=ni+1/2
      !subsonic outflow
      do j=nys,nye
         w(:,  ni+1,j) = w(:,  ni,j)
         w(4,  ni+1,j) = 1d5
         vhi(:,ni+1,j) = vhi(:,ni,j)

         w(:,  ni+2,j) = w(:,  ni+1,j)
         vhi(:,ni+2,j) = vhi(:,ni+1,j)
      end do
   end if

   call MPI_COMMUNICATIONS_I_DIRECTION

   !boundary upper and lower
   if(gy .eq. 1) then
      !j=1/2 lower wall
      !slip wall
      do i=nxs-1,nxe+1
         w(:,  i, 0) =  w(:,  i,1)
         w(2,  i, 0) =  w(2,  i,1)
         w(3,  i, 0) = -w(3,  i,1)
         vhi(:,i, 0) =  vhi(:,i,1)

         w(:,  i,-1) =  w(:,  i,0)
         vhi(:,i,-1) =  vhi(:,i,0)
      end do
   end if

   if(gy .eq.  ngy) then
      !j=nj+1/2 upper wall
      !slip wall
      do i=nxs-1,nxe+1
         w(:,  i,nj+1) =  w(:,  i,nj)
         w(2,  i,nj+1) =  w(2,  i,nj)
         w(3,  i,nj+1) = -w(3,  i,nj)
         vhi(:,i,nj+1) =  vhi(:,i,nj)

         w(:,  i,nj+2) =  w(:,  i,nj+1)
         vhi(:,i,nj+2) =  vhi(:,i,nj+1)
      end do
   end if

   call MPI_COMMUNICATIONS_J_DIRECTION
contains
subroutine MPI_COMMUNICATIONS_I_DIRECTION!{{{
   implicit none
   double precision  tmp(DLength*2*(bwmax+2))
   double precision tmps(DLength*2*(bwmax+2))
   integer ii,jj
   integer cnt,cnts

   !grid west to east
   !send
   if(gx .ne. ngx) then
      cnts=0
      do jj=nys,nye
         do ii=nxe-1,nxe
            tmps(cnts+1:cnts+nY) =vhi(:,ii,jj)
            cnts=cnts+nY
            tmps(cnts+1:cnts+dimw)    =w(  :,ii,jj)
            cnts=cnts+dimw
         end do
      end do
      call MPI_Isend(tmps,cnts,MPI_DOUBLE_PRECISION,ge,0,MPI_COMM_WORLD,ireq,ierrs)
   end if

   !receive
   if(gx .ne. 1) then
      cnt=DLength*2*bwy
      call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gw,0,MPI_COMM_WORLD, istatus, ierr)

      cnt=0
      do jj=nys,nye
         do ii=nxs-2,nxs-1
            vhi(:,ii,jj)=tmp(cnt+1:cnt+nY)
            cnt=cnt+nY
            w(  :,ii,jj)=tmp(cnt+1:cnt+dimw)
            cnt=cnt+dimw
         end do
      end do
   end if

   if(gx .ne. ngx) call MPI_Wait(ireq,istatus,ierrs)

   !grid east to west
   !send
   if(gx .ne. 1) then
      cnts=0
      do jj=nys,nye
         do ii=nxs,nxs+1
            tmps(cnts+1:cnts+nY) =vhi(:,ii,jj)
            cnts=cnts+nY
            tmps(cnts+1:cnts+dimw)    =w(  :,ii,jj)
            cnts=cnts+dimw
         end do
      end do
      call MPI_Isend(tmps,cnts,MPI_DOUBLE_PRECISION,gw,0,MPI_COMM_WORLD,ireq,ierrs)
   end if

   !receive
   if(gx .ne. ngx) then
      cnt=DLength*2*bwy
      call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,ge,0,MPI_COMM_WORLD, istatus, ierr)

      cnt=0
      do jj=nys,nye
         do ii=nxe+1,nxe+2
            vhi(:,ii,jj)=tmp(cnt+1:cnt+nY)
            cnt=cnt+nY
            w(  :,ii,jj)=tmp(cnt+1:cnt+dimw)
            cnt=cnt+dimw
         end do
      end do
   end if

   if(gx .ne. 1) call MPI_Wait(ireq,istatus,ierrs)
end subroutine MPI_COMMUNICATIONS_I_DIRECTION!}}}
subroutine MPI_COMMUNICATIONS_J_DIRECTION!{{{
   implicit none
   double precision  tmp(DLength*2*(bwmax+2))
   double precision tmps(DLength*2*(bwmax+2))
   integer ii,jj
   integer cnt,cnts

   !grid south to north
   !send
   if(gy .ne. ngy) then
      cnts=0
      do jj=nye-1,nye
         do ii=nxs-1,nxe+1
            tmps(cnts+1:cnts+nY) =vhi(:,ii,jj)
            cnts=cnts+nY
            tmps(cnts+1:cnts+dimw)    =w(  :,ii,jj)
            cnts=cnts+dimw
         end do
      end do
      call MPI_Isend(tmps,cnts,MPI_DOUBLE_PRECISION,gn,0,MPI_COMM_WORLD,ireq,ierrs)
   end if

   !receive
   if(gy .ne. 1) then
      cnt=DLength*2*(bwx+2)
      call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gs,0,MPI_COMM_WORLD, istatus, ierr)

      cnt=0
      do jj=nys-2,nys-1
         do ii=nxs-1,nxe+1
            vhi(:,ii,jj)=tmp(cnt+1:cnt+nY)
            cnt=cnt+nY
            w(  :,ii,jj)=tmp(cnt+1:cnt+dimw)
            cnt=cnt+dimw
         end do
      end do
   end if

   if(gy .ne. ngy) call MPI_Wait(ireq,istatus,ierrs)

   !grid north to south
   !send
   if(gy .ne. 1) then
      cnts=0
      do jj=nys,nys+1
         do ii=nxs-1,nxe+1
            tmps(cnts+1:cnts+nY) =vhi(:,ii,jj)
            cnts=cnts+nY
            tmps(cnts+1:cnts+dimw)    =w(  :,ii,jj)
            cnts=cnts+dimw
         end do
      end do
      call MPI_Isend(tmps,cnts,MPI_DOUBLE_PRECISION,gs,0,MPI_COMM_WORLD,ireq,ierrs)
   end if

   !receive
   if(gy .ne.  ngy) then
      cnt=DLength*2*(bwx+2)
      call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gn,0,MPI_COMM_WORLD, istatus, ierr)
   
      cnt=0
      do jj=nye+1,nye+2
         do ii=nxs-1,nxe+1
            vhi(:,ii,jj)=tmp(cnt+1:cnt+nY)
            cnt=cnt+nY
            w(  :,ii,jj)=tmp(cnt+1:cnt+dimw)
            cnt=cnt+dimw
         end do
      end do
   end if
   
   if(gy .ne. 1) call MPI_Wait(ireq,istatus,ierrs)
end subroutine MPI_COMMUNICATIONS_J_DIRECTION!}}}
end subroutine set_BC

