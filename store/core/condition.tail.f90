contains
subroutine section_exchange!{{{
   implicit none
   integer inc,i
   integer geo,width,order,p1,c1,s1,pm1,p2,c2,s2,pm2

   do i=1,num_cut_copro
      geo  = cut_copro( 1,i)
      width= cut_copro( 2,i)
      order= cut_copro( 3,i)
      p1   = cut_copro( 4,i)
      c1   = cut_copro( 5,i)
      s1   = cut_copro( 6,i)
      pm1  = cut_copro( 7,i)
      p2   = cut_copro( 8,i)
      c2   = cut_copro( 9,i)
      s2   = cut_copro(10,i)
      pm2  = cut_copro(11,i)
      select case(geo)
      case(1) !i-i
         do inc=0,width-1
            w(:,  c1+pm1,  s1+inc,      p1) = w(:,  c2,    s2+inc*order,p2)
            vhi(:,c1+pm1,  s1+inc,      p1) = vhi(:,c2,    s2+inc*order,p2)
            w(:,  c1+pm1*2,s1+inc,      p1) = w(:,  c2-pm2,s2+inc*order,p2)
            vhi(:,c1+pm1*2,s1+inc,      p1) = vhi(:,c2-pm2,s2+inc*order,p2)

            w(:,  c2+pm2,  s2+inc*order,p2) = w(:,  c1,    s1+inc,      p1)
            vhi(:,c2+pm2,  s2+inc*order,p2) = vhi(:,c1,    s1+inc,      p1)
            w(:,  c2+pm2*2,s2+inc*order,p2) = w(:,  c1-pm1,s1+inc,      p1)
            vhi(:,c2+pm2*2,s2+inc*order,p2) = vhi(:,c1-pm1,s1+inc,      p1)
         end do
      case(2) !i-j
         do inc=0,width-1
            w(:,  c1+pm1,  s1+inc,      p1) = w(:,  s2+inc*order,c2,    p2)
            vhi(:,c1+pm1,  s1+inc,      p1) = vhi(:,s2+inc*order,c2,    p2)
            w(:,  c1+pm1*2,s1+inc,      p1) = w(:,  s2+inc*order,c2-pm2,p2)
            vhi(:,c1+pm1*2,s1+inc,      p1) = vhi(:,s2+inc*order,c2-pm2,p2)

            w(:,  s2+inc*order,c2+pm2,  p2) = w(:,  c1,    s1+inc,      p1)
            vhi(:,s2+inc*order,c2+pm2,  p2) = vhi(:,c1,    s1+inc,      p1)
            w(:,  s2+inc*order,c2+pm2*2,p2) = w(:,  c1-pm1,s1+inc,      p1)
            vhi(:,s2+inc*order,c2+pm2*2,p2) = vhi(:,c1-pm1,s1+inc,      p1)
         end do
      case(3) !j-i
         do inc=0,width-1
            w(:,  s1+inc,      c1+pm1,  p1) = w(:,  c2,    s2+inc*order,p2)
            vhi(:,s1+inc,      c1+pm1,  p1) = vhi(:,c2,    s2+inc*order,p2)
            w(:,  s1+inc,      c1+pm1*2,p1) = w(:,  c2-pm2,s2+inc*order,p2)
            vhi(:,s1+inc,      c1+pm1*2,p1) = vhi(:,c2-pm2,s2+inc*order,p2)

            w(:,  c2+pm2,  s2+inc*order,p2) = w(:,  s1+inc,      c1,    p1)
            vhi(:,c2+pm2,  s2+inc*order,p2) = vhi(:,s1+inc,      c1,    p1)
            w(:,  c2+pm2*2,s2+inc*order,p2) = w(:,  s1+inc,      c1-pm1,p1)
            vhi(:,c2+pm2*2,s2+inc*order,p2) = vhi(:,s1+inc,      c1-pm1,p1)
         end do
      case(4) !j-j
         do inc=0,width-1
            w(:,  s1+inc,      c1+pm1,  p1) = w(:,  s2+inc*order,c2,    p2)
            vhi(:,s1+inc,      c1+pm1,  p1) = vhi(:,s2+inc*order,c2,    p2)
            w(:,  s1+inc,      c1+pm1*2,p1) = w(:,  s2+inc*order,c2-pm2,p2)
            vhi(:,s1+inc,      c1+pm1*2,p1) = vhi(:,s2+inc*order,c2-pm2,p2)
                                                                        
            w(:,  s2+inc*order,c2+pm2,  p2) = w(:,  s1+inc,      c1,    p1)
            vhi(:,s2+inc*order,c2+pm2,  p2) = vhi(:,s1+inc,      c1,    p1)
            w(:,  s2+inc*order,c2+pm2*2,p2) = w(:,  s1+inc,      c1-pm1,p1)
            vhi(:,s2+inc*order,c2+pm2*2,p2) = vhi(:,s1+inc,      c1-pm1,p1)
         end do
      end select
   end do
end subroutine section_exchange!}}}
subroutine MPI_COMMUNICATIONS_CUT_LINE!{{{
   implicit none
   double precision  tmp(DLength*2*bwmax)
   double precision tmps(DLength*2*bwmax,Nplane)
   integer cnt,cnts(Nplane),flag(Nplane),toProc(Nplane)
   integer inc,i
   integer cind,width,order,p,c,s,pm

   !send
   do i = 1,NumMPIComm
      cind      = MPIComm(1,i)
      width     = MPIComm(2,i)
      order     = MPIComm(3,i)
      p         = MPIComm(4,i)
      c         = MPIComm(5,i)
      s         = MPIComm(6,i)
      pm        = MPIComm(7,i)
      toProc(i) = MPIComm(8,i)

      cnts(i)=0
      if(cind .eq. 0) then
         do inc=0,width-1
            !print '(a7,4i3,9es9.1)',"send",myid,p,c,s+inc*order,w(:,  c,   s+inc*order,p)
            tmps(cnts(i)+1:cnts(i)+nY,  i) = vhi(:,c,   s+inc*order,p); cnts(i)=cnts(i)+nY
            tmps(cnts(i)+1:cnts(i)+dimw,i) = w(:,  c,   s+inc*order,p); cnts(i)=cnts(i)+dimw
            !print '(a7,4i3,9es9.1)',"send",myid,p,c-pm,s+inc*order,w(:,  c-pm,   s+inc*order,p)
            tmps(cnts(i)+1:cnts(i)+nY,  i) = vhi(:,c-pm,s+inc*order,p); cnts(i)=cnts(i)+nY
            tmps(cnts(i)+1:cnts(i)+dimw,i) = w(:,  c-pm,s+inc*order,p); cnts(i)=cnts(i)+dimw
         end do
      else
         do inc=0,width-1
            tmps(cnts(i)+1:cnts(i)+nY,  i) = vhi(:,s+inc*order,c,   p); cnts(i)=cnts(i)+nY
            tmps(cnts(i)+1:cnts(i)+dimw,i) = w(:,  s+inc*order,c,   p); cnts(i)=cnts(i)+dimw
            tmps(cnts(i)+1:cnts(i)+nY,  i) = vhi(:,s+inc*order,c-pm,p); cnts(i)=cnts(i)+nY
            tmps(cnts(i)+1:cnts(i)+dimw,i) = w(:,  s+inc*order,c-pm,p); cnts(i)=cnts(i)+dimw
         end do
      end if
      call MPI_Isend(tmps(:,i),cnts(i),MPI_DOUBLE_PRECISION,toProc(i),0,MPI_COMM_WORLD,ireq(i),ierrs(i))
   end do

   !receive
   do i = 1,NumMPIComm
      cind      = MPIComm(1,i)
      width     = MPIComm(2,i)
      order     = MPIComm(3,i)
      p         = MPIComm(4,i)
      c         = MPIComm(5,i)
      s         = MPIComm(6,i)
      pm        = MPIComm(7,i)

      cnt=DLength*2*width
      call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,toProc(i),0,MPI_COMM_WORLD, istatus, ierr)
      cnt=0
      if(cind .eq. 0) then
         do inc=0,width-1
            vhi(:,c+pm,  s+inc*order,p) = tmp(cnt+1:cnt+nY);   cnt=cnt+nY
            w(:,  c+pm,  s+inc*order,p) = tmp(cnt+1:cnt+dimw); cnt=cnt+dimw
            !print '(a7,4i3,9es9.1)',"receive",myid,p,c+pm,s+inc*order,w(:,  c+pm,  s+inc*order,p)
            vhi(:,c+pm*2,s+inc*order,p) = tmp(cnt+1:cnt+nY);   cnt=cnt+nY
            w(:,  c+pm*2,s+inc*order,p) = tmp(cnt+1:cnt+dimw); cnt=cnt+dimw
            !print '(a7,4i3,9es9.1)',"receive",myid,p,c+pm*2,s+inc*order,w(:,c+pm*2,s+inc*order,p)
         end do
      else
         do inc=0,width-1
            vhi(:,s+inc*order,c+pm,  p) = tmp(cnt+1:cnt+nY);   cnt=cnt+nY
            w(:,  s+inc*order,c+pm,  p) = tmp(cnt+1:cnt+dimw); cnt=cnt+dimw
            vhi(:,s+inc*order,c+pm*2,p) = tmp(cnt+1:cnt+nY);   cnt=cnt+nY
            w(:,  s+inc*order,c+pm*2,p) = tmp(cnt+1:cnt+dimw); cnt=cnt+dimw
         end do
      end if
   end do

   !wait
   do i = 1,NumMPIComm
      call MPI_Wait(ireq(i),istatus,ierrs(i))
   end do
end subroutine MPI_COMMUNICATIONS_CUT_LINE!}}}
subroutine MPI_COMMUNICATIONS_I_DIRECTION!{{{
   implicit none
   double precision  tmp(DLength*2*(bwmax+2))
   double precision tmps(DLength*2*(bwmax+2),Nplane)
   integer ii,jj
   integer cnt,cnts(Nplane),flag(Nplane)

   !grid west to east
   !send
   do plane = nps,npe
      if(gx(plane) .ne. ngx(plane)) then
         cnts(plane)=0
         do jj=nys(plane),nye(plane)
            do ii=nxe(plane)-1,nxe(plane)
               tmps(cnts(plane)+1:cnts(plane)+nY,  plane) =vhi(:,ii,jj,plane)
               cnts(plane)=cnts(plane)+nY
               tmps(cnts(plane)+1:cnts(plane)+dimw,plane) =w(  :,ii,jj,plane)
               cnts(plane)=cnts(plane)+dimw
            end do
         end do
         flag(plane)=plane
         call MPI_Isend(tmps(:,plane),cnts(plane),MPI_DOUBLE_PRECISION,ge(plane),flag(plane),&
                         MPI_COMM_WORLD,ireq(plane),ierrs(plane))
      end if
   end do

   !receive
   do plane = nps,npe
      if(gx(plane) .ne. 1) then
         cnt=DLength*2*bwy(plane)
         call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gw(plane),plane,&
                         MPI_COMM_WORLD, istatus, ierr)

         cnt=0
         do jj=nys(plane),nye(plane)
            do ii=nxs(plane)-2,nxs(plane)-1
               vhi(:,ii,jj,plane)=tmp(cnt+1:cnt+nY)
               cnt=cnt+nY
               w(  :,ii,jj,plane)=tmp(cnt+1:cnt+dimw)
               cnt=cnt+dimw
            end do
         end do
      end if
   end do

   do plane = nps,npe
      if(gx(plane) .ne. ngx(plane)) call MPI_Wait(ireq(plane),istatus,ierrs(plane))
   end do

   !grid east to west
   !send
   do plane = nps,npe
      if(gx(plane) .ne. 1) then
         cnts(plane)=0
         do jj=nys(plane),nye(plane)
            do ii=nxs(plane),nxs(plane)+1
               tmps(cnts(plane)+1:cnts(plane)+nY  ,plane) =vhi(:,ii,jj,plane)
               cnts(plane)=cnts(plane)+nY                                   
               tmps(cnts(plane)+1:cnts(plane)+dimw,plane) =w(  :,ii,jj,plane)
               cnts(plane)=cnts(plane)+dimw
            end do
         end do
         flag(plane)=plane
         call MPI_Isend(tmps(:,plane),cnts(plane),MPI_DOUBLE_PRECISION,gw(plane),flag(plane),&
                         MPI_COMM_WORLD,ireq(plane),ierrs(plane))
      end if
   end do

   !receive
   do plane = nps,npe
      if(gx(plane) .ne. ngx(plane)) then
         cnt=DLength*2*bwy(plane)
         call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,ge(plane),plane,&
                         MPI_COMM_WORLD, istatus, ierr)

         cnt=0
         do jj=nys(plane),nye(plane)
            do ii=nxe(plane)+1,nxe(plane)+2
               vhi(:,ii,jj,plane)=tmp(cnt+1:cnt+nY)
               cnt=cnt+nY                           
               w(  :,ii,jj,plane)=tmp(cnt+1:cnt+dimw)
               cnt=cnt+dimw
            end do
         end do
      end if
   end do

   do plane = nps,npe
      if(gx(plane) .ne. 1) call MPI_Wait(ireq(plane),istatus,ierrs(plane))
   end do
end subroutine MPI_COMMUNICATIONS_I_DIRECTION!}}}
subroutine MPI_COMMUNICATIONS_J_DIRECTION!{{{
   implicit none
   double precision  tmp(DLength*2*(bwmax+2))
   double precision tmps(DLength*2*(bwmax+2),Nplane)
   integer ii,jj
   integer cnt,cnts(Nplane),flag(Nplane)

   !grid south to north
   !send
   do plane = nps,npe
      if(gy(plane) .ne. ngy(plane)) then
         cnts(plane)=0
         do jj=nye(plane)-1,nye(plane)
            do ii=nxs(plane)-1,nxe(plane)+1
               tmps(cnts(plane)+1:cnts(plane)+nY,  plane) =vhi(:,ii,jj,plane)
               cnts(plane)=cnts(plane)+nY
               tmps(cnts(plane)+1:cnts(plane)+dimw,plane) =w(  :,ii,jj,plane)
               cnts(plane)=cnts(plane)+dimw
            end do
         end do
         flag(plane)=plane
         call MPI_Isend(tmps(:,plane),cnts(plane),MPI_DOUBLE_PRECISION,gn(plane),flag(plane),&
                         MPI_COMM_WORLD,ireq(plane),ierrs(plane))
      end if
   end do

   !receive
   do plane = nps,npe
      if(gy(plane) .ne. 1) then
         cnt=DLength*2*(bwx(plane)+2)
         call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gs(plane),plane,&
                         MPI_COMM_WORLD, istatus, ierr)

         cnt=0
         do jj=nys(plane)-2,nys(plane)-1
            do ii=nxs(plane)-1,nxe(plane)+1
               vhi(:,ii,jj,plane)=tmp(cnt+1:cnt+nY)
               cnt=cnt+nY
               w(  :,ii,jj,plane)=tmp(cnt+1:cnt+dimw)
               cnt=cnt+dimw
            end do
         end do
      end if
   end do

   do plane = nps,npe
      if(gy(plane) .ne. ngy(plane)) call MPI_Wait(ireq(plane),istatus,ierrs(plane))
   end do

   !grid north to south
   !send
   do plane = nps,npe
      if(gy(plane) .ne. 1) then
         cnts(plane)=0
         do jj=nys(plane),nys(plane)+1
            do ii=nxs(plane)-1,nxe(plane)+1
               tmps(cnts(plane)+1:cnts(plane)+nY,  plane) =vhi(:,ii,jj,plane)
               cnts(plane)=cnts(plane)+nY
               tmps(cnts(plane)+1:cnts(plane)+dimw,plane) =w(  :,ii,jj,plane)
               cnts(plane)=cnts(plane)+dimw
            end do
         end do
         flag(plane)=plane
         call MPI_Isend(tmps(:,plane),cnts(plane),MPI_DOUBLE_PRECISION,gs(plane),flag(plane),&
                         MPI_COMM_WORLD,ireq(plane),ierrs(plane))
      end if
   end do

   !receive
   do plane = nps,npe
      if(gy(plane) .ne.  ngy(plane)) then
         cnt=DLength*2*(bwx(plane)+2)
         call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,gn(plane),plane,&
                         MPI_COMM_WORLD, istatus, ierr)
      
         cnt=0
         do jj=nye(plane)+1,nye(plane)+2
            do ii=nxs(plane)-1,nxe(plane)+1
               vhi(:,ii,jj,plane)=tmp(cnt+1:cnt+nY)
               cnt=cnt+nY
               w(  :,ii,jj,plane)=tmp(cnt+1:cnt+dimw)
               cnt=cnt+dimw
            end do
         end do
      end if
   end do
   
   do plane = nps,npe
      if(gy(plane) .ne. 1) call MPI_Wait(ireq(plane),istatus,ierrs(plane))
   end do
end subroutine MPI_COMMUNICATIONS_J_DIRECTION!}}}
subroutine set_corners!{{{
   implicit none
   integer i,j,plane
   integer pi,pj

   do plane = nps,npe
      do pi = -1,1,2
         if(pi <0) then
            i = nxs(plane)
         else
            i = nxe(plane)
         end if
         do pj = -1,1,2
            if(pj <0) then
               j = nys(plane)
            else
               j = nye(plane)
            end if
            w(  :,i+pi,j+pj,plane) = &
                     1d0/3d0*(w(  :,i,j,plane)+w(  :,i,j+pj,plane)+w(  :,i+pi,j,plane))
            vhi(:,i+pi,j+pj,plane) = &
                     1d0/3d0*(vhi(:,i,j,plane)+vhi(:,i,j+pj,plane)+vhi(:,i+pi,j,plane))
         end do
      end do
   end do
end subroutine set_corners!}}}
end subroutine set_BC

