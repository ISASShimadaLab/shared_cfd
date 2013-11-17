subroutine out_bin(step)!{{{
   use variable
   use mod_mpi
   implicit none
   integer,intent(in)::step
   character*50 filename
   double precision tmp((dimq+3)*(bwmax+2)*(bwmax+2))
   integer bwxt,bwyt
   integer nxst,nyst
   integer nxet,nyet
   integer myidt,cnt
   integer i,j,k,plane
   integer ii,jj

   call set_bc_q

   if(myid .eq. 0) then
      do plane = 1, Nplane
         do j=1,ngy(plane)
            do i=1,ngx(plane)
               myidt=proc_list(i,j,plane)
               if(myidt .eq. 0) cycle !cycle if myidt=0 itself

               cnt=(dimq+3)*(nxe_mat(i,plane)-nxs_mat(i,plane)+3)&
                           *(nye_mat(j,plane)-nys_mat(j,plane)+3)

               call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,&
                             myidt,step,MPI_COMM_WORLD, istatus, ierr)

               cnt=0
               do jj=nys_mat(j,plane)-1,nye_mat(j,plane)+1
                  do ii=nxs_mat(i,plane)-1,nxe_mat(i,plane)+1
                     cnt=cnt+dimq+3
                     q(1:dimq,ii,jj,plane)=tmp(cnt-dimq-2:cnt-3)
                     w(indxR, ii,jj,plane)=tmp(cnt-2)
                     w(indxg, ii,jj,plane)=tmp(cnt-1)
                     w(     4,ii,jj,plane)=tmp(cnt) 
                  end do
               end do
            end do
         end do
      end do
      
      write(filename,'("result/result",i12.12,".bin")') step
      open(45,file=filename,form="unformatted")
      open(55,file="restart.bin",form="unformatted")
      write(45) step
      write(55) step
      write(45) tt
      write(55) tt
      write(45) ((((q(i,j,k,plane),i=1,dimq),w(indxR,j,k,plane),w(indxg,j,k,plane),w(4,j,k,plane),j=0,ni(plane)+1),k=0,nj(plane)+1),plane=1,Nplane)
      write(55) ((((q(i,j,k,plane),i=1,dimq),w(indxR,j,k,plane),w(indxg,j,k,plane),w(4,j,k,plane),j=0,ni(plane)+1),k=0,nj(plane)+1),plane=1,Nplane)
      close(45)
      close(55)
   else
      do plane=nps,npe
         cnt=0
         do j=nys(plane)-1,nye(plane)+1
            do i=nxs(plane)-1,nxe(plane)+1
               cnt=cnt+dimq+3
               tmp(cnt-dimq-2:cnt-3)=q(1:dimq,i,j,plane)
               tmp(cnt-2)           =w(indxR, i,j,plane)
               tmp(cnt-1)           =w(indxg, i,j,plane)
               tmp(cnt)             =w(     4,i,j,plane)
            end do
         end do

         call MPI_Send(tmp,cnt,MPI_DOUBLE_PRECISION,0,step,MPI_COMM_WORLD,ierr)
      end do
   end if
end subroutine out_bin!}}}

subroutine restart_bin(step_res)!{{{
   use variable
   use mod_mpi
   implicit none
   integer,intent(out)::step_res

   double precision tmp((dimq+3)*(bwmax+2)*(bwmax+2))
   double precision xi
   integer bwxt,bwyt
   integer nxst,nyst
   integer nxet,nyet
   integer myidt,cnt
   integer i,j,k,plane
   integer ii,jj

   step_res = 0
   if(myid .eq. 0) then
      open(45,file="restart.bin",form="unformatted")
      read(45) step_res
      read(45) tt
      read(45) ((((q(i,j,k,plane),i=1,dimq),w(indxR,j,k,plane),w(indxg,j,k,plane),w(4,j,k,plane),j=0,ni(plane)+1),k=0,nj(plane)+1),plane=1,Nplane)
      close(45)

      do plane=1,Nplane
         do j=1,ngy(plane)
            do i=1,ngx(plane)
               myidt=proc_list(i,j,plane)
               if(myidt .eq. 0) cycle !cycle if myidt=0 itself

               cnt=0
               do jj=nys_mat(j,plane)-1,nye_mat(j,plane)+1
                  do ii=nxs_mat(i,plane)-1,nxe_mat(i,plane)+1
                     cnt=cnt+dimq+3
                     tmp(cnt-dimq-2:cnt-3)=q(1:dimq,ii,jj,plane)
                     tmp(cnt-2)           =w(indxR, ii,jj,plane)
                     tmp(cnt-1)           =w(indxg, ii,jj,plane)
                     tmp(cnt)             =w(     4,ii,jj,plane)
                  end do
               end do

               call MPI_Send(tmp,cnt,MPI_DOUBLE_PRECISION,myidt,0,MPI_COMM_WORLD,ierr)
            end do
         end do
      end do
   else
      do plane=nps,npe
         cnt=(dimq+3)*(bwx(plane)+2)*(bwy(plane)+2)
         call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD, istatus, ierr)

         cnt=0
         do jj=nys(plane)-1,nye(plane)+1
            do ii=nxs(plane)-1,nxe(plane)+1
               cnt=cnt+dimq+3
               q(1:dimq,ii,jj,plane)=tmp(cnt-dimq-2:cnt-3)
               w(indxR, ii,jj,plane)=tmp(cnt-2)
               w(indxg, ii,jj,plane)=tmp(cnt-1)
               w(     4,ii,jj,plane)=tmp(cnt)
            end do
         end do
      end do
   end if

   call MPI_Bcast(step_res, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end subroutine restart_bin!}}}

subroutine set_bc_q!{{{
   use variable
   use mod_mpi
   implicit none
   integer i,j,k,plane
   double precision tmp

   do plane = nps,npe
      !$omp parallel do private(i,k,tmp)
      do j=nys(plane),nye(plane)
         i=nxs(plane)-1
         tmp=w(1,i,j,plane)
         do k=1,nY
            q(k,i,j,plane)= tmp*w(   k+4,i,j,plane)
         end do
         q(nY+1,i,j,plane)= tmp*w(     2,i,j,plane)
         q(nY+2,i,j,plane)= tmp*w(     3,i,j,plane)
         q(nY+3,i,j,plane)= tmp*w(indxht,i,j,plane)-w(4,i,j,plane)
         !print *,myid,plane,i,j,w(1,i,j,plane)

         i=nxe(plane)+1
         tmp=w(1,i,j,plane)
         do k=1,nY
            q(k,i,j,plane)= tmp*w(   k+4,i,j,plane)
         end do
         q(nY+1,i,j,plane)= tmp*w(     2,i,j,plane)
         q(nY+2,i,j,plane)= tmp*w(     3,i,j,plane)
         q(nY+3,i,j,plane)= tmp*w(indxht,i,j,plane)-w(4,i,j,plane)
         !print *,myid,plane,i,j,w(1,i,j,plane)
      end do
      !$omp end parallel do

      !$omp parallel do private(j,k,tmp)
      do i=nxs(plane)-1,nxe(plane)+1
         j=nys(plane)-1
         tmp=w(1,i,j,plane)
         do k=1,nY
            q(k,i,j,plane)= tmp*w(   k+4,i,j,plane)
         end do
         q(nY+1,i,j,plane)= tmp*w(     2,i,j,plane)
         q(nY+2,i,j,plane)= tmp*w(     3,i,j,plane)
         q(nY+3,i,j,plane)= tmp*w(indxht,i,j,plane)-w(4,i,j,plane)
         !print *,myid,plane,i,j,w(1,i,j,plane)

         j=nye(plane)+1
         tmp=w(1,i,j,plane)
         do k=1,nY
            q(k,i,j,plane)= tmp*w(   k+4,i,j,plane)
         end do
         q(nY+1,i,j,plane)= tmp*w(     2,i,j,plane)
         q(nY+2,i,j,plane)= tmp*w(     3,i,j,plane)
         q(nY+3,i,j,plane)= tmp*w(indxht,i,j,plane)-w(4,i,j,plane)
         !print *,myid,plane,i,j,w(1,i,j,plane)
      end do
      !$omp end parallel do
   end do
end subroutine set_bc_q!}}}
