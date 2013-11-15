subroutine out_bin(step)!{{{
   use variable
   use mod_mpi
   implicit none
   integer,intent(in)::step
   character*50 filename
   double precision tmp((dimq+3)*bwmax*bwmax)
   integer bwxt,bwyt
   integer nxst,nyst
   integer nxet,nyet
   integer myidt,cnt
   integer i,j,k
   integer ii,jj

   if(myid .eq. 0) then
      do j=1,ngy
         do i=1,ngx
            myidt=j+ngy*(i-1)-1
            if(myidt .eq. 0) cycle !cycle if myidt=0 itself

            cnt=(dimq+3)*(nxe_mat(i)-nxs_mat(i)+1)*(nye_mat(j)-nys_mat(j)+1)

            call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,&
                          myidt,step,MPI_COMM_WORLD, istatus, ierr)

            cnt=0
            do jj=nys_mat(j),nye_mat(j)
               do ii=nxs_mat(i),nxe_mat(i)
                  cnt=cnt+dimq+3
                  q(1:dimq,ii,jj)=tmp(cnt-dimq-2:cnt-3)
                  w(indxR,   ii,jj)=tmp(cnt-2)
                  w(indxg, ii,jj)=tmp(cnt-1)
                  w(     4,ii,jj)=tmp(cnt) 
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
      write(45) (((q(i,j,k),i=1,dimq),w(indxR,j,k),w(indxg,j,k),w(4,j,k),j=1,ni),k=1,nj)
      write(55) (((q(i,j,k),i=1,dimq),w(indxR,j,k),w(indxg,j,k),w(4,j,k),j=1,ni),k=1,nj)
      close(45)
      close(55)
   else
      cnt=0
      do j=nys,nye
         do i=nxs,nxe
            cnt=cnt+dimq+3
            tmp(cnt-dimq-2:cnt-3)=q(1:dimq,i,j)
            tmp(cnt-2)           =   w(indxR,i,j)
            tmp(cnt-1)           =w(indxg, i,j)
            tmp(cnt)             =w(     4,i,j)
         end do
      end do

      call MPI_Send(tmp,cnt,MPI_DOUBLE_PRECISION,0,step,MPI_COMM_WORLD,ierr)
   end if
end subroutine out_bin!}}}

subroutine restart_bin(step_res)!{{{
   use variable
   use mod_mpi
   implicit none
   integer,intent(out)::step_res

   double precision tmp((dimq+3)*bwmax*bwmax)
   double precision xi
   integer bwxt,bwyt
   integer nxst,nyst
   integer nxet,nyet
   integer myidt,cnt
   integer i,j,k
   integer ii,jj

   step_res = 0
   if(myid .eq. 0) then
      open(45,file="restart.bin",form="unformatted")
      read(45) step_res
      read(45) tt
      read(45) (((q(i,j,k),i=1,dimq),w(indxR,j,k),w(indxg,j,k),w(4,j,k),j=1,ni),k=1,nj)
      close(45)

      do j=1,ngy
         do i=1,ngx
            myidt=j+ngy*(i-1)-1
            if(myidt .eq. 0) cycle !cycle if myidt=0 itself

            cnt=0
            do jj=nys_mat(j),nye_mat(j)
               do ii=nxs_mat(i),nxe_mat(i)
                  cnt=cnt+dimq+3
                  tmp(cnt-dimq-2:cnt-3)=q(1:dimq,ii,jj)
                  tmp(cnt-2)           =   w(indxR,ii,jj)
                  tmp(cnt-1)           =w(indxg, ii,jj)
                  tmp(cnt)             =w(     4,ii,jj)
               end do
            end do

            call MPI_Send(tmp,cnt,MPI_DOUBLE_PRECISION,myidt,0,MPI_COMM_WORLD,ierr)
         end do
      end do
   else
      cnt=(dimq+3)*bwx*bwy
      call MPI_Recv(tmp,cnt,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD, istatus, ierr)

      cnt=0
      do jj=nys,nye
         do ii=nxs,nxe
            cnt=cnt+dimq+3
            q(1:dimq,ii,jj)=tmp(cnt-dimq-2:cnt-3)
               w(indxR,ii,jj)=tmp(cnt-2)
            w(indxg, ii,jj)=tmp(cnt-1)
            w(     4,ii,jj)=tmp(cnt)
         end do
      end do
   end if

   call MPI_Bcast(step_res, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end subroutine restart_bin!}}}
