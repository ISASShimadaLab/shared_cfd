module conditions
   implicit none
   double precision,allocatable,dimension(:)::list_p
   double precision,allocatable,dimension(:)::list_T
   double precision,allocatable,dimension(:)::list_Yf
   integer nump
   integer numT
   integer numYf
end module conditions

subroutine read_conditions
   use conditions
   implicit none
   character*200 bufnext,line
   double precision lb,ub
   integer i


   open(22,file='control.inp')
   bufnext=""

   line=next_word()
   read(line,*) lb
   line=next_word()
   read(line,*) ub
   line=next_word()
   read(line,*) nump
   if(lb >ub) then
      stop "lower bound is larger than upper bound."
   else if(lb.eq. ub) then
      nump=1
      allocate(list_p(nump))
      list_p=lb
   else
      if(nump<=0) stop "Number of Sample point is lower than 1."
      allocate(list_p(nump))
      do i=1,nump
         list_p(i)=lb*(ub/lb)**(dble(i-1)/dble(nump-1))
      end do
   end if

   line=next_word()
   read(line,*) lb
   line=next_word()
   read(line,*) ub
   line=next_word()
   read(line,*) numT
   if(lb >ub) then
      stop "lower bound is larger than upper bound."
   else if(lb.eq. ub) then
      numT=1
      allocate(list_T(numT))
      list_T=lb
   else
      if(numT<=0) stop "Number of Sample point is lower than 1."
      allocate(list_T(numT))
      do i=1,numT
         list_T(i)=lb+(ub-lb)*dble(i-1)/dble(numT-1)
      end do
   end if

   line=next_word()
   read(line,*) lb
   line=next_word()
   read(line,*) ub
   line=next_word()
   read(line,*) numYf
   if(lb >ub) then
      stop "lower bound is larger than upper bound."
   else if(lb.eq. ub) then
      numYf=1
      allocate(list_Yf(numYf))
      list_Yf=lb
   else
      if(numYf<=0) stop "Number of Sample point is lower than 1."
      allocate(list_Yf(numYf))
      do i=1,numYf
         list_Yf(i)=lb+(ub-lb)*dble(i-1)/dble(numYf-1)
      end do
   end if

   close(22)
contains
character*200 function next_word()
   implicit none
   integer ind
   if(bufnext .eq. "") bufnext=next_line()
   ind = index(bufnext,' ')
   if(ind == 0) then
      next_word=bufnext
      bufnext=""
   else
      next_word=bufnext(:ind-1)
      bufnext=trim(adjustl(bufnext(ind:)))
   end if
end function next_word
character*200 function next_line()
   implicit none
   integer ind,ind0
   do while(bufnext .eq. "")
      read(22,'(a)') bufnext
      ind  =       index(bufnext,'#')
      ind0 = max(1,index(bufnext,':'))
      if(ind == 1) then
         bufnext=""
      else if(ind > 1) then
         bufnext=bufnext(ind0+1:ind-1)
      else
         bufnext=bufnext(ind0+1:)
      end if
      bufnext=trim(adjustl(bufnext))
   end do
   next_line=bufnext
   bufnext=""
end function next_line
end subroutine read_conditions

program driver
   use chem
   use chem_var
   use conditions
   implicit none
   character*100    buf
   double precision p,H,T, MWave,kappa,mu
   double precision Y(2)
   double precision,dimension(max_ns)::n,vhi,Yv,Yvmax
   integer         ,dimension(max_ns)::ind
   integer i,j,k,l,Ntic

   integer          itemp
   double precision dtemp

   !read files
   call init_pack_cea('hp')
   call read_conditions

   !start calculations
   T=300d0
   n=no+initial_eps
   Yvmax=0d0
   do k=1,numYf
      Y(1)=list_Yf(k)
      Y(2)=1d0-Y(1)
      do j=1,numT
         call calcH(list_T(j),Y,H)
         do i=1,numP
            p = list_p(i)
            call cea_hp(p,Y,H, T,n, MWave,kappa,mu,Yv,vhi)
            do l=1,ns
               if(Yv(l)>Yvmax(l)) Yvmax(l)=Yv(l)
            end do
         end do
      end do
   end do

   !!! sort
   ! initialization
   do j=1,ns
      ind(j)=j
   end do

   ! sort main
   do i=1,ns-1
      do j=2,ns-i+1
         if(Yvmax(j)>Yvmax(j-1)) then
            dtemp     =Yvmax(j)
            Yvmax(j)  =Yvmax(j-1)
            Yvmax(j-1)=dtemp
            itemp     =ind(  j)
            ind(  j)  =ind(  j-1)
            ind(  j-1)=itemp
         end if
      end do
   end do
   write(20,'(i3,es15.7,a21)')  (i,Yvmax(i)," # "//SYM_SPC(ind(i)),i=1,ns)
   print *,"done."
end program driver
