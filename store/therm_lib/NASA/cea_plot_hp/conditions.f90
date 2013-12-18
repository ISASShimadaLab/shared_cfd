module conditions
   use chem
   use chem_var
   implicit none
   double precision,allocatable,dimension(:)::list_p
   double precision,allocatable,dimension(:)::list_T
   double precision,allocatable,dimension(:)::list_Yf
   integer nump
   integer numT
   integer numYf
   integer numtype
contains
subroutine read_conditions!{{{
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

   line=next_word()
   if(     trim(line) .eq. "P") then
      numtype=1
   else if(trim(line) .eq. "T") then
      numtype=2
   else if(trim(line) .eq. "Yv") then
      numtype=3
   else
      stop "Odd variable for x axis."
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
end subroutine read_conditions!}}}
subroutine calcTeq(flag_cea)!{{{
   logical,intent(out)::flag_cea

   double precision p,H,T, MWave,kappa,mu
   double precision Y(2)
   double precision,dimension(max_ns)::n,vhi,Yv
   integer i,j,k

   !start calculations
   T=To
   n=no+initial_eps
   do k=1,numYf
      Y(1)=list_Yf(k)
      Y(2)=1d0-Y(1)
      do j=1,numT
         call calcH(list_T(j),Y,H)
         write(line,'("",2(i3.3,"."),"dat")') j,k
         open(55,file=trim(line))
         write(55,*) "#T=",list_T(j)," Yf=",list_Yf(k)

         write(55,'(a)') "# MF = mass fraction, oxid = oxidizer, prod = product, &
                   &br = before reaction, &
                   &MW = Molecular Weight, kappa = specific heat ratio, &
                   &vis = viscosity coefficient"
         write(55,'(a1,a14,100a15)') "#","pressure(Pa)","br fuel MF","br oxid MF",&
                        "Temperature(K)","Enthalpy(J/kg)","MW (g/mol)",&
                        "kappa","vis(Pa*s)"

         do i=1,numP
            p = list_p(i)
            call cea_hp(p,Y,H, T,n, MWave,kappa,mu,Yv,vhi,flag_cea)
            write(55,'(100es15.7)') p,Y,T,H,MWave,kappa,mu
            if(.not.flag_cea) return
         end do
         close(55)
      end do
   end do
end subroutine calcTeq!}}}
end module conditions
