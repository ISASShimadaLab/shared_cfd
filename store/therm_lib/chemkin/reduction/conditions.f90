module conditions
   implicit none
   double precision,dimension(:),      allocatable::list_T
   double precision,dimension(:),      allocatable::list_P
   double precision,dimension(:),      allocatable::list_Yf
   double precision,dimension(:,:,:,:),allocatable::list_vrho
   double precision,dimension(:,:,:,:),allocatable::Tref
   integer numT
   integer numP
   integer numYf
   double precision allowable_limit
   logical flag_debug
contains
subroutine read_conditions!{{{
   implicit none
   character*200 bufnext,line
   double precision lb,ub
   integer i

   flag_debug=.false.

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
   read(line,*) allowable_limit

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
subroutine init_vrho!{{{
   use chem
   use chem_var
   implicit none
   double precision Y(2),vw(max_ns),sn,p,rho
   integer i,j,k
   allocate(list_vrho(max_ns,numT,numP,numYf))
   allocate(Tref(          2,numT,numP,numYf))

   do k=1,numYf
      Y(1)=list_Yf(k)
      Y(2)=1d0-Y(1)
      vw=vwf*Y(1)+vwo*Y(2)
      sn=sum(vw/MWs,1)
      do j=1,numP
         p=list_P(j)
         do i=1,numT
            rho=p/(sn*Ru*1d-4*list_T(i))
            list_vrho(:,i,j,k)=rho*vw
         end do
      end do
   end do
end subroutine init_vrho!}}}
subroutine set_Tref!{{{
   implicit none
   double precision tign,Teq
   integer i,j,k,sm
   logical flag

   open(22,file="calculatable.out")
   open(23,file="not_calculatable.out")
   open(24,file="Tref.out")
   sm=0
   do k=1,numYf
      do j=1,numP
         do i=1,numT
            call reaction_one_steps(list_T(i),list_vrho(:,i,j,k),0,Tref(:,i,j,k),flag)
            if(flag) then
               print *,i,j,k,"OK"
               write(22,*) list_P(j),list_Yf(k),list_T(i)
               write(24,'(5es15.7)') list_P(j),list_Yf(k),list_T(i),Tref(:,i,j,k)
               sm=sm+1
            else
               print *,i,j,k,"NG"
               write(23,*) list_P(j),list_Yf(k),list_T(i)
               Tref(:,i,j,k)=-1d0
            end if
         end do
      end do
   end do
   close(22)
   close(23)
   close(24)

   print *,"(calculated points)/(all points)=",sm,"/",numYf*numP*numT
   if(sm .eq. 0) stop "Error:There is no calculatable point. Stop calculation."

   call check_model(flag)
   if(.not.flag) stop "Consistency check failed at set_Tref."
end subroutine set_Tref!}}}

subroutine check_model(flag)!{{{
   use chem
   implicit none
   logical         ,intent(out)::flag

   double precision T
   double precision vrho(ns)

   double precision Tini,tign,Teq
   double precision tt,to
   double precision n(ns+1)
   double precision RWORK(LRW)
   integer          IWORK(LIW)
   integer          istate
   double precision atol

   double precision Trefn(2)

   logical flag_vode

   integer i,j,k

   flag   = .false.

   do k=1,numYf
      do j=1,numP
         do i=1,numT
            Tini=list_T(i)
            vrho=list_vrho(:,i,j,k)
            tign=Tref(     1,i,j,k)
            Teq =Tref(     2,i,j,k)
            if(tign<0) cycle

            call reaction_one_steps(list_T(i),list_vrho(:,i,j,k),0,Trefn,flag_vode)
            !if(flag_debug) write(25,'(5es15.7)') list_P(j),list_Yf(k),list_T(i),Trefn(1),Trefn(2)
            if(.not.flag_vode) return 
            if(abs(Trefn(1)-tign)>tign*allowable_limit) return
            if(abs(Trefn(2)- Teq)> Teq*allowable_limit) return
         end do
      end do
   end do
   flag=.true.
end subroutine check_model!}}}
subroutine change_nr_via_ns!{{{
   use chem
   implicit none
   integer i,k,n
   logical flag
   double precision dtmp
   integer          itmp
   nr=max_nr
   i =max_nr
   do while(i>0)
      flag=.true.
      if(NumNu(1,i)<1)  flag = .false.
      do k=1,NumNu(1,i)
         if(IndNu(k,1,i)>ns) flag=.false.
      end do
      do k=1,NumNu(2,i)
         if(IndNu(k,2,i)>ns) flag=.false.
      end do

      if(flag .and. exist_M(i)) then
         !calc original NumM
         n=1
         do
            if(IndM(n,i) .eq. 0) then
               n=n-1
               exit
            end if
            if(n .eq. max_nr) exit
            n=n+1
         end do

         k=n
         do while(k>0)
            if(IndM(k,i)>ns) then
               itmp     = IndM(n,i)
               IndM(n,i)= IndM(k,i)
               IndM(k,i)= itmp

               dtmp     = Men(n,i)
               Men(n,i) = Men(k,i)
               Men(k,i) = dtmp
               n=n-1
            end if
            k=k-1
         end do
         NumM(i)=n
      end if

      if(.not.flag) then
         call swap_reaction(nr,i)
         nr=nr-1
      end if
      i=i-1
   end do
end subroutine change_nr_via_ns!}}}
subroutine remove_remaining_reactions!{{{
   use chem
   implicit none
   integer i
   logical flag
   do i=nr+1,max_nr
      SYM_RCT(  i)=""
      duplicate(i)=.false.
      exist_M(  i)=.false.
      Rstate( i,1)=0
      Rstate( i,2)=0
      IndNu(:,:,i)=0
      NumNu(:,  i)=0
      snu(      i)=0
      IndM(:,   i)=0
      NumM(     i)=0
      Men(:,    i)=0d0
      ABE(:,    i)=0d0
      cABE(:,   i)=0d0
      TROE(:,   i)=0d0
   end do
   call check_model(flag)
   if(.not.flag) stop "Something Odd occured at remove_remaining_reactions"
end subroutine remove_remaining_reactions!}}}

subroutine sort_n(flag)!{{{
   use chem
   implicit none
   logical         ,intent(in)::flag

   double precision Ymax(ns)

   integer          i,j,k
   double precision dtemp

   logical flag_vode

   call check_model(flag_vode)
   if(.not.flag_vode) stop "Something Odd occured at sort_n"

   Ymax=0d0
   do k=1,numYf
      do j=1,numP
         do i=1,numT
            if(Tref(1,i,j,k)<0d0) cycle
            call reaction_one_steps(list_T(i),list_vrho(:,i,j,k),1,Ymax,flag_vode)
            if(.not.flag_vode) stop "Calculation failed at sort_n."
         end do
      end do
   end do

   ! sort --- bubble sort
   do i=1,ns-1-num_reac
      do j=num_reac+2,ns-i+1
         if(Ymax(j)>Ymax(j-1)) then
            dtemp    =Ymax(j)
            Ymax(j)  =Ymax(j-1)
            Ymax(j-1)=dtemp

            call swap_species(j-1,j)
         end if
      end do
   end do

   if(flag) then
      !write out
      open(20,file="sort.dat")
      write(20,'(i3,es15.7,a21)')  (i,Ymax(i)," # "//SYM_SPC(i),i=1,ns)
      close(20)
   end if
end subroutine sort_n!}}}
subroutine reduction_species!{{{
   use chem
   implicit none
   integer ns_org,i
   logical flag

   call sort_n(.true.)

   ns_org=ns
   i=1
   do
      write(*,'(10x,a,i2)') "reduct_species_in_order. try number=",i
      call reduct_species_in_order
      if(ns .eq. ns_org) exit
      ns_org=ns
      call sort_n(.false.)
      i=i+1
   end do
   call reduct_every_species

   call check_model(flag)
   if(.not.flag) stop "Something Odd occured at reduction_species"
end subroutine reduction_species!}}}
subroutine reduct_species_in_order!{{{
   use chem
   implicit none
   logical flag
   integer nsl,nsu

   nsl=num_reac
   nsu=ns
   do while(nsl.ne.nsu)
      ns=(nsl+nsu)/2
      call change_nr_via_ns
      call check_model(flag)
      if(flag) then
         nsu=ns
      else
         nsl=ns
         if(nsu-nsl .eq. 1) exit
      end if
   end do
   ns=nsu
   call change_nr_via_ns
end subroutine reduct_species_in_order!}}}
subroutine reduct_every_species!{{{
   use chem
   implicit none
   logical flag
   integer ind
   ind = ns
   ns  = ns-1
   do while(ind .ne. num_reac)
      call swap_species(ind,ns+1)
      call change_nr_via_ns
      call check_model(flag)
      if(flag) then
         write(*,'(10x,a,i3)') "reduct_every_species. current ns=",ns
         ind=ns
         ns =ns-1
         cycle
      end if
      ind=ind-1
   end do
   ns=ns+1
   call change_nr_via_ns
   call out_cheminp
end subroutine reduct_every_species!}}}

subroutine sort_r(flag)!{{{
   use chem
   implicit none
   logical,intent(in)::flag

   double precision csp(nr)

   integer          i,j,k

   double precision dtemp

   logical flag_vode

   csp=0d0
   do k=1,numYf
      do j=1,numP
         do i=1,numT
            if(Tref(1,i,j,k)<0d0) cycle
            call reaction_one_steps(list_T(i),list_vrho(:,i,j,k),2,csp,flag_vode)
            if(.not.flag_vode) stop "Calculation failed at sort_r."
         end do
      end do
   end do

   ! sort --- bubble sort
   do i=1,nr-1
      do j=2,nr-i+1
         if(csp(j)>csp(j-1)) then
            dtemp   =csp(j)
            csp(j)  =csp(j-1)
            csp(j-1)=dtemp

            call swap_reaction(j-1,j)
         end if
      end do
   end do

   if(flag) then
      !write out
      open(20,file="sort_r.dat")
      write(20,'(i3,es15.7,a83)')  (i,csp(i)," # "//SYM_RCT(i),i=1,nr)
      close(20)
   end if
end subroutine sort_r!}}}
subroutine reduction_reactions!{{{
   use chem
   implicit none
   integer nr_org,i

   call sort_r(.true.)

   nr_org=nr
   i=1
   do
      write(*,'(10x,a,i2)') "reduct_reaction_in_order. try number=",i
      call reduct_reactions_in_order
      if(nr .eq. nr_org) exit
      nr_org=nr
      call sort_r(.false.)
      i=i+1
   end do
   call reduct_every_reactions
end subroutine reduction_reactions!}}}
subroutine reduct_reactions_in_order!{{{
   use chem
   implicit none
   logical flag
   integer nrl,nru

   nrl=1
   nru=nr
   do while(nrl.ne.nru)
      nr=(nrl+nru)/2
      call check_model(flag)
      if(flag) then
         nru=nr
      else
         nrl=nr
         if(nru-nrl .eq. 1) exit
      end if
   end do
   nr=nru
end subroutine reduct_reactions_in_order!}}}
subroutine reduct_every_reactions!{{{
   use chem
   implicit none
   logical flag
   integer ind

   ind = nr
   nr  = nr-1
   do while(ind .ne. 0)
      call swap_reaction(ind,nr+1)
      call check_model(flag)
      if(flag) then
         write(*,'(10x,a,i3)') "reduct_every_reaction. current nr=",nr
         ind=nr
         nr =nr-1
         cycle
      end if
      ind=ind-1
   end do
   nr=nr+1
end subroutine reduct_every_reactions!}}}
end module conditions

