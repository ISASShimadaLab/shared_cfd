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
   double precision,dimension(max_ns)::Yvmax

   double precision,allocatable,dimension(:,:,:)::list_Teq_org

   double precision limit
   integer          num_reac
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
   read(line,*) limit

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
subroutine count_reactants!{{{
   implicit none
   double precision n(max_ns)
   integer,dimension(max_ns)::list_reac
   integer i,itmp

   num_reac=0
   n=nfini+noini
   do i=1,ns
      if(n(i).ne.0d0) then
         num_reac=num_reac+1
         list_reac(num_reac)=i
      end if
   end do
   do i=1,num_reac
      call swap_species(i,list_reac(i))
   end do
end subroutine count_reactants!}}}
subroutine swap_species(ind1,ind2)!{{{
   implicit none
   integer,intent(in)::ind1,ind2

   integer indd1,indd2
   integer i,j

   if(ind1 .eq. ind2) return

   call swap_integer(num_sctn(ind1),num_sctn(ind2))
   call swap_real(MW(ind1),MW(ind2))

   !Trange
   do j=1,6
      do i=1,2
         call swap_real(Trange(i,j,ind1),Trange(i,j,ind2))
      end do
   end do

   !co
   do j=1,6
      do i=1,9
         call swap_real(co(i,j,ind1),co(i,j,ind2))
      end do
   end do

   !SYM_SPC
   call swap_char(SYM_SPC(ind1),SYM_SPC(ind2))

   !Ac
   do i=1,ne
      call swap_real(Ac(i,ind1),Ac(i,ind2))
   end do
   !nfo
   call swap_real(nf(ind1),nf(ind2))
   call swap_real(no(ind1),no(ind2))
   !nfoini
   call swap_real(nfini(ind1),nfini(ind2))
   call swap_real(noini(ind1),noini(ind2))
   !maskfo
   call swap_real(maskf(ind1),maskf(ind2))
   call swap_real(masko(ind1),masko(ind2))

   !tr2th
   indd1=0;indd2=0
   do i=1,nt
      if(tr2th(i) .eq. ind1) indd1=i
      if(tr2th(i) .eq. ind2) indd2=i
   end do
   if(indd1 .ne. 0) tr2th(indd1)=ind2
   if(indd2 .ne. 0) tr2th(indd2)=ind1
contains
subroutine swap_integer(a,b)
   integer a,b
   integer tmp
   tmp = a
   a   = b
   b   = tmp
end subroutine swap_integer
subroutine swap_real(a,b)
   double precision a,b
   double precision tmp
   tmp = a
   a   = b
   b   = tmp
end subroutine swap_real
subroutine swap_char(a,b)
   character*18 a,b
   character*18 tmp
   tmp = a
   a   = b
   b   = tmp
end subroutine swap_char
end subroutine swap_species!}}}
subroutine set_list_Teq_org!{{{
   logical flag_cea
   allocate(list_Teq_org(numP,numT,numYf))
   call calcTeq(list_Teq_org,flag_cea)
   if(.not.flag_cea) stop "calculation does not converge at original model."
end subroutine set_list_Teq_org!}}}
subroutine calcTeq(list_Teq,flag_cea)!{{{
   double precision,dimension(numP,numT,numYf),intent(out)::list_Teq
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
         do i=1,numP
            p = list_p(i)
            call cea_hp(p,Y,H, T,n, MWave,kappa,mu,Yv,vhi,flag_cea)
            list_Teq(i,j,k)=T
            if(.not.flag_cea) return
         end do
      end do
   end do
end subroutine calcTeq!}}}
subroutine compare(maxdiff,flag)!{{{
   implicit none
   double precision,intent(out)::maxdiff
   logical         ,intent(out)::flag

   double precision list_Teq(numP,numT,numYf)
   double precision diff
   integer i,j,k

   call calcTeq(list_Teq,flag)
   if(.not.flag) return

   maxdiff=0d0
   do k=1,numYf
      do j=1,numT
         do i=1,numP
            diff=abs(list_Teq(i,j,k)-list_Teq_org(i,j,k))&
                   /(list_Teq(i,j,k)+list_Teq_org(i,j,k))&
                   /2d0
            if(diff>maxdiff) maxdiff=diff
         end do
      end do
   end do
end subroutine compare!}}}
subroutine sort_species(flag)!{{{
   logical,intent(in)::flag
   double precision p,H,T, MWave,kappa,mu
   double precision Y(2)
   double precision,dimension(max_ns)::n,vhi,Yv
   integer i,j,k,l,Ntic
   integer,dimension(max_ns)::iind

   logical flag_cea
   integer          itemp
   double precision dtemp

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
            call cea_hp(p,Y,H, T,n, MWave,kappa,mu,Yv,vhi,flag_cea)
            do l=1,ns
               if(Yv(l)>Yvmax(l)) Yvmax(l)=Yv(l)
            end do
         end do
      end do
   end do

   ! sort --- bubble sort
   do i=1,ns-1-num_reac
      do j=num_reac+2,ns-i+1
         if(Yvmax(j)>Yvmax(j-1)) then
            dtemp     =Yvmax(j)
            Yvmax(j)  =Yvmax(j-1)
            Yvmax(j-1)=dtemp

            call swap_species(j-1,j)
         end if
      end do
   end do

   if(flag) then
      !write out
      open(20,file="sort.dat")
      write(20,'(i3,es15.7,a21)')  (i,Yvmax(i)," # "//SYM_SPC(i),i=1,ns)
      close(20)
   end if
end subroutine sort_species!}}}
subroutine reduction!{{{
   implicit none
   integer ns_org
   ns_org=ns
   do
      call reduct_in_order
      if(ns .eq. ns_org) exit
      ns_org=ns
      call sort_species(.false.)
   end do
   call reduct_every_species
end subroutine reduction!}}}
subroutine reduct_in_order!{{{
   implicit none
   logical flag
   double precision maxdiff
   integer nsl,nsu
   nsl=num_reac
   nsu=ns
   do while(nsl.ne.nsu)
      ns=(nsl+nsu)/2
      call compare(maxdiff,flag)
      if(.not. flag .or. maxdiff>limit) then
         nsl=ns
         if(nsu-nsl .eq. 1) exit
      else
         nsu=ns
      end if
   end do
   ns=nsu
end subroutine reduct_in_order!}}}
subroutine reduct_in_order_plot!{{{
   implicit none
   logical flag
   double precision maxdiff
   integer ns_org

   ns_org=ns
   do ns=ns_org-1,num_reac+1,-1
      call compare(maxdiff,flag)
      if(flag) write(20,*) ns,maxdiff
   end do
end subroutine reduct_in_order_plot!}}}
subroutine reduct_every_species!{{{
   implicit none
   logical flag
   double precision maxdiff
   integer ind

   ind = ns
   ns  = ns-1
   do while(ind .ne. num_reac)
      call swap_species(ind,ns+1)
      call compare(maxdiff,flag)
      if(flag .and. maxdiff<limit) then
         ind=ns
         ns =ns-1
         cycle
      end if
      ind=ind-1
   end do
   ns=ns+1
end subroutine reduct_every_species!}}}
subroutine out_chem!{{{
   implicit none
   integer i
   double precision maxdiff
   logical flag
   call compare(maxdiff,flag)
   print '(a,f5.2)',"final error(%) : ",maxdiff*1d2
   print '(a,i5)'  ,"final ns       : ",ns

   open(33,file="chem.inp.new")
   write(33,'(a)') "elements"
   do i=1,ne
      write(33,'(x,a)') trim(SYM_ELM(i))
   end do
   write(33,'(a)') "end"
   write(33,'(a)') "species"
   do i=1,ns
      write(33,'(x,a)') trim(SYM_SPC(i))
   end do
   write(33,'(a)') "end"
   write(33,'(a)') "thermo"
   write(33,'(a)') "end"
   write(33,'(a)') "reactions cal/mole  moles"
   write(33,'(a)') "end"
   close(33)
end subroutine out_chem!}}}
end module conditions
