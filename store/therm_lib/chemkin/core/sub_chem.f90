module func_therm
contains
subroutine check_section_number(T,sec)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in) ::T
   integer,         intent(out)::sec(ns)

   integer i
   
   do i=1,ns
      if(T<=Tthre(3,i)) then
         sec(i)=1
      else
         sec(i)=2
      end if
   end do
end subroutine check_section_number!}}}
subroutine set_coeff_less(T,logT,chrt,ccpr)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in)::T
   double precision,intent(in)::logT
   double precision,intent(out)::chrt(6)
   double precision,intent(out)::ccpr(6)

   double precision T2,T3,T4,Ti1

   Ti1=1d0/T
   T2 =T  *T
   T3 =T2 *T
   T4 =T3 *T

   chrt( 1)=1d0
   chrt( 2)=T /2d0
   chrt( 3)=T2/3d0
   chrt( 4)=T3/4d0
   chrt( 5)=T4/5d0
   chrt( 6)=Ti1

   ccpr( 1)=1d0
   ccpr( 2)=T
   ccpr( 3)=T2
   ccpr( 4)=T3
   ccpr( 5)=T4
   ccpr( 6)=0d0
end subroutine set_coeff_less!}}}
subroutine set_coeff(T,logT,cmurt,chrt,ccpr)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in)::T
   double precision,intent(in)::logT
   double precision,intent(out)::cmurt(7)
   double precision,intent(out)::chrt( 7)
   double precision,intent(out)::ccpr( 7)

   double precision T2,T3,T4,Ti1

   Ti1=1d0/T
   T2 =T  *T
   T3 =T2 *T
   T4 =T3 *T
   
   cmurt(1)= (1d0-logT)
   cmurt(2)=-T /2d0
   cmurt(3)=-T2/6d0
   cmurt(4)=-T3/12d0
   cmurt(5)=-T4/20d0
   cmurt(6)= Ti1
   cmurt(7)=-1d0

   chrt( 1)=1d0
   chrt( 2)=T /2d0
   chrt( 3)=T2/3d0
   chrt( 4)=T3/4d0
   chrt( 5)=T4/5d0
   chrt( 6)=Ti1
   chrt( 7)=0d0

   ccpr( 1)=1d0
   ccpr( 2)=T
   ccpr( 3)=T2
   ccpr( 4)=T3
   ccpr( 5)=T4
   ccpr( 6)=0d0
   ccpr( 7)=0d0
end subroutine set_coeff!}}}
subroutine set_coeff_more(T,logT,cmurt,chrt,ccpr,cdmurt,cdcpr)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in)::T
   double precision,intent(in)::logT
   double precision,intent(out)::cmurt( 7)
   double precision,intent(out)::chrt(  7)
   double precision,intent(out)::ccpr(  7)
   double precision,intent(out)::cdmurt(7)
   double precision,intent(out)::cdcpr( 7)

   double precision T2,T3,T4,Ti1,Ti2

   Ti1=1d0/T
   Ti2=Ti1*Ti1
   T2 =T  *T
   T3 =T2 *T
   T4 =T3 *T
   
   cmurt( 1)= (1d0-logT)
   cmurt( 2)=-T /2d0
   cmurt( 3)=-T2/6d0
   cmurt( 4)=-T3/12d0
   cmurt( 5)=-T4/20d0
   cmurt( 6)= Ti1
   cmurt( 7)=-1d0

   chrt(  1)=1d0
   chrt(  2)=T /2d0
   chrt(  3)=T2/3d0
   chrt(  4)=T3/4d0
   chrt(  5)=T4/5d0
   chrt(  6)=Ti1
   chrt(  7)=0d0

   ccpr(  1)=1d0
   ccpr(  2)=T
   ccpr(  3)=T2
   ccpr(  4)=T3
   ccpr(  5)=T4
   ccpr(  6)=0d0
   ccpr(  7)=0d0

   cdmurt(1)=-Ti1
   cdmurt(2)=-1d0/2d0
   cdmurt(3)=-T  /3d0
   cdmurt(4)=-T2 /4d0
   cdmurt(5)=-T3 /5d0
   cdmurt(6)=-Ti2
   cdmurt(7)= 0d0

   cdcpr( 1)=0d0
   cdcpr( 2)=1d0
   cdcpr( 3)=2d0*T
   cdcpr( 4)=3d0*T2
   cdcpr( 5)=4d0*T3
   cdcpr( 6)=0d0
   cdcpr( 7)=0d0
end subroutine set_coeff_more!}}}
subroutine check_section_number_trans(T,sec)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in):: T
   integer,         intent(out)::sec(nt)

   integer i,l

   do i=1,nt
      sec(i)=-1
      if (Trange_trans(1,1,i) > T) sec(i) = 1

      do l=1,num_sctn_trans(i)
         if(Trange_trans(2,l,i)>T) then
            sec(i)=l
            exit
         end if
      end do

      if(sec(i) .eq. -1) sec(i)=num_sctn_trans(i)
   end do
end subroutine check_section_number_trans!}}}
subroutine set_coeff_mu(T,logT,cmu)!{{{
   implicit none
   double precision,intent(in)::T
   double precision,intent(in)::logT
   double precision,intent(out)::cmu(4)

   double precision invT
 
   invT=1d0/T 
   cmu(1)=logT
   cmu(2)=invT
   cmu(3)=invT*invT
   cmu(4)=1d0
end subroutine set_coeff_mu!}}}
end module func_therm

! initialization
subroutine read_cheminp!{{{
   use mod_mpi
   use const_chem
   use chem
   implicit none
   double precision   MWe(max_ne)

   character*100    linebuf,str
   character*50     line_reactants,line_products
   character*24     cha
   character        sscha
   character*2      scha(5)
   integer          inte(5)
   double precision T_def(3)
   double precision T_buf(3)

   integer ind,Nfound
   integer i,j,k,l

   integer,parameter::LMW =86
   integer,parameter::LCHM=87

   if(myid .eq. 0) then
      open(LCHM,file="chem.inp",status="old")

      !!!!!!!!!!!!!!!!!!!!!! element read !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call next_line(LCHM,linebuf)
      if(trim(linebuf) .ne. "elements") stop "elements declaration not found."
      ne=0
      do
         call next_line(LCHM,linebuf)
         call chomp(linebuf)
         if(linebuf .eq. "end" .or. linebuf .eq. "END")  exit
         do
            ne=ne+1
            if(ne > max_ne) stop "chem.inp element number exceeded assigned one."
            call read_head_char(linebuf,cha)
            SYM_ELM(ne)=cha
            if(linebuf .eq. "") exit
         end do
      end do

      !read MW.log
      Nfound=0
      open(LMW,file="MW.inp",status="old")
      do
         call next_line(LMW,linebuf)
         if(linebuf .eq. "") exit

         call read_head_char(linebuf,cha)
         do i=1,ne
            if(trim(cha) .eq. trim(SYM_ELM(i))) then
               Nfound=Nfound+1
               call read_tail_number(linebuf,MWe(i))
               exit
            end if
         end do
      end do
      close(LMW)
      if(ne .ne. Nfound) stop "Error. There is Element which does not exit in MW.inp."
      !!!!!!!!!!!!!!!!!!!!!! END of element read !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!! species read !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call next_line(LCHM,linebuf)
      if(trim(linebuf) .ne. "species")  stop "species declaration not found."
      ns=0
      do
         call next_line(LCHM,linebuf)
         call chomp(linebuf)
         if(linebuf .eq. "end" .or. linebuf .eq. "END")  exit
         do
            ns=ns+1
            if(ns > max_ns) stop "chem.inp species number exceeded assigned one."
            call read_head_char(linebuf,cha)
            SYM_SPC(ns)=cha
            if(linebuf .eq. "") exit
         end do
      end do
      !!!!!!!!!!!!!!!!!!!!!! END of species read !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!! therm read !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call next_line(LCHM,linebuf)
      if(trim(linebuf) .ne. "thermo all") stop "'thermo all' declaration not found."

      Nfound=0
      !read default common threshold temperatures
      call next_line(LCHM,linebuf)
      read(linebuf,'(3F10.0)') T_def(1),T_def(2),T_def(3)
      do
         call next_line(LCHM,linebuf)
         if(linebuf .eq. "end" .or. linebuf .eq. "END")  exit
         read(linebuf,'(a24,4(a2,i3),a1,2f10.0,f8.0,(a2,i3))') &
            cha,(scha(i),inte(i),i=1,4),sscha,(T_buf(i),i=1,3),&
            scha(5),inte(5)
         ind=index(cha," ")
         cha=cha(:ind-1)
         call search_SYM(SYM_SPC,ns,cha,k)

         if(k.gt.0) then
            if(sscha .ne. "G") stop "Error! Other than Gas phase is contained."
            Nfound=Nfound+1
            do i=1,5
               if(len_trim(scha(i)) .eq. 0) cycle
               call search_SYM(SYM_ELM,ne,scha(i),j)
               ES(j,k)=inte(i)
            end do

            MWs(k)=dot_product(ES(:ne,k),MWe(:ne))

            !Temperature range
            if(T_buf(1) .ne. 0d0) then
               Tthre(:,k)=T_buf(:)
            else
               Tthre(:,k)=T_def(:)
            end if

            !Data
            call next_line(LCHM,linebuf)
            read(linebuf,'(5E15.1)') (coeff(i,2,k),i=1,5)
            call next_line(LCHM,linebuf)
            read(linebuf,'(5E15.1)') (coeff(i,2,k),i=6,7),(coeff(i,1,k),i=1,3)
            call next_line(LCHM,linebuf)
            read(linebuf,'(4E15.1)') (coeff(i,1,k),i=4,7)
         else
            call next_line(LCHM,linebuf)
            call next_line(LCHM,linebuf)
            call next_line(LCHM,linebuf)
         end if
      end do
      if(Nfound .ne. ns) stop "There is species which does not have thermal data."
      !!!!!!!!!!!!!!!!!!!!!! END of therm read !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!! reaction read !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call next_line(LCHM,linebuf)
      if(trim(linebuf) .ne. "reactions      cal/mole  moles") stop "'reactions' declaration not found."
      nr=0
      exist_M=.false.
      Rstate =0
      do
         call next_line(LCHM,linebuf)
         call chomp(linebuf)
         if(linebuf .eq. "end" .or. linebuf .eq. "END")  exit

         if(linebuf .eq. "duplicate" .or. linebuf .eq. "DUPLICATE") then
            duplicate(nr)=.true.
            cycle
         end if

         ind=index(linebuf,"/")
         if(ind .eq. 0) then ! reaction declaration line
            nr=nr+1
            if(nr > max_nr) stop "chem.inp reaction number exceeded assigned one."

            !read default ABE
            duplicate(nr)=.false. ! initialization
            call read_tail_number_arr(linebuf,3,ABE(:,nr))

            SYM_RCT(nr)=linebuf
            ind=index(linebuf,"<=>")
            if(ind .ne. 0) then
               line_reactants=linebuf(:ind-1)
               line_products =linebuf(ind+3:)
            else
               ind=index(linebuf,"=")
               line_reactants=linebuf(:ind-1)
               line_products =linebuf(ind+1:)
            end if

            if(line_products(1:1) .eq. ">") then
               Rstate(nr,1)=1
               line_products=line_products(2:len_trim(line_products))
            end if
           
            !about M
            ind=index(line_reactants,"(+M)")
            if(ind>0) then
               if(ind+3 .ne. len_trim(line_reactants)) stop "Error. (+M) must be end of each side."
               line_reactants=line_reactants(:len_trim(line_reactants)-4)
               line_products =line_products (:len_trim(line_products )-4)
               exist_M(nr)=.true.
            else
               ind=index(line_reactants,"+M")
               if(ind>0) then
                  if(ind+1 .ne. len_trim(line_reactants)) stop "Error. +M must be end of each side."
                  line_reactants=line_reactants(:len_trim(line_reactants)-2)
                  line_products =line_products (:len_trim(line_products )-2)
                  exist_M(nr)=.true.
                  Rstate(nr,2)=4
               end if
            end if

            !reactants
            i=1
            do
               call search_SYM(SYM_SPC,ns,line_reactants,j)
               IndNu(i,1,nr)=j
               if(line_reactants .eq. "") exit
               i=i+1
            end do
            NumNu(1,nr)=i

            !products
            i=1
            do
               call search_SYM(SYM_SPC,ns,line_products,j)
               IndNu(i,2,nr)=j
               if(line_products .eq. "") exit
               i=i+1
            end do
            NumNu(2,nr)=i
         else
            !auxiliary
            cha =trim(linebuf(:ind-1))
            linebuf=linebuf(ind:)
            if     (cha .eq. "REV"  .or. cha .eq. "rev")  then
                 if(Rstate(nr,1) .ne. 0) stop "Error at rev"
                 if(Rstate(nr,2) .ne. 0 .and. &
                    Rstate(nr,2) .ne. 4) stop "Error at rev"
                 Rstate(nr,2)=Rstate(nr,2)+1

                 call fetch_slash(linebuf,str)
                 call read_tail_number_arr(str,3,cABE(:,nr))
            else if(cha .eq. "LOW"  .or. cha .eq. "low")  then
                 if(Rstate(nr,1) .ne. 0) stop "Error at low"
                 if(Rstate(nr,2) .ne. 0) stop "Error at low"
                 Rstate(nr,2)=2

                 call fetch_slash(linebuf,str)
                 call read_tail_number_arr(str,3,cABE(:,nr))
            else if(cha .eq. "TROE" .or. cha .eq. "troe") then
                 if(Rstate(nr,1) .ne. 0) stop "Error at TROE"
                 if(Rstate(nr,2) .ne. 2) stop "Error at TROE"
                 Rstate(nr,2)=3

                 call fetch_slash(linebuf,str)
                 call read_tail_number_arr(str,4,TROE(:,nr))
            else !M enhanced
                 if(.not. exist_M(nr)) stop "Error at M enhanced"

                 Men(:,nr)=0d0
                 j=0
                 do
                    j=j+1
                    call search_SYM(SYM_SPC,ns,cha,IndM(j,nr))
                    call fetch_slash(linebuf,str)

                    call read_tail_number(str,Men(j,nr))
                    Men(j,nr)=Men(j,nr)-1d0
                    if(linebuf .eq. "") exit
                    ind =index(linebuf,"/")
                    cha    =linebuf(:ind-1)
                    linebuf=linebuf(ind:)
                 end do
                 NumM(nr)=j
            end if
         end if
      end do
      !!!!!!!!!!!!!!!!!!!!!! END of reaction read !!!!!!!!!!!!!!!!!!!!!!!!!!!!

95    close(LCHM)
     
      !!!!!!!!!!! post process !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,nr
         snu(i)=NumNu(2,i)-NumNu(1,i)
      end do
      do j=1,ns
         invMW(j)=1d0/MWs(j)
      end do
      nr=max(nr,1)
      !!!!!!!!!!! post process !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end if

   call MPI_Bcast(     ne,      1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(     ns,      1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(     nr,      1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  Tthre,   3*ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  coeff, 7*2*ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(    MWs,     ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  invMW,     ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast( Rstate,   2*nr,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  IndNu, 3*2*nr,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  NumNu,   2*nr,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(    snu,     nr,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(exist_M,     nr,          MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(   IndM,  ns*nr,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(   NumM,     nr,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(    Men,  ns*nr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(    ABE,   3*nr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(   cABE,   3*nr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(   TROE,   4*nr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
contains
subroutine read_tail_number_arr(line,num,arr)!{{{
   implicit none
   character(*),intent(inout)::line
   integer,intent(in)::num
   double precision,intent(out)::arr(num)

   integer ind,i

   do i=num,1,-1
      ind=index(trim(line)," ",.true.)
      read(line(ind+1:),*) arr(i)
      line=trim(line(:ind))
   end do
   if(num .eq. 3) then
      arr(1)=log(arr(1))
      arr(3)=arr(3)*invRc
   else if(num .eq. 4) then
      arr(2)=1d0/arr(2)
      arr(3)=1d0/arr(3)
   end if
end subroutine read_tail_number_arr!}}}
subroutine next_line(NUM,linebuf)!{{{
   implicit none
   integer,intent(in)::NUM
   character(100),intent(out)::linebuf

   integer ind

   do
      read(NUM,'(a100)',end=99) linebuf

      ind=index(linebuf,"!")
      if(ind .gt. 1) then
         linebuf=linebuf(:ind-1)
         linebuf=trim(linebuf)
      else if(ind .eq. 1) then
         cycle
      else
         linebuf=trim(linebuf)
      end if

      if(len_trim(linebuf) .ne. 0) exit
   end do
   return
99 linebuf=""
end subroutine next_line!}}}
subroutine chomp(line)!{{{
   character(*),intent(inout)::line
   integer ind

   line=trim(adjustl(line))
end subroutine chomp!}}}
subroutine read_tail_number(line,num)!{{{
   character(*),intent(inout)::line
   double precision,intent(out)::num

   integer ind

   ind=index(trim(line)," ",.true.)
   read(line(ind+1:),*) num
   line=trim(line(:ind))
end subroutine read_tail_number!}}}
subroutine read_head_char(line,cha)!{{{
   character(*),intent(inout)::line
   character(*),intent(out)::cha

   integer ind

   ind=index(trim(line)," ")
   if(ind .ne. 0) then
      cha=line(:ind-1)
      line=trim(adjustl(line(ind:)))
   else
      cha=line
      line=""
   end if
end subroutine read_head_char!}}}
subroutine search_SYM(SYM,NUM_SYM,str,i)!{{{ 
   implicit none
   character(*),dimension(*),intent(in)::SYM
   integer,intent(in)                  ::NUM_SYM
   character(*),intent(inout)          ::str
   integer,intent(out)                 ::i
   integer ind
   character*18 buf

   ind = index(str,"+")
   if(ind .eq. 0) then
      buf=trim(str)
      str=""
   else
      buf=trim(str(:ind-1))
      str=trim(str(ind+1:))
   end if

   i=1
   do while(i<=NUM_SYM)
      if(trim(SYM(i)) .eq. buf) exit
      i=i+1
   end do

   if(i>NUM_SYM) i=-1
end subroutine search_SYM!}}}
subroutine fetch_slash(line,cha)!{{{
   character(*),intent(inout)::line
   character(*),intent(out)::cha

   integer ind1,ind2

   ind1=index(line,         "/")
   ind2=index(line(ind1+1:),"/")
   ind2=ind2+ind1

   cha=line(ind1+1:ind2-1)
   line=trim(adjustl(line(ind2+1:)))
end subroutine fetch_slash!}}}
end subroutine read_cheminp!}}}
subroutine out_cheminp!{{{
   use mod_mpi
   use const_chem
   use chem
   implicit none
   character*100 line
   character*50  word1,word2
   character*80  NAME_RCT
   integer i,j,k

   integer ind_dup(nr),Ndup
   integer itmp

   integer,parameter::LCHM=87

   if(myid .eq. 0) then
      !check duplicate
      Ndup=0
      do k=1,nr
         if(duplicate(k)) then
            Ndup=Ndup+1
            ind_dup(Ndup)=k
         end if
      end do

      dup:do while(Ndup>0)
         k=Ndup
         NAME_RCT=SYM_RCT(ind_dup(k))
         do k=1,Ndup-1
            if(trim(NAME_RCT) .eq. trim(SYM_RCT(ind_dup(k)))) then
               Ndup=Ndup-1
               itmp          = ind_dup(k)
               ind_dup(k)    = ind_dup(Ndup)
               ind_dup(Ndup) = itmp
               Ndup=Ndup-1
               cycle dup
            end if
         end do
         duplicate(ind_dup(Ndup))=.false.
         Ndup=Ndup-1
      end do dup



      open(LCHM,file="chem.inp.new",status="replace")
      !credit
      write(LCHM,'(a)') "!"
      write(LCHM,'(a)') "! generated by SHIMADA CODE."
      write(LCHM,'(a)') "! written by Shota Yamanaka. 2013."
      write(LCHM,'(a)') "!"

      !elements
      write(LCHM,'(a)') "elements"
      line=" "
      do i=1,ne
         call append_line(SYM_ELM(i))
      end do
      if(len_trim(line)>0) write(LCHM,'(a)') trim(line)
      write(LCHM,'(a)') "end"

      !species
      write(LCHM,'(a)') "species"
      line=" "
      do i=1,ns
         call append_line(SYM_SPC(i))
      end do
      if(len_trim(line)>0) write(LCHM,'(a)') trim(line)
      write(LCHM,'(a)') "end"

      !thermo
      write(LCHM,'(a)') "thermo all"
      write(LCHM,'(a)') "   300.000  1000.000  5000.000"
      do j=1,ns
         write(line,'(a)') SYM_SPC(j)
         write(LCHM,'(a24)',advance='no') line
         line=""
         do i=1,ne
            if(ES(i,j) .ne. 0d0) then
               write(word1,'(a2,i3)') SYM_ELM(i),int(ES(i,j))
               line=trim(line)//word1(1:5)
            end if
         end do
         write(LCHM,'(a20)',advance='no') line
         write(LCHM,'("G",2f10.2,f8.2,i07)') Tthre(:,j),1
         write(LCHM,'(5E15.8,i05)')  (coeff(i,2,j),i=1,5),2
         write(LCHM,'(5E15.8,i05)')  (coeff(i,2,j),i=6,7),(coeff(i,1,j),i=1,3),3
         write(LCHM,'(4E15.8,i020)') (coeff(i,1,j),i=4,7),4
      end do
      write(LCHM,'(a)') "end"

      !reaction
      write(LCHM,'(a)') "reactions      cal/mole  moles"
      do i=1,nr
         if(len_trim(SYM_RCT(i))<25) then
            write(line,'(a)') SYM_RCT(i)
            write(LCHM,'(x,a25,2x,es13.3e3,f8.3,f10.1)') line,&
                          exp(ABE(1,i)),ABE(2,i),ABE(3,i)/invRc
         else
            write(LCHM,'(x,a,2x,es13.3e3,f8.3,f10.1)') trim(SYM_RCT(i)),&
                          exp(ABE(1,i)),ABE(2,i),ABE(3,i)/invRc
         end if

         if(duplicate(i)) write(LCHM,'(23x,"DUPLICATE")')

         itmp=Rstate(i,2)
         if(itmp .eq. 1 .or. itmp .eq. 5) then
            write(LCHM,'(23x,"REV /",es13.3e3,f8.3,f10.1," /")') &
                          exp(cABE(1,i)),cABE(2,i),cABE(3,i)/invRc
         else if(itmp .eq. 2 .or. itmp .eq. 3) then
            write(LCHM,'(23x,"LOW /",es13.3e3,f8.3,f10.1," /")') &
                          exp(cABE(1,i)),cABE(2,i),cABE(3,i)/invRc
            if(itmp .eq. 3) then
               line="TROE /"
               call dble2word1(TROE(1,i))
               line=trim(line)//" "//trim(word1)

               call dble2word1(1d0/TROE(2,i))
               line=trim(line)//" "//trim(word1)

               call dble2word1(1d0/TROE(3,i))
               line=trim(line)//" "//trim(word1)

               call dble2word1(TROE(4,i))
               line=trim(line)//" "//trim(word1)//"/"

               if(len_trim(line)<60) then
                  write(LCHM,'(a61)')  trim(line)
               else
                  write(LCHM,'(a)') trim(line)
               end if
            end if
         end if

         if(NumM(i)>0) then
            line=""
            do k=1,NumM(i)
               call dble2word1(Men(k,i)+1d0)
               write(word2,'(a,"/",a,"/")') trim(SYM_SPC(IndM(k,i))),trim(word1)
               line=trim(line)//" "//trim(word2)
            end do
            if(len_trim(line)<60) then
               write(LCHM,'(a61)')  trim(line)
            else
               write(LCHM,'(a)') trim(line)
            end if
         end if
      end do
      write(LCHM,'(a)') "end"
      close(LCHM)
   end if
contains
subroutine append_line(word)
   implicit none
   character(*) word
   line=trim(line)//" "//trim(word)
   if(len_trim(line)>40) then
      write(LCHM,'(a)') trim(line)
      line=" "
   end if
end subroutine append_line
subroutine dble2word1(val)
   implicit none
   double precision val
   character*50 pro2
   if((1d-1 < abs(val) .and. abs(val)<1d8) .or. val .eq. 0d0) then
      write(word1,'(f15.5)') val
      word1=trim(adjustl(word1))
      do while(word1(len_trim(word1):len_trim(word1)) .eq. "0")
         word1=word1(:len_trim(word1)-1)
      end do
      if(word1(len_trim(word1):len_trim(word1)) .eq. ".") &
         word1=word1(:len_trim(word1)-1)
   else
      write(word1,'(es15.6e3)') val
      pro2=word1(11:)
      if(pro2(3:3).eq.'0') pro2=pro2(1:2)//pro2(4:)

      word1=word1(:10)
      word1=trim(word1)
      do while(word1(len_trim(word1):len_trim(word1)) .eq. "0")
         word1=word1(:len_trim(word1)-1)
      end do
      if(word1(len_trim(word1):len_trim(word1)) .eq. ".") &
         word1=word1(:len_trim(word1)-1)
      word1=trim(word1)//pro2
   end if
   word1=trim(adjustl(word1))
end subroutine dble2word1
end subroutine out_cheminp!}}}
subroutine set_trans_data!{{{
   use mod_mpi
   use const_chem
   use chem
   implicit none
   character*18     sname,comment
   character*300    comment_long
   character  v, c
   integer   iv,ic
   integer   itr2th
   integer i,j


   if(myid .eq. 0) then
      open(10,file='trans.inp',status='old')

      !for first line comment
      read(10,'(A30)',END=999) comment_long

      nt=0
      do
            !fetch line
            read(10,'(a)',END=999) comment_long

            !for comments or END line
            if(comment_long(1:1) .eq. '!') cycle
            if(comment_long(1:3) .eq. 'END' .or. comment_long(1:3) .eq. 'end') exit

            !process first line
            read(comment_long,'(A15,x,A15,3x,A1,I1,A1,I1)') sname, comment,v,iv,c,ic

            !search species from thermo data
            itr2th=search_species(sname)

            !!! for binary interaction statements
            if(len(trim(comment)) .ne. 0 .or. itr2th .eq. -1) then
               do i=1,iv+ic
                  read(10,'(a)') comment_long
               end do
               cycle
            end if
            !!! end for binary interaction statements

            !!! the statements below are for binary interaction
            !increment nt and save sname as species_name_trans
            nt=nt+1
            species_name_trans(nt)=sname
            num_sctn_trans(nt)    =iv
            tr2th(nt)             =itr2th

            !for V
            do i=1,iv
               read(10,'(1x,a1,2(f7.1,2x),4e15.8)') v,(Trange_trans(j,i,nt),j=1,2),(trans(j,i,nt),j=1,4)
               if(v .ne. "V") stop "ERROR : ODD FORMAT AT trans.inp"
            end do

            !for C --- which are not used now.
            do i=1,ic
               read(10,'(a)') comment_long
               if(comment_long(2:2) .ne. "C") stop "ERROR : ODD FORMAT AT trans.inp"
            end do
      end do
999   continue

      !print *,"the number of species of trans.inp:",nt

      close(10)
   end if

   call MPI_Bcast(            nt,      1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(num_sctn_trans,     ns,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(         trans, 4*3*ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  Trange_trans, 2*3*ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(         tr2th,     ns,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
contains
   integer function search_species(str)
      use chem
      character(*),intent(in)::str
      integer i
   
      do i=1,ns
         if(trim(SYM_SPC(i)) .eq. trim(str)) then
            search_species=i
            return
         end if
      end do
      search_species=-1
   end function search_species
   integer function search_species_trans(str)
      use chem
      character(*),intent(in)::str
      integer i
   
      do i=1,nt
         if(trim(species_name_trans(i)) .eq. trim(str)) then
            search_species_trans=i
            return
         end if
      end do
      search_species_trans=-1
   end function search_species_trans
end subroutine set_trans_data!}}}
subroutine read_fo_composition!{{{
   use chem
   use mod_mpi
   use chem_var
   implicit none
   character*100 buf,buf2
   double precision amount
   integer ind

   double precision no(max_ns),nf(max_ns),vrho(max_ns)
   double precision DHit(max_ns)
   integer,dimension(max_ns)::list_reac
   integer i

   if(myid .eq. 0) then
      no=0d0
      nf=0d0

      open(8,file="control_chem.inp")
      read(8,'()')
      read(8,'(25x,es15.7)') po
      read(8,'(25x,es15.7)') To
      read(8,'()')
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit
         buf=adjustl(trim(buf)) !delete front and back spaces
         ind=index(buf,' ')
         if(ind .eq. 0) stop "Bad format at control_chem.inp while reading Oxidizer Compositions"
         read(buf(ind+1:),*) amount
         no(search_species(buf(1:ind-1))) = amount
      end do
      read(8,'(25x,es15.7)') pf
      read(8,'(25x,es15.7)') Tf
      read(8,'()')
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit
         buf=adjustl(trim(buf)) !delete front and back spaces
         ind=index(buf,' ')
         if(ind .eq. 0) stop "Bad format at control_chem.inp while reading Fuel Compositions"
         read(buf(ind+1:),*) amount
         nf(search_species(buf(1:ind-1))) = amount
      end do

      read(8,'(25x,es15.7)') rtol_user
      read(8,'(25x,es15.7)') atol_user
      close(8)
   end if

   !print *,"ns_tocalc = ",ns_tocalc

   !!! MPI COMMUNICATIONS
   call MPI_Bcast(po,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(To,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(vrhoo,   ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   call MPI_Bcast(pf,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(Tf,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(vrhof,   ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   call MPI_Bcast(ns_tocalc,1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   ns_tocalc=ns! temporal one

   of=sum(no*MWs,1)/sum(nf*MWs,1)
   vrhoo=no*MWs
   vrhof=nf*MWs

   !!! calculate other properties for oxidizer and fuel
   call rho_rel2abs(po,To, vrhoo, rhoo,Eo)
   call calc_T(vrhoo,rhoo,Eo, To, kappao,MWo,DHit,vhio,muo)
   vwo(1:nY) = vrhoo(1:nY)/rhoo

   call rho_rel2abs(pf,Tf, vrhof, rhof,Ef)
   call calc_T(vrhof,rhof,Ef, Tf, kappaf,MWf,DHit,vhif,muf)
   vwf(1:nY) = vrhof(1:nY)/rhof

   !move reactants to top of the list
   num_reac=0
   vrho=vrhof+vrhoo
   do i=1,ns
      if(vrho(i).ne.0d0) then
         num_reac=num_reac+1
         list_reac(num_reac)=i
      end if
   end do
   do i=1,num_reac
      call swap_species(i,list_reac(i))
   end do
   if(isColdFlow) then
      ns_tocalc=num_reac
   else
      ns_tocalc=ns
   end if
   call set_static_qw(pf,rhof,Tf,Ef,kappaf,muf,vrhof,qf,wf)
   call set_static_qw(po,rhoo,To,Eo,kappao,muo,vrhoo,qo,wo)
contains
   integer function search_species(str)
      use chem
      character(*),intent(in)::str
      integer i
   
      do i=1,ns
         if(trim(SYM_SPC(i)) .eq. trim(str)) then
            search_species=i
            return
         end if
      end do
      search_species=-1
   end function search_species
   subroutine set_static_qw(p,rho,T,E,kappa,mu,vrho,q_local,w_local)
      implicit none
      double precision,intent(in)::p
      double precision,intent(in)::rho
      double precision,intent(in)::T
      double precision,intent(in)::E
      double precision,intent(in)::kappa
      double precision,intent(in)::mu
      double precision,intent(in)::vrho(max_ns)
      double precision,intent(out)::q_local(dimq)
      double precision,intent(out)::w_local(dimw)
      integer i
      do i=1,ns
         q_local(i)   = vrho(i)
         w_local(i+4) = vrho(i)/rho
      end do
      q_local(nY+1)=0d0
      q_local(nY+2)=0d0
      q_local(nY+3)=rho*E

      w_local(1)=rho
      w_local(2)=0d0
      w_local(3)=0d0
      w_local(4)=p
      w_local(indxg )=kappa
      w_local(indxht)=E+p/rho
      w_local(indxR )=p/(rho*T)
      w_local(indxMu)=mu
   end subroutine set_static_qw
end subroutine read_fo_composition!}}}
subroutine init_therm!{{{
   call read_cheminp
   call set_trans_data
   call read_fo_composition
end subroutine init_therm!}}}
subroutine calc_vrho(p,T,deg, rho,vrho,E)!{{{
   use const_chem
   use chem
   use chem_var
   use func_therm
   implicit none
   double precision,intent(in)::p
   double precision,intent(in)::T
   double precision,intent(in)::deg
   double precision,intent(out)::rho
   double precision,intent(out)::vrho(ns)
   double precision,intent(out)::E

   vrho = rhof * deg + rhoo*(1d0-deg)

   call rho_rel2abs(p,T, vrho, rho,E)
end subroutine calc_vrho!}}}
subroutine rho_rel2abs(p,T, vrho, rho,E)!{{{
   use const_chem
   use chem
   use func_therm
   implicit none
   double precision,intent(in)::p
   double precision,intent(in)::T
   double precision,intent(inout)::vrho(ns)
   double precision,intent(out)::rho
   double precision,intent(out)::E

   double precision n(ns)
   double precision,dimension(7)::cmurt,chrt,ccpr
 
   double precision dh,tmp
   integer i,j
   integer sec(ns)

   !vrho to n
   do j=1,ns
      n(j)=vrho(j)*invMW(j)
   end do

   !calc n
   tmp=sum(n,1)*Ru*T
   tmp=p*1d1/tmp
   n  =n*tmp

   !calc rho
   rho=dot_product(n,MWs(1:ns))
   !convert unit from g/cm^3 to kg/m^3
   rho=rho*1d3

   !calc vrho
   do j=1,ns
      vrho(j)=n(j)*MWs(j)*1d3
   end do

   call check_section_number(T,sec)
   call set_coeff(T,log(T),cmurt,chrt,ccpr)
   E = 0d0
   do j=1,ns
      dh= 0d0
      do i=1,7
         dh=dh+coeff(i,sec(j),j)*chrt(i)
      end do
      E=E+(dh-1d0)*n(j)
   end do
   E=E*Ru*T/rho/1d1
end subroutine rho_rel2abs!}}}

! swapping
subroutine swap_species(ind1,ind2)!{{{
   use chem
   use chem_var
   implicit none
   integer,intent(in)::ind1,ind2

   integer i,j,k

   if(ind1 .eq. ind2) return

   call swap_char(SYM_SPC(ind1),SYM_SPC(ind2))
   do i=1,max_ne
      call swap_real(ES(i,ind1),ES(i,ind2))
   end do
   do i=1,3
      call swap_real(Tthre(i,ind1),Tthre(i,ind2))
   end do
   do j=1,2
      do i=1,7
         call swap_real(coeff(i,j,ind1),coeff(i,j,ind2))
      end do
   end do
   call swap_real(  MWs(ind1),  MWs(ind2))
   call swap_real(invMW(ind1),invMW(ind2))
   call swap_real(vrhof(ind1),vrhof(ind2))
   call swap_real(vrhoo(ind1),vrhoo(ind2))
   call swap_real(  vwf(ind1),  vwf(ind2))
   call swap_real(  vwo(ind1),  vwo(ind2))
   call swap_real( vhif(ind1), vhif(ind2))
   call swap_real( vhio(ind1), vhio(ind2))

   do k=1,max_nr
      do j=1,2
         do i=1,NumNu(j,k)
            call swap_integer_var(ind1,ind2,IndNu(i,j,k))
         end do
      end do

      do i=1,max_ns
         call swap_integer_var(ind1,ind2,IndM(i,k))
      end do
   end do

   do i=1,max_ns
      call swap_integer_var(ind1,ind2,tr2th(i))
   end do
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
subroutine swap_integer_var(ind1,ind2,var)
   integer ind1,ind2
   integer var

   if(     var .eq. ind1) then
      var=ind2
   else if(var .eq. ind2) then
      var=ind1
   end if
end subroutine swap_integer_var
end subroutine swap_species!}}}
subroutine swap_reaction(ind1,ind2)!{{{
   use chem
   use chem_var
   implicit none
   integer,intent(in)::ind1,ind2

   integer i,j,k

   if(ind1 .eq. ind2) return

   call swap_char(SYM_RCT(ind1),SYM_RCT(ind2))
   do j=1,2
      do i=1,3
         call swap_integer(IndNu(i,j,ind1),IndNu(i,j,ind2))
      end do
   end do
   do j=1,max_ns
      call swap_integer(IndM(j,ind1),IndM(j,ind2))
      call swap_real(    Men(j,ind1), Men(j,ind2))
   end do
   call swap_integer(Rstate(ind1,1),Rstate(ind2,1))
   call swap_integer(Rstate(ind1,2),Rstate(ind2,2))
   call swap_integer(NumNu(1,ind1),NumNu(1,ind2))
   call swap_integer(NumNu(2,ind1),NumNu(2,ind2))
   call swap_integer(  snu(  ind1),  snu(  ind2))
   call swap_integer( NumM(  ind1), NumM(  ind2))

   call swap_logical(duplicate(ind1),duplicate(ind2))

   call swap_logical(exist_M(ind1),exist_M(ind2))

   do i=1,3
      call swap_real( ABE(i,ind1), ABE(i,ind2))
      call swap_real(cABE(i,ind1),cABE(i,ind2))
      call swap_real(TROE(i,ind1),TROE(i,ind2))
   end do
   call swap_real(TROE(4,ind1),TROE(4,ind2))
contains
subroutine swap_logical(a,b)
   logical a,b
   logical tmp
   tmp = a
   a   = b
   b   = tmp
end subroutine swap_logical
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
   character(80) a,b
   character(80) tmp
   tmp = a
   a   = b
   b   = tmp
end subroutine swap_char
end subroutine swap_reaction!}}}

! reaction
subroutine reaction(T,vrho,tout,flag)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in)   ::T
   double precision,intent(inout)::vrho(ns)
   double precision,intent(in)   ::tout
   logical         ,intent(out)  ::flag

   double precision n(ns+1)
   double precision tt

   integer istate
   double precision,save::RWORK(LRW)
   integer         ,save::IWORK(LIW)
   !$omp threadprivate(RWORK,IWORK)
   integer          neq
   double precision rtol,atol
   integer          ipar(1)
   double precision rpar(1)
   external Fex,Jex

   integer,parameter::itol   =1  ! scalar atol
   integer,parameter::itask  =1  ! normal output
   integer,parameter::iopt   =0  ! optional input off
   integer,parameter::mf     =21 ! full matrix and direct jac.
   !integer,parameter::mf     =22 ! full matrix and non-direct jac.

   integer Nrecalc,j

   flag=.false.
   !vrho to n
   do j=1,ns
      n(j)=vrho(j)*invMW(j)*1d-3
      n(j)=max(n(j),n_eps)
   end do

   if(T    < 0d0) stop "negative temperature"
   if(tout < 0d0) stop "negative dt"

   !set parameters
   n(ns+1)= T
   tt      = 0d0
   istate  = 1
   neq     = ns+1
   rtol    = rtol_user
   atol    = 0d0
   Nrecalc = 0

   do
      call dvode (Fex, neq, n(1:ns+1), tt, tout, itol, rtol, atol, itask,  &
                  istate, iopt, rwork, lrw, iwork, liw, jex, mf,&
                  rpar, ipar)

      if(istate > 0) then
         exit
      else if(istate .eq. -1) then
         istate = 1
         atol=atol_user
         Nrecalc = Nrecalc+1
         if(Nrecalc>max_recalc) then
            print *,"Exceed max_recalc."
            return
         end if
      else
         print *,"odd istate=",istate
         return
      end if
   end do

   !!set T
   !T=n(ns+1)

   !n to vrho
   do j=1,ns
      n(j) = max(n(j),n_eps)
      vrho(j)=n(j)*MWs(j)*1d3
   end do
   flag=.true.
end subroutine reaction!}}}
subroutine reaction_cont(T,tt,tout,vrho,n,istate,RWORK,IWORK,atol,flag)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(inout)::T
   double precision,intent(inout)::tt
   double precision,intent(in)   ::tout
   double precision,intent(inout)::vrho(ns)
   double precision,intent(inout)::n(ns+1)
   integer         ,intent(inout)::istate
   double precision,intent(inout)::RWORK(LRW)
   integer         ,intent(inout)::IWORK(LIW)
   double precision,intent(inout)::atol
   logical         ,intent(out)  ::flag

   integer          neq,Nrecalc,j
   double precision rtol
   integer          ipar(1)
   double precision rpar(1)
   external Fex,Jex

   integer,parameter::itol   = 1  ! scalar atol
   integer,parameter::itask  = 1  ! normal output
   integer,parameter::iopt   = 1  ! optional input on
   integer,parameter::mf     = 21 ! full matrix and direct jac.


   flag=.false.
   if(istate .eq. 1) then
      !vrho to n
      do j=1,ns
         n(j)=vrho(j)*invMW(j)*1d-3
         n(j)=max(n(j),n_eps)
      end do
      n(ns+1) = T
      !atol    = 0d0 !cannot use to stiff problem
      atol    = atol_user
      IWORK   = 0
      RWORK   = 0d0
      IWORK(7)= 1
   end if

   !set parameters
   neq     = ns+1
   rtol    = rtol_user
   Nrecalc = 0
   do
      call dvode (Fex, neq, n(1:ns+1), tt, tout, itol, rtol, atol, itask,  &
                  istate, iopt, rwork, lrw, iwork, liw, jex, mf,&
                  rpar, ipar)
      if(IWORK(21)>0 .and. istate <0) then
         print *,"convergence test failed at VODE."
         return
      end if

      if(istate > 0) then
         exit
      else if(istate .eq. -1) then
         istate = 1
         atol=atol_user
         Nrecalc = Nrecalc+1
         if(Nrecalc>max_recalc) then
            print *,"Exceed max_recalc."
            return
         end if
      else
         stop
         print *,"Odd istate=",istate
         return
      end if
   end do

   !n to vrho
   do j=1,ns
      n(j)=max(n(j),n_eps)
      vrho(j)=n(j)*MWs(j)*1d3
   end do
   T=n(ns+1)
   flag=.true.
end subroutine reaction_cont!}}}
subroutine reaction_one_steps(T,vrho,mode,outarr,flag)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in)   ::T
   double precision,intent(in)   ::vrho(ns)
   integer         ,intent(in)   ::mode     ! 0:ignition / 1:Ymax / 2:CSP
   double precision,intent(inout)::outarr(*)
   logical         ,intent(out)  ::flag

   double precision n(ns+1)
   double precision tt,tout,ttold,Told,Tnow
   double precision tign,Teq

   double precision Y(ns),sm
   double precision atol,rtol

   integer istate
   double precision RWORK(LRW)
   integer          IWORK(LIW)
   integer          neq
   integer          ipar(1)
   double precision rpar(1)
   external Fex,Jex

   integer,parameter::itol   =1  ! scalar atol
   integer,parameter::itask  =5  ! ***************** only one step***********************
   integer,parameter::iopt   =0  ! optional input off
   integer,parameter::mf     =21 ! full matrix and direct jac.

   integer j,Nrecalc

   logical ignited

   flag=.false.
   !vrho to n
   do j=1,ns
      n(j)=vrho(j)*invMW(j)*1d-3
      n(j)=max(n(j),n_eps)
   end do
   !set parameters
   tout     = 100d0
   neq      = ns+1
   RWORK(1) = tout
   n(ns+1)  = T
   tt       = 0d0
   istate   = 1
   ignited  = .false.
   atol     = atol_user
   rtol     = rtol_user
   Nrecalc  = 0

   open(20,file="mono.out")
   do while(tt<tout)
      call dvode (Fex, neq, n(1:ns+1), tt, tout, itol, rtol, atol, itask,  &
                  istate, iopt, rwork, lrw, iwork, liw, jex, mf,&
                  rpar, ipar)
      write(20,*) tt,n(ns+1)

      if(istate < 0) then
         !print *, "istate=",istate
         close(20)
         return
      end if

      if(tt .eq. ttold) then
         !print *,"dt becomes too small to resolve by double precision."
         close(20)
         return
      end if

      if(.not. ignited .and. n(ns+1)>=T+400d0) then
         Tnow    = n(ns+1)
         tign    = tt-(Tnow-(T+400d0))/(Tnow-Told)*(tt-ttold)
         tout    = tign*10d0
         ignited = .true.
      end if
      Told =n(ns+1)
      ttold=tt

      Nrecalc=Nrecalc+1
      if(Nrecalc>500*max_recalc) then
         print *,"Exceed max_recalc."
         close(20)
         return
      end if

      if(mode .eq. 1) then
         Y=n(1:ns)*MWs(1:ns)
         sm=sum(Y,1)
         Y=Y/sm
         do j=1,ns
            if(Y(j)>outarr(j)) outarr(j)=Y(j)
         end do
      else if(mode .eq. 2) then
         call calc_CSP(n, outarr)
      end if
   end do
   close(20)

   if(.not.ignited) return
   Teq=n(ns+1)

   if(mode .eq. 0) then
      outarr(1)=tign
      outarr(2)=Teq
   end if
   flag=.true.
end subroutine reaction_one_steps!}}}

! cold flow
subroutine calc_T(vrho,rho,E, T, kappa,MWave,DHi,vhi,mu)!{{{
   use const_chem
   use chem
   use func_therm
   implicit none
   double precision,intent(in)::vrho(ns_tocalc)
   double precision,intent(in)::rho
   double precision,intent(in)::E
   double precision,intent(inout)::T
   double precision,intent(out)::kappa
   double precision,intent(out)::MWave
   double precision,intent(out)::DHi(ns)
   double precision,intent(out)::vhi(ns)
   double precision,intent(out)::mu

   double precision n(ns_tocalc)

   double precision Eobj,Enow,cvnow,de,dcv,dh,dT
   double precision sn,tmp
   double precision,dimension(6)::chrt,ccpr
   double precision,dimension(4)::cmu
   integer sec(ns)

   double precision muN(nt),denom,phi
   integer it,jt

   integer,parameter::max_T_loop=1000
   double precision,parameter::omega=1d-2

   integer i,j,k

   !vrho to n
   tmp = rho * 1d-11
   do j=1,ns_tocalc
      if(vrho(j)>tmp) then
         n(j)=vrho(j)*invMW(j)*1d-3
      else
         n(j)=0d0
      end if
   end do

   !temperature determination iteration
   Eobj = E*1d1*rho/Ru
   k=1
   do
      call check_section_number(T,sec)
      call set_coeff_less(T,log(T),chrt,ccpr)
      Enow  = 0d0
      cvnow = 0d0
      do j=1,ns_tocalc
         if(n(j)>0d0) then
            de  = -1d0
            dcv = -1d0
            do i=1,6
               de =de +coeff(i,sec(j),j)*chrt(i)
               dcv=dcv+coeff(i,sec(j),j)*ccpr(i)
            end do
            Enow =Enow +de *n(j)
            cvnow=cvnow+dcv*n(j)
         end if
      end do
      dT=(Eobj-Enow*T)/cvnow
      if(abs(dT)<1d-8 .or. (abs(dT)<1d-4 .and. abs(T-1d3) < 1d-4) .or. k>max_T_loop) exit
      T=max(100d0,T+dT)
      k=k+1
   end do

   !error detection of temperature determination iteration
   if(k>max_T_loop) then
      print '(a)',"Error: Exceed max_T_loop at temperature determination iteration"
      print '(a,es15.7)',"T=",T
      call exit(1)
   end if

   !specific heat ratio calculation
   sn=0d0
   do j=1,ns_tocalc
      sn=sn+n(j)
   end do
   kappa=1d0+sn/cvnow

   !calc molecular weight
   MWave=rho/sn*1d-3

   !calc DH
   tmp=Ru*1d-4*T
   do j=1,ns
      dh= 0d0
      do i=1,6
         dh =dh +coeff(i,sec(j),j)*chrt(i)
      end do
      DHi(j)=tmp*invMW(j)*(kappa-(kappa-1d0)*dh)
      vhi(j)=tmp*invMW(j)*dh
   end do

   !calc mu
   call check_section_number_trans(T,sec)
   call set_coeff_mu(T,log(T),cmu)
   do i=1,nt
      if(n(tr2th(i))>0d0) then
         tmp=0d0
         do j=1,4
            tmp=tmp+cmu(j)*trans(j,sec(i),i)
         end do
         muN(i) = exp(tmp)
      else
         muN(i) = 0d0
      end if
   end do

   mu = 0d0
   do i=1,nt
      if(muN(i)>0d0) then
         it=tr2th(i)
         denom = 0d0
         do j=1,nt
            if(muN(j)>0d0) then
               jt=tr2th(j)
               phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MWs(jt) / MWs(it) )**0.25d0 )**2 &
                            * sqrt(2d0 *MWs(jt)/(MWs(it)+MWs(jt)))
               denom = denom + n(jt) * phi
            end if
         end do
         mu = mu + n(it) * muN(i) / denom
      end if
   end do
   mu = mu * 1d-7

   !debug print
   !print '(a10, i15)',  "N itr=",k
   !print '(a10, f15.7)',"T=",T
   !print '(a10,es15.7)',"kappa=",kappa
   !print '(a10, f15.7)',"MWave=",MWave
   !print '(a10,es15.7)',"DH N2=",DHi(2)
   !print '(a10,es15.7)',"DH O2=",DHi(3)
   !print '(a10, f15.7)',"mu(mP)=",mu*1d4
   !print '(a10, f15.7)',"mu N2(mP)=",muN(10)*1d-3
   !print '(a10, f15.7)',"mu O2(mP)=",muN(13)*1d-3
end subroutine calc_T!}}}

!!for post process
!calculate dTdt
subroutine calc_dTdt(T,vrho,dTdt)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in)::T
   double precision,intent(in)::vrho(ns)
   double precision,intent(out)::dTdt

   double precision n(ns+1)
   double precision dndt(ns+1)
   double precision tt

   integer          neq
   integer          ipar(1)
   double precision rpar(1)
   external Fex

   integer j

   !vrho to n
   do j=1,ns
      n(j)=vrho(j)*invMW(j)*1d-3
      n(j) = max(n(j),n_eps)
   end do

   n(ns+1)= T
   neq         = ns+1
   tt=0d0

   call Fex(neq, tt, n, dndt, rpar, ipar)

   dTdt=dndt(ns+1)
end subroutine calc_dTdt!}}}

!chemeq2
subroutine reaction_chemeq2(T,vrho,dt,tout)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in)::T
   double precision,intent(inout)::vrho(ns)
   double precision,intent(in)::dt
   double precision,intent(in)::tout

   double precision n(ns+1)
   double precision tt,to
   integer num_recalc,i,j,nloop

   !vrho to n
   do j=1,ns
      n(j)=vrho(j)*invMW(j)*1d-3
      n(j)=max(n(j),n_eps)
   end do

   if(T    < 0d0) stop "negative temperature"
   if(tout < 0d0) stop "negative dt"

   !set parameters
   n(ns+1) = T
   tt=0d0;to=0d0

   nloop=tout/dt
   do i=1,nloop
      to=to+dt
      if(i .eq. nloop) to=tout
      call chemeq2(n,tt,to,atol_user,rtol_user)
   end do

   !n to vrho
   do j=1,ns
      n(j)=max(n(j),n_eps)
      vrho(j)=n(j)*MWs(j)*1d3
   end do
end subroutine reaction_chemeq2!}}}
subroutine chemeq2(y,tt,tout,atol,rtol)!{{{
   use chem
   implicit none
   double precision,intent(inout)  ::y(ns+1)
   double precision,intent(inout)  ::tt
   double precision,intent(in)     ::tout
   double precision,intent(in)     ::atol
   double precision,intent(in)     ::rtol

   !for chemeq2
   double precision,dimension(ns+1)::y0,yp,p,q,p0,q0
   double precision dt,ratio
   double precision var,var2,var3,alp
   double precision counter
   integer i,countersum
   external Fex_var

   !for dvode
   integer istate,Nrecalc
   double precision RWORK(LRW)
   integer          IWORK(LIW)
   integer          ipar(1)
   double precision rpar(1)
   external Fex,Jex
   integer,parameter::itol  = 1  ! scalar atol
   integer,parameter::itask = 1  ! normal output
   integer,parameter::iopt  = 0  ! optional input off
   integer,parameter::mf    = 21 ! full matrix and direct jac.


   dt=tout-tt
   countersum=0
   do
      y0=y
      counter=0d0
      if(tout-tt<dt) dt=tout-tt
      do
         y=y0
         countersum=countersum+1

         !move to dvode
         if(countersum>10) then
            istate=1
            Nrecalc=0
            do
               call dvode (Fex, ns+1, y, tt, tout, itol, rtol, atol, itask,  &
                        istate, iopt, rwork, lrw, iwork, liw, jex, mf,&
                        rpar, ipar)
               if(istate > 0) then
                  exit
               else if(istate .eq. -1) then
                  istate = 1
                  Nrecalc = Nrecalc+1
                  if(Nrecalc>max_recalc) stop "Exceed max_recalc."
               else
                  print *,"odd istate=",istate
                  stop
               end if
            end do
            dt=0d0
            exit
         end if

         !predictor
         call Fex_var(y,q0,p0)
         do i=1,ns
            p0(i)=-p0(i)/(y(i)+1d-300)

            var=p0(i)*dt
            var2=var*var;var3=var*var2
            alp=(180d0+60d0*var+11d0*var2+var3)/(360d0+60d0*var+12d0*var2+var3)
            yp(i)=y0(i)+dt*(q0(i)-p0(i)*y0(i))/(1d0+alp*p0(i)*dt)
         end do
         yp(ns+1) =y0(ns+1)+dt*(q0(ns+1)+p0(ns+1))
         
         !corrector
         ratio=1d300 !init residual
         call Fex_var(yp,q,p)
         do i=1,ns
            p(i) =-p(i)/(yp(i)+1d-300)
            p(i) =0.5d0*(p0(i)+p(i))

            var=p(i)*dt
            var2=var*var; var3=var*var2
            alp=(180d0+60d0*var+11d0*var2+var3)/(360d0+60d0*var+12d0*var2+var3)
            q(i) =alp*q(i)+(1d0-alp)*q0(i)
            y(i) =y0(i)+dt*(q(i) -p(i) *y0(i))/(1d0+alp*p(i) *dt)

            !calc residual
            var = (y(i)*rtol+rtol*1d-5)/abs(y(i)-yp(i)+1d-300)
            ratio=min(ratio,var)
         end do
         y(ns+1) =y0(ns+1)+0.5d0*dt*(q(ns+1)+p(ns+1)+q0(ns+1)+p0(ns+1))

         !finalize
         if(1d0<ratio) exit
         counter=counter+1d0
         dt=dt*sqrt(ratio)/counter
         if(counter>10d0) stop "not converted at CHEMEQ2"
      end do
      tt=tt+dt
      write(20,'(2es15.7,es9.1,i5)') tt,y(ns+1),dt,countersum
      if(tout <= tt) exit
      dt=dt*sqrt(ratio)
   end do
end subroutine chemeq2!}}}
