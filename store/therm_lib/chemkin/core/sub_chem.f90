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

         if(linebuf .eq. "duplicate" .or. linebuf .eq. "DUPLICATE") cycle

         ind=index(linebuf,"/")
         if(ind .eq. 0) then ! reaction declaration line
            nr=nr+1
            if(nr > max_nr) stop "chem.inp reaction number exceeded assigned one."

            !read default ABE
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
   
      do i=1,ns_tocalc
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

      do i=1,NumM(k)
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
subroutine calc_CSP(n, csp_out)!{{{
   !variables
   use const_chem
   use chem
   use func_therm
   implicit none
   double precision,intent(in) :: n(ns+1)
   double precision,intent(inout):: csp_out(nr)

   double precision csp(ns,nr)

   double precision,parameter::dd=0.14d0
   double precision,parameter::loge210=0.43429448190325!log10(exp(1d0))

   double precision T
   double precision,dimension(ns)::vmu0rt
   double precision vlogM,Dr,r
   double precision vlogk(2)
   double precision de,cv,sn,logsn
   double precision logT,Tinv
   double precision logkinf,logk0
   double precision Pr
   double precision logFcent,logF,cc,nn,aa,log_Pr,logPrc,logPRT
   double precision temp1,temp2,temp3,tmp
   double precision,dimension(7)::cmurt,chrt,ccpr
   integer          sec(ns)
   integer i,j,k,itmp

   csp=0d0

   T    =n(ns+1)
   logT=log(T)
   Tinv=1d0/T
   logPRT=log(pst/Ru)-logT
   sn   =sum(n(1:ns))
   logsn=log(sn)

   !set thermodynamical values
   call check_section_number(T,sec)
   call set_coeff(T,logT,cmurt,chrt,ccpr)
   do j=1,ns
      temp1=0d0
      do i=1,7
         temp1=temp1+coeff(i,sec(j),j)*cmurt(i)
      end do
      vmu0rt(j)=temp1
   end do

   do i=1,nr
      !calc M
      if(exist_M(i)) then
         tmp =sn
         do j=1,NumM(i)
            tmp=tmp+Men(j,i)*n(IndM(j,i))
         end do

         if(tmp .eq. sn) then
            vlogM=logsn
         else
            vlogM=log(tmp)
         end if
      end if

      !calc k
      vlogk(1)=ABE(1,i)+ABE(2,i)*logT-ABE(3,i)*Tinv
      vlogk(2)=0d0
      if(Rstate(i,2) >= 4) vlogk(1)=vlogk(1)+vlogM

      if(Rstate(i,1) .eq. 0) then !bi-direction
         if(Rstate(i,2) .eq. 1 .or. Rstate(i,2) .eq. 5) then !rev
            vlogk(2)=cABE(1,i)+cABE(2,i)*logT-cABE(3,i)*Tinv !read cABE
            if(Rstate(i,2) .eq. 5) vlogk(2)=vlogk(2)+vlogM
         else !none, low, troe
            if(Rstate(i,2) .eq. 2 .or. Rstate(i,2) .eq. 3) then!low, troe
               !Pr
               logk0  =cABE(1,i)+cABE(2,i)*logT-cABE(3,i)*Tinv
               logkinf=vlogk(1)
               log_Pr=logk0-logkinf+vlogM
               Pr=exp(log_Pr)
               vlogk(1)=logkinf+log(Pr/(1d0+Pr))

               !Troe
               if(Rstate(i,2) .eq. 3) then
                  aa=TROE(1,i)
                  logFcent=log10((1d0-aa)*exp(-T   *TROE(2,i))&
                                +aa      *exp(-T   *TROE(3,i))&
                                +         exp(-Tinv*TROE(4,i)))
                  cc=-0.4d0-0.67d0*logFcent
                  nn=0.75d0-1.27d0*logFcent
                  logPrc=log_Pr*loge210+cc
                  logF=(1d0+(logPrc/(nn-dd*logPrc))**2)**(-1)*logFcent
                  vlogk(1)=vlogk(1)+logF/loge210
               end if
            end if

            !calc numu0rt
            temp1=0d0
            do k=1,NumNu(1,i)
               temp1=temp1-vmu0rt(IndNu(k,1,i))
            end do
            do k=1,NumNu(2,i)
               temp1=temp1+vmu0rt(IndNu(k,2,i))
            end do

            vlogk(2)=vlogk(1)-snu(i)*logPRT+temp1
         end if
      end if

      !calc r
      r=exp(vlogk(1))
      do k=1,NumNu(1,i)
         r=r*n(IndNu(k,1,i))
      end do
      Dr=r
      if(Rstate(i,1) .eq. 0) then
         r=exp(vlogk(2))
         do k=1,NumNu(2,i)
            r=r*n(IndNu(k,2,i))
         end do
         Dr=Dr-r
      end if
      Dr=abs(Dr)

      !calc csp
      do k=1,NumNu(1,i)
         itmp=IndNu(k,1,i)
         csp(itmp,i)=csp(itmp,i)+Dr
      end do
      do k=1,NumNu(2,i)
         itmp=IndNu(k,2,i)
         csp(itmp,i)=csp(itmp,i)+Dr
      end do
   end do

   do j=1,ns
      tmp=0d0
      do i=1,nr
         tmp=tmp+csp(j,i)
      end do
      tmp=1d0/(tmp+1d-300)
      do i=1,nr
         csp(j,i)=csp(j,i)*tmp
      end do
   end do
   do i=1,nr
      do j=1,ns
         csp_out(i)=max(csp_out(i),csp(j,i))
      end do
   end do
end subroutine calc_CSP!}}}
subroutine Fex(neq_outer, tt, n, dndt, rpar, ipar)!{{{
   !variables
   use const_chem
   use chem
   use func_therm
   implicit none
   integer         ,intent(in) :: neq_outer   !not used now.
   double precision,intent(in) :: tt
   double precision,intent(in) :: n(ns+1)
   double precision,intent(out):: dndt(ns+1)
   double precision,intent(inout):: rpar(*)
   integer         ,intent(inout):: ipar(*)

   double precision,parameter::dd=0.14d0
   double precision,parameter::loge210=0.43429448190325!log10(exp(1d0))

   double precision T
   double precision,dimension(ns)::vmu0rt,vert,vcvr
   double precision vlogM,Dr,r
   double precision vlogk(2)
   double precision de,cv,sn,logsn
   double precision logT,Tinv
   double precision logkinf,logk0
   double precision Pr
   double precision logFcent,logF,cc,nn,aa,log_Pr,logPrc,logPRT
   double precision temp1,temp2,temp3,tmp
   double precision,dimension(7)::cmurt,chrt,ccpr
   integer          sec(ns)
   integer i,j,k,itmp

   T    =n(ns+1)
   logT=log(T)
   Tinv=1d0/T
   logPRT=log(pst/Ru)-logT
   sn   =sum(n(1:ns))
   logsn=log(sn)

   !set thermodynamical values
   call check_section_number(T,sec)
   call set_coeff(T,logT,cmurt,chrt,ccpr)
   do j=1,ns
      temp1=0d0
      temp2=-1d0
      temp3=-1d0
      do i=1,7
         tmp=coeff(i,sec(j),j)
         temp1=temp1+tmp*cmurt(i)
         temp2=temp2+tmp*chrt( i)
         temp3=temp3+tmp*ccpr( i)
      end do
      vmu0rt(j)=temp1
      vert(  j)=temp2
      vcvr(  j)=temp3

      dndt(  j)=0d0
   end do

   do i=1,nr
      !calc M
      if(exist_M(i)) then
         tmp =sn
         do j=1,NumM(i)
            tmp=tmp+Men(j,i)*n(IndM(j,i))
         end do

         if(tmp .eq. sn) then
            vlogM=logsn
         else
            vlogM=log(tmp)
         end if
      end if

      !calc k
      vlogk(1)=ABE(1,i)+ABE(2,i)*logT-ABE(3,i)*Tinv
      vlogk(2)=0d0
      if(Rstate(i,2) >= 4) vlogk(1)=vlogk(1)+vlogM

      if(Rstate(i,1) .eq. 0) then !bi-direction
         if(Rstate(i,2) .eq. 1 .or. Rstate(i,2) .eq. 5) then !rev
            vlogk(2)=cABE(1,i)+cABE(2,i)*logT-cABE(3,i)*Tinv !read cABE
            if(Rstate(i,2) .eq. 5) vlogk(2)=vlogk(2)+vlogM
         else !none, low, troe
            if(Rstate(i,2) .eq. 2 .or. Rstate(i,2) .eq. 3) then!low, troe
               !Pr
               logk0  =cABE(1,i)+cABE(2,i)*logT-cABE(3,i)*Tinv
               logkinf=vlogk(1)
               log_Pr=logk0-logkinf+vlogM
               Pr=exp(log_Pr)
               vlogk(1)=logkinf+log(Pr/(1d0+Pr))

               !Troe
               if(Rstate(i,2) .eq. 3) then
                  aa=TROE(1,i)
                  logFcent=log10((1d0-aa)*exp(-T   *TROE(2,i))&
                                +aa      *exp(-T   *TROE(3,i))&
                                +         exp(-Tinv*TROE(4,i)))
                  cc=-0.4d0-0.67d0*logFcent
                  nn=0.75d0-1.27d0*logFcent
                  logPrc=log_Pr*loge210+cc
                  logF=(1d0+(logPrc/(nn-dd*logPrc))**2)**(-1)*logFcent
                  vlogk(1)=vlogk(1)+logF/loge210
               end if
            end if

            !calc numu0rt
            temp1=0d0
            do k=1,NumNu(1,i)
               temp1=temp1-vmu0rt(IndNu(k,1,i))
            end do
            do k=1,NumNu(2,i)
               temp1=temp1+vmu0rt(IndNu(k,2,i))
            end do

            vlogk(2)=vlogk(1)-snu(i)*logPRT+temp1
         end if
      end if

      !calc r
      r=exp(vlogk(1))
      do k=1,NumNu(1,i)
         r=r*n(IndNu(k,1,i))
      end do
      Dr=r
      if(Rstate(i,1) .eq. 0) then
         r=exp(vlogk(2))
         do k=1,NumNu(2,i)
            r=r*n(IndNu(k,2,i))
         end do
         Dr=Dr-r
      end if

      !calc dndt
      do k=1,NumNu(1,i)
         itmp=IndNu(k,1,i)
         dndt(itmp)=dndt(itmp)-Dr
      end do
      do k=1,NumNu(2,i)
         itmp=IndNu(k,2,i)
         dndt(itmp)=dndt(itmp)+Dr
      end do
   end do

   de=0d0; cv=0d0
   do j=1,ns
      de=de+vert(j)*dndt(j)
      cv=cv+vcvr(j)*n(   j)
   end do
   dndt(ns+1)=-de*T/cv
end subroutine Fex!}}}
subroutine Fex_var(n, dndtp, dndtm)!{{{
   !variables
   use const_chem
   use chem
   use func_therm
   implicit none
   double precision,intent(in) :: n(ns+1)
   double precision,intent(out):: dndtp(ns+1)
   double precision,intent(out):: dndtm(ns+1)

   double precision,parameter::dd=0.14d0
   double precision,parameter::loge210=0.43429448190325!log10(exp(1d0))

   double precision T
   double precision,dimension(ns)::vmu0rt,vert,vcvr
   double precision vlogM,Dr,r
   double precision vlogk(2)
   double precision de,cv,sn,logsn
   double precision logT,Tinv
   double precision logkinf,logk0
   double precision Pr
   double precision logFcent,logF,cc,nn,aa,log_Pr,logPrc,logPRT
   double precision temp1,temp2,temp3,tmp
   double precision,dimension(7)::cmurt,chrt,ccpr
   integer          sec(ns)
   integer i,j,k,itmp

   T    =n(ns+1)
   logT=log(T)
   Tinv=1d0/T
   logPRT=log(pst/Ru)-logT
   sn   =sum(n(1:ns))
   logsn=log(sn)

   !set thermodynamical values
   call check_section_number(T,sec)
   call set_coeff(T,logT,cmurt,chrt,ccpr)
   do j=1,ns
      temp1=0d0
      temp2=-1d0
      temp3=-1d0
      do i=1,7
         tmp=coeff(i,sec(j),j)
         temp1=temp1+tmp*cmurt(i)
         temp2=temp2+tmp*chrt( i)
         temp3=temp3+tmp*ccpr( i)
      end do
      vmu0rt(j)=temp1
      vert(  j)=temp2
      vcvr(  j)=temp3

      dndtp(  j)=0d0
      dndtm(  j)=0d0
   end do

   do i=1,nr
      !calc M
      if(exist_M(i)) then
         tmp =sn
         do j=1,NumM(i)
            tmp=tmp+Men(j,i)*n(IndM(j,i))
         end do

         if(tmp .eq. sn) then
            vlogM=logsn
         else
            vlogM=log(tmp)
         end if
      end if

      !calc k
      vlogk(1)=ABE(1,i)+ABE(2,i)*logT-ABE(3,i)*Tinv
      vlogk(2)=0d0
      if(Rstate(i,2) >= 4) vlogk(1)=vlogk(1)+vlogM

      if(Rstate(i,1) .eq. 0) then !bi-direction
         if(Rstate(i,2) .eq. 1 .or. Rstate(i,2) .eq. 5) then !rev
            vlogk(2)=cABE(1,i)+cABE(2,i)*logT-cABE(3,i)*Tinv !read cABE
            if(Rstate(i,2) .eq. 5) vlogk(2)=vlogk(2)+vlogM
         else !none, low, troe
            if(Rstate(i,2) .eq. 2 .or. Rstate(i,2) .eq. 3) then!low, troe
               !Pr
               logk0  =cABE(1,i)+cABE(2,i)*logT-cABE(3,i)*Tinv
               logkinf=vlogk(1)
               log_Pr=logk0-logkinf+vlogM
               Pr=exp(log_Pr)
               vlogk(1)=logkinf+log(Pr/(1d0+Pr))

               !Troe
               if(Rstate(i,2) .eq. 3) then
                  aa=TROE(1,i)
                  logFcent=log10((1d0-aa)*exp(-T   *TROE(2,i))&
                                +aa      *exp(-T   *TROE(3,i))&
                                +         exp(-Tinv*TROE(4,i)))
                  cc=-0.4d0-0.67d0*logFcent
                  nn=0.75d0-1.27d0*logFcent
                  logPrc=log_Pr*loge210+cc
                  logF=(1d0+(logPrc/(nn-dd*logPrc))**2)**(-1)*logFcent
                  vlogk(1)=vlogk(1)+logF/loge210
               end if
            end if

            !calc numu0rt
            temp1=0d0
            do k=1,NumNu(1,i)
               temp1=temp1-vmu0rt(IndNu(k,1,i))
            end do
            do k=1,NumNu(2,i)
               temp1=temp1+vmu0rt(IndNu(k,2,i))
            end do

            vlogk(2)=vlogk(1)-snu(i)*logPRT+temp1
         end if
      end if

      !calc r
      r=exp(vlogk(1))
      do k=1,NumNu(1,i)
         r=r*n(IndNu(k,1,i))
      end do
      do k=1,NumNu(1,i)
         itmp=IndNu(k,1,i)
         dndtm(itmp)=dndtm(itmp)-r
      end do
      do k=1,NumNu(2,i)
         itmp=IndNu(k,2,i)
         dndtp(itmp)=dndtp(itmp)+r
      end do

      if(Rstate(i,1) .eq. 0) then
         r=exp(vlogk(2))
         do k=1,NumNu(2,i)
            r=r*n(IndNu(k,2,i))
         end do
         do k=1,NumNu(1,i)
            itmp=IndNu(k,1,i)
            dndtp(itmp)=dndtp(itmp)+r
         end do
         do k=1,NumNu(2,i)
            itmp=IndNu(k,2,i)
            dndtm(itmp)=dndtm(itmp)-r
         end do
      end if
   end do

   de=0d0; cv=0d0
   do j=1,ns
      de=de+vert(j)*(dndtp(j)+dndtm(j))
      cv=cv+vcvr(j)*n(   j)
   end do
   dndtp(ns+1)=max(-de*T/cv,0d0)
   dndtm(ns+1)=min(-de*T/cv,0d0)
end subroutine Fex_var!}}}
subroutine Jex(neq_outer, tt, n, ML, MU, dndn, nrpd, rpar, ipar)!{{{
   !variables
   use const_chem
   use chem
   use func_therm
   implicit none
   integer         ,intent(in) :: neq_outer   !not used now.
   double precision,intent(in) :: tt
   double precision,intent(in) :: n(ns+1)
   integer         ,intent(in) :: ML          !not used now.
   integer         ,intent(in) :: MU          !not used now.
   double precision,intent(out):: dndn(ns+1,ns+1)
   integer         ,intent(in) :: nrpd        !not used now.
   double precision,intent(inout):: rpar(*)
   integer         ,intent(inout):: ipar(*)

   double precision,parameter::dd=0.14d0
   double precision,parameter::loge210=0.43429448190325!log10(exp(1d0))

   double precision T
   double precision dndt(ns+1)
   double precision drdn(ns+1)
   double precision,dimension(ns)::vmu0rt,vdmu0rt,vert,vcvr,vdcvrdT
   double precision numu0rt,nudmu0rt,vlogM,vM,dlnkdlnM
   double precision vlogk(2),dlnkdT(2),r
   double precision Dr
   double precision de,cv,dcv,logsn,sn
   double precision logT,Tinv
   double precision logkinf,logk0
   double precision dlnkinfdT,dlnk0dT
   double precision Pr,dlnPrdT
   double precision logFcent,logF,Fcent,cc,nn,aa,log_Pr,logPrc,logPRT,dlnFcentdT,dFcentdT
   double precision temp1,temp2,temp3,temp4,temp5,tmp
   integer i,j,k,itmp
   double precision,dimension(7)::cmurt,chrt,ccpr,cdmurt,cdcpr
   integer          sec(ns)

   T     =n(ns+1)
   logT  =log(T)
   Tinv  =1d0/T
   logPRT=log(pst/Ru)-logT
   sn    =max(sum(n(1:ns)),1d-300)
   logsn =log(sn)

   !set thermodynamical values
   call check_section_number(T,sec)
   call set_coeff_more(T,logT,cmurt,chrt,ccpr,cdmurt,cdcpr)
   do j=1,ns
      temp1= 0d0
      temp2= 0d0
      temp3=-1d0
      temp4=-1d0
      temp5= 0d0
      do i=1,7
         tmp=coeff(i,sec(j),j)
         temp1=temp1+tmp* cmurt(i)
         temp2=temp2+tmp*cdmurt(i)
         temp3=temp3+tmp*  chrt(i)
         temp4=temp4+tmp*  ccpr(i)
         temp5=temp5+tmp* cdcpr(i)
      end do
      vmu0rt( j)=temp1
      vdmu0rt(j)=temp2
      vert(   j)=temp3
      vcvr(   j)=temp4
      vdcvrdT(j)=temp5

      dndt(j)=0d0
      do i=1,ns+1
         dndn(j,i)=0d0
      end do
   end do

   do i=1,nr
      !calc M
      if(exist_M(i)) then
         tmp =sn
         do j=1,NumM(i)
            tmp=tmp+Men(j,i)*n(IndM(j,i))
         end do

         vM   =tmp
         if(tmp .eq. sn) then
            vlogM=logsn
         else
            vlogM=log(tmp)
         end if
      end if

      !set k
      vlogk (1)= ABE(1,i)+ABE(2,i)*logT-ABE(3,i)*Tinv
      dlnkdT(1)=(ABE(2,i)+ABE(3,i)*Tinv)*Tinv
      vlogk (2)=0d0
      dlnkdT(2)=0d0
      if(Rstate(i,2) >= 4) then
         vlogk(1)=vlogk(1)+vlogM
         dlnkdlnM=1d0
      end if

      if(Rstate(i,1) .eq. 0) then !bi-direction
         if(Rstate(i,2) .eq. 1 .or. Rstate(i,2) .eq. 5) then !rev
            vlogk (2)=cABE(1,i)+cABE(2,i)*logT-cABE(3,i)*Tinv !read cABE
            dlnkdT(2)=(cABE(2,i)+cABE(3,i)*Tinv)*Tinv
            if(Rstate(i,2) .eq. 5) vlogk(2)=vlogk(2)+vlogM
         else !none, low, troe
            if(Rstate(i,2) .eq. 2 .or. Rstate(i,2) .eq. 3) then!low, troe
               !Pr
               logk0  =cABE(1,i)+cABE(2,i)*logT-cABE(3,i)*Tinv
               dlnk0dT=(cABE(2,i)+cABE(3,i)*Tinv)*Tinv
               logkinf=vlogk(1)
               dlnkinfdT=dlnkdT(1)
               log_Pr=logk0-logkinf+vlogM
               Pr=exp(log_Pr)
               vlogk (1)=logkinf+log(Pr/(1d0+Pr))
               dlnkdT(1)=(dlnk0dT+Pr*dlnkinfdT)/(1d0+Pr)
               dlnkdlnM=1d0/(1d0+Pr)

               !Troe
               if(Rstate(i,2) .eq. 3) then
                  dlnPrdT=dlnk0dT-dlnkinfdT
                  aa=TROE(1,i)
                  Fcent   = (1d0-aa)          *exp(-T*TROE(2,i))&
                           +aa                *exp(-T*TROE(3,i))&
                           +                   exp(-TROE(4,i)*Tinv)
                  dFcentdT=-(1d0-aa)*TROE(2,i)*exp(-T*TROE(2,i))&
                           -aa*TROE(3,i)      *exp(-T*TROE(3,i))&
                           +TROE(4,i)*Tinv**2 *exp(-TROE(4,i)*Tinv)
                  dlnFcentdT=dFcentdT/Fcent
                  logFcent=log10(Fcent)
                  cc=-0.4d0-0.67d0*logFcent
                  nn=0.75d0-1.27d0*logFcent
                  logPrc=log_Pr*loge210+cc
                  temp1=logPrc/(nn-dd*logPrc)
                  temp2=(1d0+temp1**2)**(-1)
                  logF=temp2*logFcent
                  vlogk( 1)=vlogk( 1)+logF/loge210
                  dlnkdT(1)=dlnkdT(1)&
                           +temp2*(-2d0*logF*temp1/(nn-dd*logPrc)**2 &
                                   *(nn*dlnPrdT+(1.27d0*logPrc-0.67d0*nn)*dlnFcentdT)&         
                                   +dlnFcentdT)                                                
                  dlnkdlnM =dlnkdlnM-2d0*temp2*logF*temp1/(nn-dd*logPrc)**2*nn
               end if
            end if

            !set vnumu0rt
            numu0rt =0d0
            nudmu0rt=0d0
            do k=1,NumNu(1,i)
               itmp=IndNu(k,1,i)
               numu0rt =numu0rt - vmu0rt(itmp)
               nudmu0rt=nudmu0rt-vdmu0rt(itmp)
            end do
            do k=1,NumNu(2,i)
               itmp=IndNu(k,2,i)
               numu0rt =numu0rt + vmu0rt(itmp)
               nudmu0rt=nudmu0rt+vdmu0rt(itmp)
            end do

            vlogk (2)=vlogk (1)-snu(i)*logPRT+numu0rt
            dlnkdT(2)=dlnkdT(1)+snu(i)*Tinv  +nudmu0rt
         end if
      end if

      !!!!!! set r and drdn !!!!!!
      !initialize
      do k=1,ns
         drdn(k)=0d0
      end do
      !forward
      r=exp(vlogk(1))
      do k=1,NumNu(1,i)
         r=r*n(IndNu(k,1,i))
      end do
      Dr        =          r
      drdn(ns+1)=dlnkdT(1)*r
      do k=1,NumNu(1,i)
         itmp=IndNu(k,1,i)
         drdn(itmp)=drdn(itmp)+r/n(itmp)
      end do

      !reverse
      if(Rstate(i,1) .eq. 0) then
         r=exp(vlogk(2))
         do k=1,NumNu(2,i)
            r=r*n(IndNu(k,2,i))
         end do
         Dr        =Dr        -r
         drdn(ns+1)=drdn(ns+1)-r*dlnkdT(2)
         do k=1,NumNu(2,i)
            itmp=IndNu(k,2,i)
            drdn(itmp)=drdn(itmp)-r/n(itmp)
         end do
      end if

      !about M
      if(exist_M(i)) then
         temp1=Dr*dlnkdlnM/vM
         do k=1,ns
            drdn(k)=drdn(k)+temp1
         end do
         do k=1,NumM(i)
            itmp=IndM(k,i)
            drdn(itmp)=drdn(itmp)+temp1*Men(k,i)
         end do
      end if

      !calc dndt,dndn
      do k=1,NumNu(1,i)
         itmp=IndNu(k,1,i)
         dndt(itmp)=dndt(itmp)-Dr
         do j=1,ns+1
            dndn(itmp,j)=dndn(itmp,j)-drdn(j)
         end do
      end do

      do k=1,NumNu(2,i)
         itmp=IndNu(k,2,i)
         dndt(itmp)=dndt(itmp)+Dr
         do j=1,ns+1
            dndn(itmp,j)=dndn(itmp,j)+drdn(j)
         end do
      end do
   end do !!!!!! END OF I LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !set dTdt and dTdT
   de=0d0; cv=0d0; dcv=0d0
   temp1=0d0;temp2=0d0
   do j=1,ns
      temp1=temp1+   vert(j)*dndn(j,ns+1)
      de   =de   +   vert(j)*dndt(j)
      dcv  =dcv  +   vcvr(j)*dndt(j)
      cv   =cv   +   vcvr(j)*   n(j)
      temp2=temp2+vdcvrdT(j)*   n(j)
   end do
   dndt(ns+1)=-de*T/cv
   dndn(ns+1,ns+1)=-(T*temp1+dcv+dndt(ns+1)*temp2)/cv

   !dTdn
   do j=1,ns
      temp1=0d0
      do k=1,ns
         temp1=temp1+vert(k)*dndn(k,j)
      end do
      dndn(ns+1,j)=-(T*temp1+dndt(ns+1)*vcvr(j))/cv
   end do

   !do j=1,ns+1
   !   do i=1,ns+1
   !      if(isNaN(dndn(i,j))) then
   !         print *,i,j,"NaN!"
   !         call exit(1)
   !      end if
   !   end do
   !end do
end subroutine Jex!}}}

subroutine reaction(T,vrho,tout)!{{{
   use const_chem
   use chem
   implicit none
   double precision,intent(in)::T
   double precision,intent(inout)::vrho(ns)
   double precision,intent(in)::tout

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

   integer num_recalc,j

   !vrho to n
   do j=1,ns
      n(j)=vrho(j)*invMW(j)*1d-3
      n(j)=max(n(j),n_eps)
   end do

   if(T    < 0d0) stop "negative temperature"
   if(tout < 0d0) stop "negative dt"

   !set parameters
   n(ns+1)= T
   tt          = 0d0
   istate      = 1
   neq         = ns+1
   rtol        = rtol_user
   atol        = 0d0
   num_recalc  = 0

   do
      call dvode (Fex, neq, n(1:ns+1), tt, tout, itol, rtol, atol, itask,  &
                  istate, iopt, rwork, lrw, iwork, liw, jex, mf,&
                  rpar, ipar)

      if(istate > 0) then
         exit
      else if(istate < 0) then
         istate = 1
         atol=atol_user
         num_recalc = num_recalc+1
         if(num_recalc>max_recalc) then
            print *,"Exceed max_recalc=",max_recalc," istate = ",istate
            call exit(1)
         end if
      end if
   end do

   !!set T
   !T=n(ns+1)

   !n to vrho
   do j=1,ns
      n(j) = max(n(j),n_eps)
      vrho(j)=n(j)*MWs(j)*1d3
   end do
end subroutine reaction!}}}
subroutine reaction_cont(T,tt,tout,vrho,n,istate,RWORK,IWORK,atol)!{{{
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

   integer          neq,Nrecalc,j
   double precision rtol
   integer          ipar(1)
   double precision rpar(1)
   external Fex,Jex

   integer,parameter::itol   =1  ! scalar atol
   integer,parameter::itask  =1  ! normal output
   integer,parameter::iopt   =0  ! optional input off
   integer,parameter::mf     =21 ! full matrix and direct jac.


   if(istate .eq. 1) then
      !vrho to n
      do j=1,ns
         n(j)=vrho(j)*invMW(j)*1d-3
         n(j)=max(n(j),n_eps)
      end do
      n(ns+1) = T
      atol    = 0d0
   end if

   !set parameters
   neq     = ns+1
   rtol    = rtol_user
   Nrecalc = 0
   do
      call dvode (Fex, neq, n(1:ns+1), tt, tout, itol, rtol, atol, itask,  &
                  istate, iopt, rwork, lrw, iwork, liw, jex, mf,&
                  rpar, ipar)

      if(istate > 0) then
         exit
      else if(istate < 0) then
         istate = 1
         atol=atol_user
         Nrecalc = Nrecalc+1
         if(Nrecalc>max_recalc) then
            print *,"Exceed max_recalc=",max_recalc," istate = ",istate
            call exit(1)
         end if
      end if
   end do

   !n to vrho
   do j=1,ns
      n(j)=max(n(j),n_eps)
      vrho(j)=n(j)*MWs(j)*1d3
   end do
   T=n(ns+1)
end subroutine reaction_cont!}}}
subroutine plot_time_history(T,vrho,tout,Nloop)!{{{
   use const_chem
   implicit none
   double precision,intent(inout)::T
   double precision,intent(inout)::vrho(ns)
   double precision,intent(in)   ::tout
   integer,         intent(inout)::Nloop

   double precision tt,to
   double precision n(ns+1)
   double precision RWORK(LRW)
   integer          IWORK(LIW)
   integer          istate,i
   double precision atol

   tt    =0d0
   istate=1

   open(22,file="plt.dat")
   write(22,'(2es15.7)') tt,T
   do i=1,Nloop
      to=dble(i)/dble(Nloop)*tout
      call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
      write(22,'(2es15.7)') tt,T
   end do
   close(22)
end subroutine plot_time_history!}}}
subroutine calc_Tign_Teq(Tini,vrhoini,allowable_limit,tign,Teq)!{{{
   use const_chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::allowable_limit
   double precision,intent(out)::tign
   double precision,intent(out)::Teq

   double precision T
   double precision vrho(ns)

   double precision csp(nr)

   double precision tt,to,tu,tl
   double precision n(ns+1)
   double precision RWORK(LRW)
   integer          IWORK(LIW)
   integer          istate,i
   double precision atol

   logical flag

   T   =Tini
   vrho=vrhoini

   tt    =0d0
   istate=1
   to=1d-10
   call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
   if(tt<to) stop "Ignition occurred in 1e-10 sec. Too small."
   do
      to=to*1d2
      call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
      if(T>Tini+400d0) exit
      if(to>1d2) stop "Ignition delay is larger than 1e2 sec. Too large."
   end do
   tu=to
   tl=to*1d-2

   tt    =0d0
   istate=1
   T   =Tini
   vrho=vrhoini
   to=exp((log(tu)+log(tl))/2d0)
   do while((tu-tl)/(tu+tl)*2d0>allowable_limit*1d-1)
      call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
      if(T>Tini+400d0) then
         tu=to
         tt    =0d0
         istate=1
         T   =Tini
         vrho=vrhoini
      else
         tl=to
      end if
      to=exp((log(tu)+log(tl))/2d0)
   end do
   tign=tu

   tt    =0d0
   istate=1
   T   =Tini
   vrho=vrhoini
   to=10d0*tign
   call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
   Teq=T
   call check_model(Tini,vrhoini,tign,Teq,allowable_limit,flag)
   if(.not.flag) stop "Consistency check failed at calc_Tign_Teq."
end subroutine calc_Tign_Teq!}}}
subroutine check_model(Tini,vrhoini,tign,Teq,allowable_limit,flag)!{{{
   use chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::tign
   double precision,intent(in)::Teq
   double precision,intent(in)::allowable_limit
   logical         ,intent(out)::flag

   double precision T
   double precision vrho(ns)

   double precision tt,to
   double precision n(ns+1)
   double precision RWORK(LRW)
   integer          IWORK(LIW)
   integer          istate
   double precision atol

   T   =Tini
   vrho=vrhoini

   flag   = .false.
   tt     = 0d0
   istate = 1
   to=tign*(1d0-allowable_limit)
   call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
   if(T>=Tini+400d0) return
   to=tign*(1d0+allowable_limit)
   call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
   if(T<Tini+400d0) return
   to=tign*10d0
   call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
   if(abs(T-Teq)>Teq*allowable_limit) return
   flag=.true.
end subroutine check_model!}}}
subroutine sort_n(Tini,vrhoini,tign,flag)!{{{
   use chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::tign
   logical         ,intent(in)::flag

   double precision T
   double precision vrho(ns)
   double precision Y(ns),Ymax(ns)

   double precision tt,to
   double precision n(ns+1)
   double precision RWORK(LRW)
   integer          IWORK(LIW)
   integer          istate,i,j
   double precision atol

   double precision sm,dtemp

   T   =Tini
   vrho=vrhoini

   tt    =0d0
   istate=1
   to=tign/2d0**3
   Ymax=0d0
   do i=1,7
      to=to*2d0
      call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
      Y=n(1:ns)*MWs(1:ns)
      sm=sum(Y,1)
      Y=Y/sm
      do j=1,ns
         if(Y(j)>Ymax(j)) Ymax(j)=Y(j)
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
subroutine reduction_species(Tini,vrhoini,tign,Teq,allowable_limit)!{{{
   use chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::tign
   double precision,intent(in)::Teq
   double precision,intent(in)::allowable_limit

   integer ns_org

   call sort_n(Tini,vrhoini,tign,.true.)

   ns_org=ns
   do
      call reduct_species_in_order(Tini,vrhoini,tign,Teq,allowable_limit)
      if(ns .eq. ns_org) exit
      ns_org=ns
      call sort_n(Tini,vrhoini,tign,.false.)
   end do
   call reduct_every_species(Tini,vrhoini,tign,Teq,allowable_limit)
end subroutine reduction_species!}}}
subroutine reduct_species_in_order(Tini,vrhoini,tign,Teq,allowable_limit)!{{{
   use chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::tign
   double precision,intent(in)::Teq
   double precision,intent(in)::allowable_limit

   logical flag
   integer nsl,nsu

   nsl=num_reac
   nsu=ns
   do while(nsl.ne.nsu)
      ns=(nsl+nsu)/2
      call change_nr_via_ns
      call check_model(Tini,vrhoini,tign,Teq,allowable_limit,flag)
      if(flag) then
         nsu=ns
      else
         nsl=ns
         if(nsu-nsl .eq. 1) exit
      end if
   end do
   ns=nsu
end subroutine reduct_species_in_order!}}}
subroutine reduct_every_species(Tini,vrhoini,tign,Teq,allowable_limit)!{{{
   use chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::tign
   double precision,intent(in)::Teq
   double precision,intent(in)::allowable_limit

   logical flag
   integer ind

   ind = ns
   ns  = ns-1
   do while(ind .ne. num_reac)
      call swap_species(ind,ns+1)
      call change_nr_via_ns
      call check_model(Tini,vrhoini,tign,Teq,allowable_limit,flag)
      if(flag) then
         ind=ns
         ns =ns-1
         cycle
      end if
      ind=ind-1
   end do
   ns=ns+1
   call change_nr_via_ns
end subroutine reduct_every_species!}}}
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
subroutine sort_r(Tini,vrhoini,tign,flag)!{{{
   use chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::tign
   logical         ,intent(in)::flag

   double precision T
   double precision vrho(ns)
   double precision csp(nr)

   double precision tt,to
   double precision n(ns+1)
   double precision RWORK(LRW)
   integer          IWORK(LIW)
   integer          istate,i,j
   double precision atol

   double precision sm,dtemp

   T   =Tini
   vrho=vrhoini

   tt    =0d0
   istate=1
   to=tign/2d0**3
   csp=0d0
   do i=1,7
      to=to*2d0
      call reaction_cont(T,tt,to,vrho,n,istate,RWORK,IWORK,atol)
      call calc_CSP(n, csp)
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
subroutine reduction_reactions(Tini,vrhoini,tign,Teq,allowable_limit)!{{{
   use chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::tign
   double precision,intent(in)::Teq
   double precision,intent(in)::allowable_limit

   integer nr_org

   call sort_r(Tini,vrhoini,tign,.true.)

   nr_org=nr
   do
      call reduct_reactions_in_order(Tini,vrhoini,tign,Teq,allowable_limit)
      if(nr .eq. nr_org) exit
      nr_org=nr
      call sort_r(Tini,vrhoini,tign,.false.)
   end do
   call reduct_every_reactions(Tini,vrhoini,tign,Teq,allowable_limit)
end subroutine reduction_reactions!}}}
subroutine reduct_reactions_in_order(Tini,vrhoini,tign,Teq,allowable_limit)!{{{
   use chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::tign
   double precision,intent(in)::Teq
   double precision,intent(in)::allowable_limit

   logical flag
   integer nrl,nru

   nrl=1
   nru=nr
   do while(nrl.ne.nru)
      nr=(nrl+nru)/2
      call check_model(Tini,vrhoini,tign,Teq,allowable_limit,flag)
      if(flag) then
         nru=nr
      else
         nrl=nr
         if(nru-nrl .eq. 1) exit
      end if
   end do
   nr=nru
end subroutine reduct_reactions_in_order!}}}
subroutine reduct_every_reactions(Tini,vrhoini,tign,Teq,allowable_limit)!{{{
   use chem
   implicit none
   double precision,intent(in)::Tini
   double precision,intent(in)::vrhoini(ns)
   double precision,intent(in)::tign
   double precision,intent(in)::Teq
   double precision,intent(in)::allowable_limit

   logical flag
   integer ind

   ind = nr
   nr  = nr-1
   do while(ind .ne. 0)
      call swap_reaction(ind,nr+1)
      call check_model(Tini,vrhoini,tign,Teq,allowable_limit,flag)
      if(flag) then
         ind=nr
         nr =nr-1
         cycle
      end if
      ind=ind-1
   end do
   nr=nr+1
end subroutine reduct_every_reactions!}}}

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

!boundary condition
subroutine calc_boundary(p,T,deg, wt,vhi)!{{{
   use grbl_prmtr
   use const_chem
   implicit none
   double precision,intent(in)::p
   double precision,intent(inout)::T
   double precision,intent(in)::deg

   double precision,intent(out)::wt(dimw)
   double precision,intent(out)::vhi(ns)

   double precision,dimension(ns)::vrho,DHi
   double precision rho,E,MWave
   integer j

   call calc_vrho(p,T,deg, rho,vrho,E)
   call calc_T(vrho,rho,E, T, wt(indxg),MWave,DHi,vhi,wt(indxMu))

   !set wt
   wt(1)=rho
   wt(2)=0d0
   wt(3)=0d0
   wt(4)=p
   do j=1,ns
      wt(4+j)=vrho(j)/rho
   end do
   wt(indxht)=E+p/rho
   wt(indxR) =Ru*1d-4/MWave
end subroutine calc_boundary!}}}

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

   integer  neq
   external Fex_var

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
   neq = ns+1

   nloop=tout/dt
   do i=1,nloop
      to=to+dt
      if(i .eq. nloop) to=tout
      call chemeq2(Fex_var,neq,n,to,rtol_user,atol_user,tt)
   end do

   !n to vrho
   do j=1,ns
      n(j)=max(n(j),n_eps)
      vrho(j)=n(j)*MWs(j)*1d3
   end do
end subroutine reaction_chemeq2!}}}
subroutine chemeq2(Fex,neq,y,tout,rtol,atol,tt)!{{{
   use chem
   implicit none
   external Fex
   integer,intent(in)::neq
   double precision,intent(inout) ::y(neq)
   double precision,intent(in)    ::tout
   double precision,intent(in)    ::rtol
   double precision,intent(in)    ::atol
   double precision,dimension(neq)::y0,yp,p,q,p0,q0,alp
   double precision tt,dt,ratio,var,ex
   integer i,counter
   integer imin

   !j=0
   dt=tout-tt
   imin=1
   do
      y0=y
      if(tout-tt<dt) dt=tout-tt
      counter=0
      do
         y=y0
         !predictor
         call Fex(y,q0,p0)
         p0=-p0/(y+1d-300)
         do i=1,neq-1
            var=p0(i)*dt
            if(var < 1d0) then
               alp(i)=(180d0+60d0*var+11d0*var**2+var**3)/(360d0+60d0*var+12d0*var**2+var**3)
            else
               var=1d0/var
               alp(i)=(180d0*var**3+60d0*var**2+11d0*var+1d0)/(360d0*var**3+60d0*var**2+12d0*var+1d0)
            end if
         end do
         do i=1,neq-1
            yp(i)=y0(i)+dt*(q0(i)-p0(i)*y0(i))/(1d0+alp(i)*p0(i)*dt)
            if(isNaN(yp(i))) then
               print *,"flag1"
               print *,i,yp(i),y0(i)
               print *,q0(i),p0(i),alp(i)
               stop
            end if
         end do
         yp(neq) =y0(neq)+dt*(q0(neq)-p0(neq)*y0(neq))
         
         !corrector
         call Fex(yp,q,p)
         p=-p/(yp+1d-300)
         p  =0.5d0*(p0+p)
         do i=1,neq-1
            var=p(i)*dt
            if(var < 1d0) then
               alp(i)=(180d0+60d0*var+11d0*var**2+var**3)/(360d0+60d0*var+12d0*var**2+var**3)
            else
               var=1d0/var
               alp(i)=(180d0*var**3+60d0*var**2+11d0*var+1d0)/(360d0*var**3+60d0*var**2+12d0*var+1d0)
            end if
         end do
         q  =alp*q+(1d0-alp)*q0
         do i=1,neq-1
            y(i) =y0(i)+dt*(q(i) -p(i) *y0(i))/(1d0+alp(i)*p(i) *dt)
            if(isNaN(y(i))) then
               print *,"flag2"
               print *,i,y(i)
               print *,q(i),p(i),p0(i),alp(i)
               stop
            end if
         end do
         y(neq) =y0(neq)+dt*(q(neq)-p(neq)*y0(neq))

         !calc residual
         ratio=1d300
         do i=1,neq-1
            var = (y(i)*rtol+rtol*1d-5)/abs(y(i)-yp(i)+1d-300)
            if(y(i)>y0(i) .and. ratio>var) then
               ratio=var
               imin = i
            end if
            !ratio=min(ratio,(y(i)*rtol+atol)/abs(y(i)-yp(i)+1d-300))
         end do

         !print *,"ratio",ratio
         !j=j+1
         if(1d0<ratio) exit
         !print *,"recalc"
         !dt=dt*sqrt(ratio)*0.8d0
         dt=dt*sqrt(ratio)/dble(counter+1)
         counter=counter+1
         write(40,*) dt,ratio
         if(counter>10) stop "not converted at CHEMEQ2"
      end do
      close(40)
      tt=tt+dt
      write(20,'(2es15.7,3es9.1,2i5,x,a8,es9.1)') tt,y(neq),y(imin),y0(imin),dt,counter,&
                                              imin,trim(SYM_SPC(imin)),alp(imin)
      if(tout <= tt) exit
      dt=dt*sqrt(ratio)
   end do
   !print *,j
end subroutine chemeq2!}}}

