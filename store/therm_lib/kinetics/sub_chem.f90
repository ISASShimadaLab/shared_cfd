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


! data_read routines
subroutine set_reac_and_therm_data!{{{
   use mod_mpi
   use const_chem
   use chem
   !results
   double precision   MWe(ne)
   double precision T_def(3)

   !buffers
   character*100    linebuf
   character*50     line_reactants
   character*50     line_products
   character*100    str
   character*24     cha(2)
   character*2      scha(5)
   character        sscha
   integer          inte(6)
   double precision num
   integer          inum
   double precision T_buf(3)

   !flag etc.
   integer state !0...unormal, 1...elements, 2...species, 3...thermo,4...reactions
   integer ind
   integer i,j,k,l
   logical flag
   logical thermo_flag
   integer thermo_line
   integer elm_found
   integer spc_found

   integer NUM_ELM
   integer NUM_SPC
   integer NUM_RCT

   integer,parameter::LMW =86
   integer,parameter::LCHM=87
   integer,parameter::LOUT=88

   if(myid .eq. 0) then
      state=0
      open(LCHM,file="chem.inp",status="old")

      do
         !header{{{
         read(LCHM,'(a100)',end=95) linebuf
         call rm_comment(linebuf)
         if(len_trim(linebuf) .eq. 0) cycle
         if(state .ne. 3) call chomp(linebuf)
         !if(state .eq. 4) print '(a)',trim(linebuf)
         !}}}

         if("end" .eq. trim(linebuf)) then
            !finalization of each status
            if     (state .eq. 1) then
               NUM_ELM=i
               !read MW.log{{{
               elm_found=0

               open(LMW,file="MW.inp",status="old")
               do!{{{
                  read(LMW,'(a100)',end=91) linebuf
                  call rm_comment(linebuf)
                  if(len_trim(linebuf) .eq. 0) cycle

                  call read_head_char(linebuf,cha(1),flag)
                  do i=1,NUM_ELM
                     if(trim(cha(1)) .eq. trim(SYM_ELM(i))) then
                        elm_found=elm_found+1
                        call read_tail_number(linebuf,MWe(i))
                        exit
                     end if
                  end do
               end do!}}}
               91 close(LMW)

               if( NUM_ELM .ne. elm_found) go to 99
               !}}}
            else if(state .eq. 2) then
               NUM_SPC=i
            else if(state .eq. 3) then
               if(spc_found .ne. NUM_SPC) then!{{{
                  print *,"Error: spc_found:",spc_found," is not equal to NUM_SPC:",NUM_SPC
                  do i=1,NUM_SPC
                     if(Tthre(1,i) .eq. 0d0) then
                        print '(a,3f15.7)',SYM_SPC(i),Tthre(1,i),Tthre(2,i),Tthre(3,i)
                     end if
                  end do
                  stop
               end if!}}}
            end if
            state=0
         else if(state .eq. 0) then
            !data type config
            if     (trim(linebuf) .eq. "elements") then
               state=1
            else if(trim(linebuf) .eq. "species") then
               state=2
            else if(trim(linebuf) .eq. "thermo all") then
               state=3
               spc_found=0
               thermo_flag=.true.
               thermo_line=1
            else if(trim(linebuf) .eq. "reactions      cal/mole  moles") then
               state=4
               NUM_RCT=0
            else
               stop "Error!"
            end if
            i=0
         else if(state .eq. 1) then
            !elements data
            do
               i=i+1
               if(i > ne) stop "chem.inp element number exceeded assigned one."
               call read_head_char(linebuf,cha(1),flag)
               SYM_ELM(i)=cha(1)
               if(flag) exit
            end do
         else if(state .eq. 2) then
            !species data
            do
               i=i+1
               if(i > ns) stop "chem.inp species number exceeded assigned one."
               call read_head_char(linebuf,cha(1),flag)
               SYM_SPC(i)=cha(1)
               if(flag) exit
            end do
         else if(state .eq. 3) then
            !thermo data
            if(thermo_flag) then
               !read default common threshold temperatures{{{
               read(linebuf,'(3F10.0)') T_def(1),T_def(2),T_def(3)
               thermo_flag=.false.
               !}}}
            else
               if(thermo_line .eq. 1) then!{{{
                  !line 1{{{
                  read(linebuf,'(a24,4(a2,i3),a1,2f10.0,f8.0,a2,i3,1x,i1)') &
                     cha(2),(scha(i),inte(i),i=1,4),sscha,T_buf(1),T_buf(2),T_buf(3),&
                     scha(5),inte(5),inte(6)
                  ind=index(cha(2)," ")
                  cha(1)=cha(2)(:ind-1)
                  call search_SYM(SYM_SPC,NUM_SPC,cha(1),k)

                  if(k.gt.0) then!{{{
                     spc_found=spc_found+1
                     do i=1,5!{{{
                        if(len_trim(scha(i)) .eq. 0) cycle
                        call search_SYM(SYM_ELM,NUM_ELM,scha(i),j)
                        ES(j,k)=inte(i)
                     end do!}}}

                     !set MW{{{
                     MWs(k)=dot_product(ES(:NUM_ELM,k),MWe(:NUM_ELM))
                     !}}}
                    
                     if(sscha .ne. "G") then!{{{
                        print *,"Error! Odd phase.",sscha
                        stop
                     end if!}}}

                     !Temperature{{{
                     do i=1,3
                        if(T_buf(i) .ne. 0d0) then
                           Tthre(i,k)=T_buf(i)
                        else
                           Tthre(i,k)=T_def(i)
                        end if
                     end do
                     !}}}
                  end if!}}}
                  thermo_line=thermo_line+1
                  !}}}
               else if(thermo_line .eq. 2) then
                  !line 2{{{
                  thermo_line=thermo_line+1
                  if(k .le. 0) cycle
                  read(linebuf,'(5E15.0)') (coeff(i,2,k),i=1,5)
                  !}}}
               else if(thermo_line .eq. 3) then
                  !line 3{{{
                  thermo_line=thermo_line+1
                  if(k .le. 0) cycle
                  read(linebuf,'(5E15.0)') (coeff(i,2,k),i=6,7),(coeff(i,1,k),i=1,3)
                  !}}}
               else if(thermo_line .eq. 4) then
                  !line 4{{{
                  thermo_line=1
                  if(k .le. 0) cycle
                  read(linebuf,'(4E15.0)') (coeff(i,1,k),i=4,7)
                  !}}}
               else
                  !catch error{{{
                  print *,"Error"
                  !print *,"k=",k
                  !print '(a,3es7.1)',"ES=",(ES(i,k),i=1,NUM_ELM)
                  !print '(a,3es15.7)',"Tthre=",(Tthre(k,i),i=1,3)
                  !print '(a,7es15.7)',"Coeff(1)=",(coeff(i,1,k),i=1,7)
                  !print '(a,7es15.7)',"Coeff(2)=",(coeff(i,2,k),i=1,7)
                  stop
                  !}}}
               end if!}}}
            end if
         else if(state .eq. 4) then
            !reactions
            if(trim(linebuf) .eq. "duplicate" .or. trim(linebuf) .eq. "DUPLICATE") cycle

            ind=index(linebuf,"/")
            if(ind .eq. 0) then
               !main{{{
               NUM_RCT=NUM_RCT+1
               if(NUM_RCT > nr) stop "chem.inp reaction number exceeded assigned one."

               !read default ABE{{{
               call read_tail_number(linebuf,NUM)
               ABE(3,NUM_RCT)=NUM*invRc
               call read_tail_number(linebuf,NUM)
               ABE(2,NUM_RCT)=NUM
               call read_tail_number(linebuf,NUM)
               ABE(1,NUM_RCT)=log(NUM)
               !}}}

               SYM_RCT(NUM_RCT)=linebuf
               ind=index(linebuf,"=")
               line_reactants=linebuf(:ind-1)
               line_products=linebuf(ind+1:)

               if(line_products(1:1) .ne. ">") then!{{{
                  Rstate(NUM_RCT,1)=0
               else
                  Rstate(NUM_RCT,1)=1
                  inum=len_trim(line_products)
                  line_products(:inum-1)=line_products(2:inum)
                  line_products(inum:inum)=""
               end if!}}}
              
               !about M{{{
               ind=index(line_reactants,"(+M)")
               if(ind>0) then
                  if(ind+3 .ne. len_trim(line_reactants)) go to 99
                  line_reactants=line_reactants(:len_trim(line_reactants)-4)
                  line_products =line_products (:len_trim(line_products )-4)
                  exist_M(NUM_RCT)=.true.
                  NumM(   NUM_RCT)=0
                  Rstate(NUM_RCT,2)=0
               else
                  ind=index(line_reactants,"+M")
                  if(ind>0) then
                     if(ind+1 .ne. len_trim(line_reactants)) go to 99
                     line_reactants=line_reactants(:len_trim(line_reactants)-2)
                     line_products =line_products (:len_trim(line_products )-2)
                     exist_M(NUM_RCT)=.true.
                     NumM(   NUM_RCT)=0
                     Rstate(NUM_RCT,2)=4
                  else
                     exist_M(NUM_RCT)=.false.
                     Rstate(NUM_RCT,2)=0
                  end if
               end if
               !}}}

               !reactants{{{
               i=1
               do
                  ind=index(line_reactants,"+")
                  if(ind .eq. 0) then
                     call search_SYM(SYM_SPC,NUM_SPC,line_reactants,j)
                     IndNu(i,1,NUM_RCT)=j
                     exit
                  else
                     call search_SYM(SYM_SPC,NUM_SPC,line_reactants(:ind-1),j)
                     IndNu(i,1,NUM_RCT)=j
                     line_reactants=line_reactants(ind+1:)
                     i=i+1
                  end if
               end do
               NumNu(1,NUM_RCT)=i
               !}}}

               !products{{{
               i=1
               do
                  ind=index(line_products,"+")
                  if(ind .eq. 0) then
                     call search_SYM(SYM_SPC,NUM_SPC,line_products,j)
                     IndNu(i,2,NUM_RCT)=j
                     exit
                  else
                     call search_SYM(SYM_SPC,NUM_SPC,line_products(:ind-1),j)
                     IndNu(i,2,NUM_RCT)=j
                     line_products=line_products(ind+1:)
                     i=i+1
                  end if
               end do
               NumNu(2,NUM_RCT)=i
               !}}}
               !}}}
            else
               !auxiliary{{{
               cha(1) =linebuf(:ind-1)
               linebuf=linebuf(ind:)
               if     (trim(cha(1)) .eq. "REV"       .or. trim(cha(1)) .eq. "rev")  then
                    !rev{{{
                    if(Rstate(NUM_RCT,1) .ne. 0) go to 99

                    select case(Rstate(NUM_RCT,2))
                       case(2,3)
                          go to 99
                       case(0)
                          Rstate(NUM_RCT,2)=1
                       case(4)
                          Rstate(NUM_RCT,2)=5
                       case default
                          print *,"Odd rev! Rstate2=",Rstate(NUM_RCT,2)
                          go to 99
                    end select

                    call fetch_slash(linebuf,str)
                    call read_tail_number(str,NUM)
                    cABE(3,NUM_RCT)=NUM*invRc
                    call read_tail_number(str,NUM)
                    cABE(2,NUM_RCT)=NUM
                    call read_tail_number(str,NUM)
                    cABE(1,NUM_RCT)=log(NUM)

                    if(len_trim(str) .ne. 0) go to 99
                    !}}}
               else if(trim(cha(1)) .eq. "LOW"       .or. trim(cha(1)) .eq. "low")  then
                    !low{{{
                    if(Rstate(NUM_RCT,1) .ne. 0) go to 99
                    select case(Rstate(NUM_RCT,2))
                       case(1,3)
                          go to 99
                       case(0)
                          Rstate(NUM_RCT,2)=2
                       case default
                          print *,"Odd low! Rstate2=",Rstate(NUM_RCT,2)
                          go to 99
                    end select

                    if(Rstate(NUM_RCT,2) .eq. 1 .or. Rstate(NUM_RCT,2) .eq. 3) go to 99
                    Rstate(NUM_RCT,2)=2

                    call fetch_slash(linebuf,str)
                    call read_tail_number(str,NUM)
                    cABE(3,NUM_RCT)=NUM*invRc
                    call read_tail_number(str,NUM)
                    cABE(2,NUM_RCT)=NUM
                    call read_tail_number(str,NUM)
                    cABE(1,NUM_RCT)=log(NUM)

                    if(len_trim(str) .ne. 0) go to 99
                    !}}}
               else if(trim(cha(1)) .eq. "TROE"      .or. trim(cha(1)) .eq. "troe") then
                    !troe{{{
                    if(Rstate(NUM_RCT,1) .ne. 0) go to 99
                    if(Rstate(NUM_RCT,2) .ne. 2) go to 99
                    Rstate(NUM_RCT,2)=3

                    call fetch_slash(linebuf,str)
                    call read_tail_number(str,NUM)
                    TROE(4,NUM_RCT)=NUM
                    call read_tail_number(str,NUM)
                    TROE(3,NUM_RCT)=1d0/NUM
                    call read_tail_number(str,NUM)
                    TROE(2,NUM_RCT)=1d0/NUM
                    call read_tail_number(str,NUM)
                    TROE(1,NUM_RCT)=NUM

                    if(len_trim(str) .ne. 0) go to 99
                    !}}}
               else
                    !M enhanced{{{
                    if(.not. exist_M(NUM_RCT)) go to 99
                    Men(:,NUM_RCT)=0d0

                    do j=1,NUM_SPC
                       call search_SYM(SYM_SPC,NUM_SPC,cha(1),IndM(j,NUM_RCT))
                       call fetch_slash(linebuf,str)

                       call read_tail_number(str,Men(j,NUM_RCT))
                       Men(j,NUM_RCT)=Men(j,NUM_RCT)-1d0
                       if(len_trim(linebuf) .eq. 0) exit
                       ind=index(linebuf,"/")
                       cha(1) =linebuf(:ind-1)
                       linebuf=linebuf(ind:)
                    end do
                    NumM(   NUM_RCT)=j
                    !}}}
               end if
               !}}}
            end if
         else
            stop "Error."
         end if
      end do
95    close(LCHM)

      !check consistency of the numbers
      if(NUM_ELM .ne. ne) then
         print *,"The Number of Elements is different between input file and hard-code: ",NUM_ELM,ne," each other."
         stop
      end if
      if(NUM_SPC .ne. ns) then
         print *,"The Number of Specices is different between input file and hard-code: ",NUM_SPC,ns," each other."
         stop
      end if
      if(NUM_RCT .ne. nr .and. .not. (nr .eq. 1 .and. NUM_RCT .eq. 0)) then
         print *,"The Number of Reactions is different between input file and hard-code: ",NUM_RCT,nr," each other."
         stop
      end if

      !!!!!!!!!!! post process !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,NUM_RCT
         snu(i)=NumNu(2,i)-NumNu(1,i)
      end do
      do j=1,NUM_SPC
         invMW(j)=1d0/MWs(j)
      end do
      !!!!!!!!!!! post process !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end if

   call MPI_Bcast(NUM_ELM,      1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(NUM_SPC,      1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(NUM_RCT,      1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
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

   return

99 print *,"Error"
   stop
contains
subroutine rm_comment(line)!{{{
   character(*),intent(inout)::line
   integer ind

   ind=index(line,"!")
   if(ind .gt. 1) then
      line=line(:ind-1)
      line=trim(line)
   else if(ind .eq. 1) then
      line=""
   else
      line=trim(line)
   end if
end subroutine rm_comment!}}}

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

subroutine read_tail_char(line,cha,flag)!{{{
   character(*),intent(inout)::line
   character(*),intent(out)::cha
   logical,intent(out)::flag

   integer ind

   flag=.false.


   ind=index(trim(line)," ",.true.)
   if(ind .ne. 0) then
      cha=line(ind+1:)
      line=trim(line(:ind))
   else
      cha=line
      flag=.true.
   end if
end subroutine read_tail_char!}}}

subroutine read_head_char(line,cha,flag)!{{{
   character(*),intent(inout)::line
   character(*),intent(out)::cha
   logical,intent(out)::flag

   integer ind

   flag=.false.

   ind=index(trim(line)," ")
   if(ind .ne. 0) then
      cha=line(:ind-1)
      line=trim(adjustl(line(ind:)))
   else
      cha=line
      flag=.true.
   end if
end subroutine read_head_char!}}}

subroutine search_SYM(SYM,NUM_SYM,str,i)!{{{ 
   character(*),dimension(*),intent(in)::SYM
   integer,intent(in)                  ::NUM_SYM
   character(*),intent(in)             ::str
   integer,intent(out)                 ::i

   i=1
   do while(i<=NUM_SYM)
      !print *,i,":",trim(SYM(i)),":",trim(str)
      if(trim(SYM(i)) .eq. trim(str)) exit
      i=i+1
   end do

   if(i>NUM_SYM) i=-1
end subroutine search_SYM!}}}

subroutine increment_RS(SYM_SPC,NUM_SPC,NUM_RCT,str,RS)!{{{
   use const_chem
   character(*),dimension(*),intent(in)   ::SYM_SPC
   integer                  ,intent(in)   ::NUM_SPC
   integer                  ,intent(in)   ::NUM_RCT
   character(*)             ,intent(in)   ::str
   integer                  ,intent(inout)::RS(nr,*)

   integer j

   call search_SYM(SYM_SPC,NUM_SPC,str,j)
   RS(NUM_RCT,j)=RS(NUM_RCT,j)+1
end subroutine increment_RS!}}}

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
end subroutine set_reac_and_therm_data!}}}
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
            !print *,nt,trim(sname)
            num_sctn_trans(nt)    =iv
            tr2th(nt)             =itr2th
            !print *,nt,sname

            !for V
            do i=1,iv
               read(10,'(1x,a1,2(f7.1,2x),4e15.8)') v,(Trange_trans(j,i,nt),j=1,2),(trans(j,i,nt),j=1,4)
               if(v .ne. "V") then
                  print *,"ERROR : ODD FORMAT AT trans.inp"
                  call exit(1)
               end if
            end do

            !for C --- which are not used now.
            do i=1,ic
               read(10,'(a)') comment_long
               if(comment_long(2:2) .ne. "C") then
                  print *,"ERROR : ODD FORMAT AT trans.inp"
                  call exit(1)
               end if
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

! reaction
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
      n(j) = max(n(j),1d-20)
   end do

   n(ns+1)= T
   neq         = ns+1
   tt=0d0

   call Fex(neq, tt, n, dndt, rpar, ipar)

   dTdt=dndt(ns+1)
end subroutine calc_dTdt!}}}

