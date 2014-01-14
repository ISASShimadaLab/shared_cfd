! read data files
subroutine read_cheminp!{{{
   use mod_mpi
   use chem
   implicit none
   character*200 bufnext,line
   integer i

   if(myid .eq. 0) then
      open(22,file="chem.inp")
      bufnext="" !initialize bufnext

      !elements
      line=next_line()
      i=0
      do
         line=next_word()
         if(trim(line) .eq. "end") exit
         i=i+1
         SYM_ELM(i)=trim(line)
      end do

      !species
      line=next_line()
      ns=0
      do
         line=next_word()
         if(trim(line) .eq. "end") exit
         ns=ns+1
         SYM_SPC(ns)=trim(line)
      end do
      print *,"the number of species of chem.inp:",ns

      close(22)
   end if
   call MPI_Bcast(      ns,       1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast( SYM_SPC,   ns*18,        MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast( SYM_ELM,    ne*2,        MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
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
   integer ind
   do while(bufnext .eq. "")
      read(22,'(a)') bufnext
      bufnext=trim(adjustl(bufnext))
      ind = index(bufnext,'!')
      if(ind == 1) then
         bufnext=""
      else if(ind > 1) then
         bufnext=bufnext(:ind-1)
      end if
   end do
   next_line=bufnext
   bufnext=""
end function next_line
end subroutine read_cheminp!}}}
subroutine set_therm_data!{{{
   use mod_mpi
   use chem
   implicit none
   character*18     sname
   character*300    comment
   character*2      name_elm(5)
   double precision Nelm(5)
   integer,external::search_species
   integer,external::search_elements

   integer i,j,k,l,Nfound
   double precision tmp,MWe(ne)

   if(myid .eq. 0) then
      open(10,file='MW.inp',status='old')
      Nfound=0
      do
         !first line --- name, comment
         read(10,'(A2,f10.5)',end=01) sname, tmp
         j = search_elements(sname)
         if(j>0) then
            MWe(j)=tmp
            Nfound = Nfound+1
            if(Nfound==ne) exit
         end if
      end do
01    close(10)
      if(Nfound .ne. ne) stop "Some elements have not found in MW.inp."

      open(10,file='thermo.inp',status='old')
      do
         read(10,'(A)') comment
         if(comment(1:1) .ne. '!') exit
      end do
      read(10,'()')!default ranges

      Nfound=0
      do
         !first line --- name, comment
         read(10,'(A18,A62)') sname, comment
         j = search_species(sname)
         if(sname(1:3) .eq. 'END') then
            exit
         else if(j<1)then
            read(10,'(i2)') j
            do l=1, j*3
               read(10,'()')
            end do
            cycle
         end if
         Nfound = Nfound+1

         !second line --- other property
         !read other informations
         read(10,'(i2,8x,5(a2,f6.2),2x,f13.7)') &
            num_sctn(j),(name_elm(k), Nelm(k), k=1,5), mw(j)

         !set number of elements
         Ac(1:ne,j)=0d0
         outer:do k=1,5 
            if(name_elm(k) .eq. "  ") exit
            do i=1,ne
               if(SYM_ELM(i) .eq. name_elm(k)) then
                  Ac(i,j)=Nelm(k)
                  cycle outer
               end if
            end do
         end do outer

         !other line ---coefficients of thermodynamical functions
         do l=1, num_sctn(j)
            read(10,'(2(1x ,f10.3))')    Trange(:,l,j)
            read(10,'(5D16.9)')            co(1:5,l,j)
            read(10,'(2D16.8,16x,2D16.8)') co(6:9,l,j)
         end do
         if(Nfound==ns) exit
      end do
      close(10)
      if(Nfound .ne. ns) stop "Some species have not found in thermo.inp."
   end if

   call MPI_Bcast(     MWe,      ne, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(num_sctn,      ns,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(      MW,      ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  Trange,  ns*6*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(      co,  ns*6*9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(      Ac,   ne*ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   ! calculate mass fraction of elements in each species.
   do j=1,ns
      tmp=0d0
      do i=1,ne
         YAc(i,j)=MWe(i)*Ac(i,j)
         tmp=tmp+YAc(i,j)
      end do
      tmp=1d0/tmp
      do i=1,ne
         YAc(i,j)=YAc(i,j)*tmp
      end do
   end do
end subroutine set_therm_data!}}}
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
   integer,external::search_species


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

      print *,"the number of species of trans.inp:",nt

      close(10)
   end if

   call MPI_Bcast(            nt,      1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(num_sctn_trans,     ns,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(         trans, 4*3*ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  Trange_trans, 2*3*ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(         tr2th,     ns,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end subroutine set_trans_data!}}}
subroutine read_fo_composition!{{{
   use mod_mpi
   use chem_var
   implicit none
   character*100    buf
   double precision amount
   integer          ind
   integer,external::search_species

   if(myid .eq. 0) then
      open(8,file="control_chem.inp")
      read(8,'()')
      read(8,'(25x,es15.7)') po
      read(8,'(25x,es15.7)') To
      read(8,'()')
      no = 0d0
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit
         buf=adjustl(trim(buf));ind=index(buf,' ')
         if(ind .eq. 0) stop "Bad format at control_chem.inp while reading Oxidizer Compositions"
         read(buf(ind+1:),*) amount
         no(search_species(buf(1:ind-1)))=amount
      end do

      read(8,'(25x,es15.7)') pf
      read(8,'(25x,es15.7)') Tf
      read(8,'()')
      nf = 0d0
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit
         buf=adjustl(trim(buf));ind=index(buf,' ')
         if(ind .eq. 0) stop "Bad format at control_chem.inp while reading Fuel Compositions"
         read(buf(ind+1:),*) amount
         nf(search_species(buf(1:ind-1)))=amount
      end do

      close(8)
   end if

   !!! MPI COMMUNICATIONS
   call MPI_Bcast(po,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(To,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(no, ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   call MPI_Bcast(pf,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(Tf,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(nf, ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end subroutine read_fo_composition!}}}
subroutine read_p_composition!{{{
   use mod_mpi
   use chem_var
   implicit none
   character*100    buf
   double precision amount
   integer          ind
   integer,external::search_species

   if(myid .eq. 0) then
      open(8,file="control_chem.inp")
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit ! oxidizer
      end do
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit ! fuel
      end do
      read(8,'()')
      np = 0d0
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit
         buf=adjustl(trim(buf));ind=index(buf,' ')
         if(ind .eq. 0) stop "Bad format at control_chem.inp while reading Products Compositions"
         read(buf(ind+1:),*) amount
         np(search_species(buf(1:ind-1)))=amount
      end do
      close(8)
   end if

   !!! MPI COMMUNICATIONS
   call MPI_Bcast(np, ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
end subroutine read_p_composition!}}}

! utilities
integer function search_elements(str)!{{{
   use chem
   character(*),intent(in)::str
   integer i

   do i=1,ne
      if(trim(SYM_ELM(i)) .eq. trim(str)) then
         search_elements=i
         return
      end if
   end do
   search_elements=-1
end function search_elements!}}}
integer function search_species(str)!{{{
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
end function search_species!}}}
integer function search_species_trans(str)!{{{
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
end function search_species_trans!}}}

!cold flow
subroutine init_pack_cold_flow!{{{
   call read_cheminp
   call set_therm_data
   call set_trans_data
   call read_fo_composition

   call initialize_cold_flow
end subroutine init_pack_cold_flow!}}}
subroutine initialize_cold_flow!{{{
   use chem
   use chem_var
   implicit none
   double precision sumo,sumf,sm
   double precision Y(2),DHi(2)
   integer i,j

   !set o/f ratio
   sumo=0d0
   sumf=0d0
   do j=1,ns
      sumo=sumo+MW(j)*no(j)
      sumf=sumf+MW(j)*nf(j)
   end do
   of = sumo/sumf

   !calc n,E
   Y(1)=1d0;Y(2)=0d0
   call calc_ini(pf,Tf, nf, Ef)
   nfini=nf
   call cold_flow(Y,Ef,Tf, MWf,kappaf,muf,DHi,vhif)
   rhof=pf/(Ru*1d3/MWf*Tf)
   Hf  =Ef+pf/rhof
   call set_static_qw(pf,rhof,Tf,Ef,kappaf,muf,Y,qf,wf)

   Y(1)=0d0;Y(2)=1d0
   call calc_ini(po,To, no, Eo)
   noini=no
   call cold_flow(Y,Eo,To, MWo,kappao,muo,DHi,vhio)
   rhoo=po/(Ru*1d3/MWo*To)
   Ho  =Eo+po/rhoo
   call set_static_qw(po,rhoo,To,Eo,kappao,muo,Y,qo,wo)
contains
   subroutine calc_ini(p,T, n, E)
      use func_therm
      implicit none
      double precision,intent(in) ::p
      double precision,intent(in) ::T
      double precision,intent(inout)::n(ns)
      double precision,intent(out)::E
   
      integer j,k,nsc
      double precision vThrt(9)
      double precision sn,sm,tmp
      double precision MWave
   
      sn=0d0
      sm=0d0
      do j=1,ns
         sn=sn+n(j)
         sm=sm+n(j)*MW(j)
      end do
      MWave=sm/sn
      n=n/sm*1d3 !mole/kg
   
      call calc_vThrt(T,log(T),vThrt)
      do j=1,ns
         call check_section_number(j,T,nsc)
         tmp=-1d0
         do k=1,9
            tmp=tmp+vThrt(k)*co(k,nsc,j)
         end do
         E=E+tmp*n(j)
      end do
      E=E*Ru*T
   end subroutine calc_ini
   subroutine set_static_qw(p,rho,T,E,kappa,mu,Y,q_local,w_local)
      implicit none
      double precision,intent(in)::p
      double precision,intent(in)::rho
      double precision,intent(in)::T
      double precision,intent(in)::E
      double precision,intent(in)::kappa
      double precision,intent(in)::mu
      double precision,intent(in)::Y(nY)
      double precision,intent(out)::q_local(dimq)
      double precision,intent(out)::w_local(dimw)
      integer i
      do i=1,nY
         q_local(i)   = Y(i)*rho
         w_local(i+4) = Y(i)
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
end subroutine initialize_cold_flow!}}}
subroutine cold_flow(Y,E, T, MWave,kappa,mu,DHi,vhi)!{{{
   use chem
   use chem_var
   use func_therm
   implicit none
   double precision,intent(in)   ::Y(2)
   double precision,intent(in)   ::E
   double precision,intent(inout)::T
   double precision,intent(out)  ::MWave
   double precision,intent(out)  ::kappa
   double precision,intent(out)  ::mu
   double precision,intent(out)  ::DHi(2)
   double precision,intent(out)  ::vhi(2)

   double precision n(  ns)
   double precision muN(ns)

   double precision T_old

   integer secN,sect
   double precision denom
   double precision cprnow,cvrnow,ERTnow,logT,sn,phi
   double precision,dimension(9)::vThrt,vTcpr
   double precision,dimension(4)::vTmu
   double precision MWi,MWj
   double precision,dimension(ns)::vhi_s,DHi_s

   double precision tmp,tmp2
   integer i,j,k

   !set n
   n=Y(1)*nf+Y(2)*no !mole/kg

   sn=sum(n,1)
   MWave=1d3/sn !(g/kg)/(mole/kg)=g/mole

   !calc T
   logT=log(T)
   ERTnow=0d0;cvrnow=0d0
   call calc_vThrt( T,logT,vThrt)
   call calc_vTcpr( T,     vTcpr)
   do i=1,ns
      call check_section_number(i,T,secN)
      tmp =-1;tmp2=-1
      do j=1,9
         tmp  = tmp +co(j,secN,i)*vThrt(j)
         tmp2 = tmp2+co(j,secN,i)*vTcpr(j)
      end do
      ERTnow = ERTnow + tmp *n(i)
      cvrnow = cvrnow + tmp2*n(i)
   end do

   k=0
   T_old = T*0.5d0
   do while(abs(E/Ru-ERTnow*T)>cvrnow*eps*1d-2 .and. abs(T-T_old)>eps .and. k<100)
      T_old =T
      T     =T+omega*(E/Ru-ERTnow*T)/cvrnow
      logT=log(T)

      ERTnow=0d0;cvrnow=0d0
      call calc_vThrt( T,logT,vThrt)
      call calc_vTcpr( T,     vTcpr)
      do i=1,ns
         call check_section_number(i,T,secN)
         tmp =-1;tmp2=-1
         do j=1,9
            tmp  = tmp +co(j,secN,i)*vThrt(j)
            tmp2 = tmp2+co(j,secN,i)*vTcpr(j)
         end do
         ERTnow = ERTnow + tmp *n(i)
         cvrnow = cvrnow + tmp2*n(i)
      end do

      k=k+1
   end do

   if(k >= 100) then
      print *, "Not Converted at cold flow"
      !call exit(1)
   end if

   !calc kappa
   cprnow= cvrnow+sn
   kappa=1d0+sn/cvrnow

   !calc DHi and vhi
   do i=1,ns
      call check_section_number(i,T,secN)
      tmp= 0d0
      do j=1,9
         tmp =tmp +co(j,secN,i)*vThrt(j)
      end do
      DHi_s(i)=1d0*(kappa-(kappa-1d0)*tmp)
      vhi_s(i)=1d0*tmp
   end do
   tmp=Ru*T
   DHi(1)=dot_product(nf,DHi_s)*tmp
   DHi(2)=dot_product(no,DHi_s)*tmp
   vhi(1)=dot_product(nf,vhi_s)*tmp
   vhi(2)=dot_product(no,vhi_s)*tmp

   !calc mu
   call calc_vTmu(T,logT,vTmu)
   do i=1,nt
      call check_section_number_trans(i,T,sect)
      muN(i)=exp(dot_product(trans(:,sect,i),vTmu))
   end do

   mu = 0d0
   do i=1,nt
      denom = 0d0
      MWi=MW(tr2th(i))
      do j=1,nt
         MWj=MW(tr2th(j))
         phi = 0.25d0*(1d0+sqrt(muN(i)/muN(j))*(MWj/MWi)**0.25d0 )**2 &
                     *sqrt(2d0*MWj/(MWi+MWj))
         denom = denom + n(tr2th(j)) * phi
      end do
      mu = mu + n(tr2th(i)) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine cold_flow!}}}

!flame sheet
subroutine init_pack_flame_sheet!{{{
   call read_cheminp
   call set_therm_data
   call set_trans_data
   call read_fo_composition
   call read_p_composition

   call initialize_flame_sheet
end subroutine init_pack_flame_sheet!}}}
subroutine initialize_flame_sheet!{{{
   use chem
   use chem_var
   implicit none
   double precision sumo,sumf,sm
   double precision Y(2),DHi(2)
   double precision,dimension(ne)::ner,nep
   integer i,j
   !check balance
   ner=0d0
   nep=0d0

   ! calc the amount of each element
   do j=1,ns
      do i=1,ne
         ner(i)=ner(i)+Ac(i,j)*(no(j)+nf(j))
         nep(i)=nep(i)+Ac(i,j)* np(j)
      end do
   end do

   ! check consistency
   do i=1,ne
      if(abs((ner(i)-nep(i))/(ner(i)+nep(i)))>1d-3) then
         print *,"ERROR: The amount of element: ",SYM_ELM(i)
         print *,"ERROR: does not balance between reactants and products."
         print *,"ERROR: Please modify 'control_chem.inp'."
         call exit(1)
      end if
   end do

   !set o/f ratio
   sumo=0d0
   sumf=0d0
   do j=1,ns
      sumo=sumo+MW(j)*no(j)
      sumf=sumf+MW(j)*nf(j)
   end do
   of = sumo/sumf

   !calc n,E
   Y(1)=1d0;Y(2)=0d0
   call calc_ini(pf,Tf, nf, Ef)
   nfini=nf
   call flame_sheet(Y,Ef,Tf, MWf,kappaf,muf,DHi,Yvf,vhif)
   rhof=pf/(Ru*1d3/MWf*Tf)
   Hf  =Ef+pf/rhof
   call set_static_qw(pf,rhof,Tf,Ef,kappaf,muf,Y,qf,wf)

   Y(1)=0d0;Y(2)=1d0
   call calc_ini(po,To, no, Eo)
   noini=no
   call flame_sheet(Y,Eo,To, MWo,kappao,muo,DHi,Yvo,vhio)
   rhoo=po/(Ru*1d3/MWo*To)
   Ho  =Eo+po/rhoo
   call set_static_qw(po,rhoo,To,Eo,kappao,muo,Y,qo,wo)
 
   sm=0d0
   do j=1,ns
      sm=sm+np(j)*MW(j)
   end do
   np=np/sm*1d3 !mole/kg
contains
   subroutine calc_ini(p,T, n, E)
      use func_therm
      implicit none
      double precision,intent(in) ::p
      double precision,intent(in) ::T
      double precision,intent(inout)::n(ns)
      double precision,intent(out)::E
   
      integer j,k,nsc
      double precision vThrt(9)
      double precision sn,sm,tmp
      double precision MWave
   
      sn=0d0
      sm=0d0
      do j=1,ns
         sn=sn+n(j)
         sm=sm+n(j)*MW(j)
      end do
      MWave=sm/sn
      n=n/sm*1d3 !mole/kg
   
      call calc_vThrt(T,log(T),vThrt)
      do j=1,ns
         call check_section_number(j,T,nsc)
         tmp=-1d0
         do k=1,9
            tmp=tmp+vThrt(k)*co(k,nsc,j)
         end do
         E=E+tmp*n(j)
      end do
      E=E*Ru*T
   end subroutine calc_ini
   subroutine set_static_qw(p,rho,T,E,kappa,mu,Y,q_local,w_local)
      implicit none
      double precision,intent(in)::p
      double precision,intent(in)::rho
      double precision,intent(in)::T
      double precision,intent(in)::E
      double precision,intent(in)::kappa
      double precision,intent(in)::mu
      double precision,intent(in)::Y(nY)
      double precision,intent(out)::q_local(dimq)
      double precision,intent(out)::w_local(dimw)
      integer i
      do i=1,nY
         q_local(i)   = Y(i)*rho
         w_local(i+4) = Y(i)
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
end subroutine initialize_flame_sheet!}}}
subroutine flame_sheet(Y,E, T, MWave,kappa,mu,DHi,Yv,vhi)!{{{
   use chem
   use chem_var
   use func_therm
   implicit none
   double precision,intent(in)   ::Y(2)
   double precision,intent(in)   ::E
   double precision,intent(inout)::T
   double precision,intent(out)  ::MWave
   double precision,intent(out)  ::kappa
   double precision,intent(out)  ::mu
   double precision,intent(out)  ::DHi(2)
   double precision,intent(out)  ::Yv(3)
   double precision,intent(out)  ::vhi(3)

   double precision Yf,Yo,Yp
   double precision n(  ns)
   double precision muN(ns)

   double precision T_old

   integer secN,sect
   double precision denom
   double precision cprnow,cvrnow,ERTnow,logT,sn,phi
   double precision,dimension(9)::vThrt,vTcpr
   double precision,dimension(4)::vTmu
   double precision MWi,MWj
   double precision,dimension(ns)::vhi_s,DHi_s

   double precision tmp,tmp2
   integer i,j,k

   !set n
   Yo = Y(2)-of*Y(1)
   if(Yo>0d0) then
      Yf=0d0
      Yp=Y(1)*(1d0+of)
   else
      Yf=-Yo/of
      Yo=0d0
      Yp=Y(2)*(1d0+1d0/of)
   end if
   n=Yo*no+Yf*nf+Yp*np !mole/kg
   Yv(1)=Yf
   Yv(2)=Yo
   Yv(3)=Yp

   sn=sum(n,1)
   MWave=1d3/sn !(g/kg)/(mole/kg)=g/mole

   !calc T
   logT=log(T)
   ERTnow=0d0;cvrnow=0d0
   call calc_vThrt( T,logT,vThrt)
   call calc_vTcpr( T,     vTcpr)
   do i=1,ns
      call check_section_number(i,T,secN)
      tmp =-1;tmp2=-1
      do j=1,9
         tmp  = tmp +co(j,secN,i)*vThrt(j)
         tmp2 = tmp2+co(j,secN,i)*vTcpr(j)
      end do
      ERTnow = ERTnow + tmp *n(i)
      cvrnow = cvrnow + tmp2*n(i)
   end do

   k=0
   T_old = T*0.5d0
   do while(abs(E/Ru-ERTnow*T)>cvrnow*eps*1d-2 .and. abs(T-T_old)>eps .and. k<100)
      T_old =T
      T     =T+omega*(E/Ru-ERTnow*T)/cvrnow
      logT=log(T)

      ERTnow=0d0;cvrnow=0d0
      call calc_vThrt( T,logT,vThrt)
      call calc_vTcpr( T,     vTcpr)
      do i=1,ns
         call check_section_number(i,T,secN)
         tmp =-1;tmp2=-1
         do j=1,9
            tmp  = tmp +co(j,secN,i)*vThrt(j)
            tmp2 = tmp2+co(j,secN,i)*vTcpr(j)
         end do
         ERTnow = ERTnow + tmp *n(i)
         cvrnow = cvrnow + tmp2*n(i)
      end do

      k=k+1
   end do

   if(k >= 100) then
      print *, "Not Converted at flame_sheet"
      !call exit(1)
   end if

   !calc kappa
   cprnow= cvrnow+sn
   kappa=1d0+sn/cvrnow

   !calc DHi and vhi
   do i=1,ns
      call check_section_number(i,T,secN)
      tmp= 0d0
      do j=1,9
         tmp =tmp +co(j,secN,i)*vThrt(j)
      end do
      DHi_s(i)=(kappa-(kappa-1d0)*tmp)
      vhi_s(i)=tmp
   end do
   tmp=Ru*T
   DHi(1)=dot_product(nf,DHi_s)*tmp
   DHi(2)=dot_product(no,DHi_s)*tmp

   vhi(1)=dot_product(nf,vhi_s)*tmp
   vhi(2)=dot_product(no,vhi_s)*tmp
   vhi(3)=dot_product(np,vhi_s)*tmp

   !calc mu
   call calc_vTmu(T,logT,vTmu)
   do i=1,nt
      call check_section_number_trans(i,T,sect)
      muN(i)=exp(dot_product(trans(:,sect,i),vTmu))
   end do

   mu = 0d0
   do i=1,nt
      denom = 0d0
      MWi=MW(tr2th(i))
      do j=1,nt
         MWj=MW(tr2th(j))
         phi = 0.25d0*(1d0+sqrt(muN(i)/muN(j))*(MWj/MWi)**0.25d0 )**2 &
                     *sqrt(2d0*MWj/(MWi+MWj))
         denom = denom + n(tr2th(j)) * phi
      end do
      mu = mu + n(tr2th(i)) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine flame_sheet!}}}
subroutine flame_sheet_hp(Y,H, T, MWave,kappa,mu,DHi,Yv,vhi)!{{{
   use chem
   use chem_var
   use func_therm
   implicit none
   double precision,intent(in)   ::Y(2)
   double precision,intent(in)   ::H
   double precision,intent(inout)::T
   double precision,intent(out)  ::MWave
   double precision,intent(out)  ::kappa
   double precision,intent(out)  ::mu
   double precision,intent(out)  ::DHi(2)
   double precision,intent(out)  ::Yv(3)
   double precision,intent(out)  ::vhi(3)

   double precision Yf,Yo,Yp
   double precision n(  ns)
   double precision muN(ns)

   double precision T_old

   integer secN,sect
   double precision denom
   double precision cprnow,HRTnow,logT,sn,phi
   double precision,dimension(9)::vThrt,vTcpr
   double precision,dimension(4)::vTmu
   double precision MWi,MWj
   double precision,dimension(ns)::vhi_s,DHi_s

   double precision tmp,tmp2
   integer i,j,k

   !set n
   Yo = Y(2)-of*Y(1)
   if(Yo>0d0) then
      Yf=0d0
      Yp=Y(1)*(1d0+of)
   else
      Yf=-Yo/of
      Yo=0d0
      Yp=Y(2)*(1d0+1d0/of)
   end if
   n=Yo*no+Yf*nf+Yp*np !mole/kg
   Yv(1)=Yf
   Yv(2)=Yo
   Yv(3)=Yp

   sn=sum(n,1)
   MWave=1d3/sn !(g/kg)/(mole/kg)=g/mole

   !calc T
   logT=log(T)
   HRTnow=0d0;cprnow=0d0
   call calc_vThrt( T,logT,vThrt)
   call calc_vTcpr( T,     vTcpr)
   do i=1,ns
      call check_section_number(i,T,secN)
      tmp =0d0;tmp2=0d0
      do j=1,9
         tmp  = tmp +co(j,secN,i)*vThrt(j)
         tmp2 = tmp2+co(j,secN,i)*vTcpr(j)
      end do
      HRTnow = HRTnow + tmp *n(i)
      cprnow = cprnow + tmp2*n(i)
   end do

   k=0
   T_old = T*0.5d0
   do while(abs(H/Ru-HRTnow*T)>cprnow*eps*1d-2 .and. abs(T-T_old)>eps .and. k<100)
      T_old =T
      T     =T+omega*(H/Ru-HRTnow*T)/cprnow
      logT=log(T)

      HRTnow=0d0;cprnow=0d0
      call calc_vThrt( T,logT,vThrt)
      call calc_vTcpr( T,     vTcpr)
      do i=1,ns
         call check_section_number(i,T,secN)
         tmp =0d0;tmp2=0d0
         do j=1,9
            tmp  = tmp +co(j,secN,i)*vThrt(j)
            tmp2 = tmp2+co(j,secN,i)*vTcpr(j)
         end do
         HRTnow = HRTnow + tmp *n(i)
         cprnow = cprnow + tmp2*n(i)
      end do

      k=k+1
   end do

   if(k >= 100) then
      print *, "Not Converted at flame_sheet"
      !call exit(1)
   end if

   !calc kappa
   kappa=cprnow/(cprnow-sn)

   !calc DHi and vhi
   do i=1,ns
      call check_section_number(i,T,secN)
      tmp= 0d0
      do j=1,9
         tmp =tmp +co(j,secN,i)*vThrt(j)
      end do
      DHi_s(i)=(kappa-(kappa-1d0)*tmp)
      vhi_s(i)=tmp
   end do
   tmp=Ru*T
   DHi(1)=dot_product(nf,DHi_s)*tmp
   DHi(2)=dot_product(no,DHi_s)*tmp

   vhi(1)=dot_product(nf,vhi_s)*tmp
   vhi(2)=dot_product(no,vhi_s)*tmp
   vhi(3)=dot_product(np,vhi_s)*tmp

   !calc mu
   call calc_vTmu(T,logT,vTmu)
   do i=1,nt
      call check_section_number_trans(i,T,sect)
      muN(i)=exp(dot_product(trans(:,sect,i),vTmu))
   end do

   mu = 0d0
   do i=1,nt
      denom = 0d0
      MWi=MW(tr2th(i))
      do j=1,nt
         MWj=MW(tr2th(j))
         phi = 0.25d0*(1d0+sqrt(muN(i)/muN(j))*(MWj/MWi)**0.25d0 )**2 &
                     *sqrt(2d0*MWj/(MWi+MWj))
         denom = denom + n(tr2th(j)) * phi
      end do
      mu = mu + n(tr2th(i)) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine flame_sheet_hp!}}}
subroutine flame_sheet_PTgiven(Y,T, E,MWave,kappa,mu,Yv,vhi)!{{{
   use chem
   use chem_var
   use func_therm
   implicit none
   double precision,intent(in) ::Y(2)
   double precision,intent(in) ::T
   double precision,intent(out)::E
   double precision,intent(out)::MWave
   double precision,intent(out)::kappa
   double precision,intent(out)::mu
   double precision,intent(out)::Yv(3)
   double precision,intent(out)::vhi(3)

   double precision Yf,Yo,Yp
   double precision n(  ns)
   double precision muN(ns)

   integer secN,sect
   double precision denom
   double precision cprnow,cvrnow,ERTnow,logT,sn,phi
   double precision,dimension(9)::vThrt,vTcpr
   double precision,dimension(4)::vTmu
   double precision MWi,MWj
   double precision,dimension(ns)::vhi_s

   double precision tmp,tmp2
   integer i,j,k

   !set n
   Yo = Y(2)-of*Y(1)
   if(Yo>0d0) then
      Yf=0d0
      Yp=Y(1)*(1d0+of)
   else
      Yf=-Yo/of
      Yo=0d0
      Yp=Y(2)*(1d0+1d0/of)
   end if
   n=Yo*no+Yf*nf+Yp*np !mole/kg
   Yv(1)=Yf
   Yv(2)=Yo
   Yv(3)=Yp

   sn=sum(n,1)
   MWave=1d3/sn !(g/kg)/(mole/kg)=g/mole

   !calc T
   logT=log(T)
   ERTnow=0d0;cvrnow=0d0
   call calc_vThrt( T,logT,vThrt)
   call calc_vTcpr( T,     vTcpr)
   do i=1,ns
      call check_section_number(i,T,secN)
      tmp =-1;tmp2=-1
      do j=1,9
         tmp  = tmp +co(j,secN,i)*vThrt(j)
         tmp2 = tmp2+co(j,secN,i)*vTcpr(j)
      end do
      ERTnow = ERTnow + tmp *n(i)
      cvrnow = cvrnow + tmp2*n(i)
   end do
   E=ERTnow*Ru*T

   !calc kappa
   cprnow= cvrnow+sn
   kappa=1d0+sn/cvrnow

   !calc DHi and vhi
   do i=1,ns
      call check_section_number(i,T,secN)
      tmp= 0d0
      do j=1,9
         tmp =tmp +co(j,secN,i)*vThrt(j)
      end do
      vhi_s(i)=1d0*tmp
   end do
   tmp=Ru*T
   vhi(1)=dot_product(nf,vhi_s)*tmp
   vhi(2)=dot_product(no,vhi_s)*tmp
   vhi(3)=dot_product(np,vhi_s)*tmp

   !calc mu
   call calc_vTmu(T,logT,vTmu)
   do i=1,nt
      call check_section_number_trans(i,T,sect)
      muN(i)=exp(dot_product(trans(:,sect,i),vTmu))
   end do

   mu = 0d0
   do i=1,nt
      denom = 0d0
      MWi=MW(tr2th(i))
      do j=1,nt
         MWj=MW(tr2th(j))
         phi = 0.25d0*(1d0+sqrt(muN(i)/muN(j))*(MWj/MWi)**0.25d0 )**2 &
                     *sqrt(2d0*MWj/(MWi+MWj))
         denom = denom + n(tr2th(j)) * phi
      end do
      mu = mu + n(tr2th(i)) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine flame_sheet_PTgiven!}}}

!cea
subroutine init_pack_cea(flag)!{{{
   character*2,intent(in)::flag
   call read_cheminp
   call set_therm_data
   call set_trans_data
   call read_fo_composition

   call initialize_cea(flag)
end subroutine init_pack_cea!}}}
subroutine initialize_cea(flag)!{{{
   use chem
   use chem_var
   implicit none
   character*2,intent(in)::flag
   double precision sumo,sumf,sm
   double precision Y(2),DHi(2)
   integer i,j
   integer nshifto,nshiftf
   logical flag_cea

   !set o/f ratio
   sumo=0d0
   sumf=0d0
   do j=1,ns
      sumo=sumo+MW(j)*no(j)
      sumf=sumf+MW(j)*nf(j)
   end do
   of = sumo/sumf

   !calc ini
   call calc_ini(pf,Tf, nf, Ef,MWf,Yfe)
   call calc_ini(po,To, no, Eo,MWo,Yoe)
   nfini=nf
   noini=no
 
   ! calc the amount of each element
   do j=1,ns
      do i=1,ne
         b0f(i)=b0f(i)+Ac(i,j)*nf(j)
         b0o(i)=b0o(i)+Ac(i,j)*no(j)
      end do
   end do

   !set shift and mask
   nshiftf=0
   nshifto=0
   maskbf =1d0
   maskbo =1d0
   do i=1,ne
      if(b0f(i) .eq. 0d0) then
         maskbf(i) = 0d0
         nshiftf=nshiftf+1
         nelistf(nshiftf)=i
      else
         elistf(i-nshiftf)=i
      end if
      if(b0o(i) .eq. 0d0) then
         maskbo(i) = 0d0
         nshifto=nshifto+1
         nelisto(nshifto)=i
      else
         elisto(i-nshifto)=i
      end if
   end do
   elisto(ne+1-nshifto)=ne+1
   elisto(ne+2-nshifto)=ne+2
   neo=ne+1-nshifto
   elistf(ne+1-nshiftf)=ne+1
   elistf(ne+2-nshiftf)=ne+2
   nef=ne+1-nshiftf
   do j=1,ns
      masko(j)=1d0
      maskf(j)=1d0
      do i=1,nshifto
         if(Ac(nelisto(i),j) .ne. 0d0) then
            masko(j)=0d0
         end if
      end do
      do i=1,nshiftf
         if(Ac(nelistf(i),j) .ne. 0d0) then
            maskf(j)=0d0
         end if
      end do
      if(masko(j) .ne. 1d0) masko(j)=1d-100
      if(maskf(j) .ne. 1d0) maskf(j)=1d-100
   end do

   !calc n,E
   Y(1)=1d0;Y(2)=0d0
   rhof=pf/(Ru*1d3/MWf*Tf)
   nf=nf+initial_eps
   if(flag .eq. 'uv') then
      call cea(rhof,Y,Ef, Tf,nf, MWf,kappaf,muf,Yvf,vhif,DHif)
      pf=rhof*Ru*1d3/MWf*Tf
      Hf=Ef+pf/rhof
   else
      Hf=Ef+pf/rhof
      call cea_hp(pf,Y,Hf, Tf,nf, MWf,kappaf,muf,Yvf,vhif,flag_cea)
      rhof=pf/(Ru*1d3/MWf*Tf)
      Ef=Hf-pf/rhof
   end if
   call set_static_qw(pf,rhof,Tf,Ef,kappaf,muf,Y,qf,wf)

   Y(1)=0d0;Y(2)=1d0
   rhoo=po/(Ru*1d3/MWo*To)
   no=no+initial_eps
   if(flag .eq. 'uv') then
      call cea(rhoo,Y,Eo, To,no, MWo,kappao,muo,Yvo,vhio,DHio)
      po=rhoo*Ru*1d3/MWo*To
      Ho=Eo+po/rhoo
   else
      Ho=Eo+po/rhoo
      call cea_hp(po,Y,Ho, To,no, MWo,kappao,muo,Yvo,vhio,flag_cea)
      rhoo=po/(Ru*1d3/MWo*To)
      Eo=Ho-po/rhoo
   end if
   call set_static_qw(po,rhoo,To,Eo,kappao,muo,Y,qo,wo)

   sm=sum(no(1:ns)+nf(1:ns),1)/ns
   n_save=sm
contains
   subroutine calc_ini(p,T, n, E,MWave,Yfoe)
      use func_therm
      implicit none
      double precision,intent(in) ::p
      double precision,intent(in) ::T
      double precision,intent(inout)::n(ns)
      double precision,intent(out)::E
      double precision,intent(out)::MWave
      double precision,intent(out)::Yfoe(ne)
   
      integer i,j,k,nsc
      double precision vThrt(9)
      double precision sn,sm,tmp
   
      sn=0d0
      sm=0d0
      do j=1,ns
         sn=sn+n(j)
         sm=sm+n(j)*MW(j)
      end do
      MWave=sm/sn
      n=n/sm*1d3 !mole/kg
   
      call calc_vThrt(T,log(T),vThrt)
      do j=1,ns
         call check_section_number(j,T,nsc)
         tmp=-1d0
         do k=1,9
            tmp=tmp+vThrt(k)*co(k,nsc,j)
         end do
         E=E+tmp*n(j)
      end do
      E=E*Ru*T

      Yfoe=0d0
      do j=1,ns
         tmp=n(j)*MW(j)*1d-3
         do i=1,ne
            Yfoe(i) =Yfoe(i)+YAc(i,j)*tmp
         end do
      end do
   end subroutine calc_ini
   subroutine set_static_qw(p,rho,T,E,kappa,mu,Y,q_local,w_local)
      implicit none
      double precision,intent(in)::p
      double precision,intent(in)::rho
      double precision,intent(in)::T
      double precision,intent(in)::E
      double precision,intent(in)::kappa
      double precision,intent(in)::mu
      double precision,intent(in)::Y(nY)
      double precision,intent(out)::q_local(dimq)
      double precision,intent(out)::w_local(dimw)
      integer i
      do i=1,nY
         q_local(i)   = Y(i)*rho
         w_local(i+4) = Y(i)
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
end subroutine initialize_cea!}}}
subroutine cea(rho,Y,E, T,n, MWave,kappa,mu,Yv,vhi,DHi)!{{{
   use const_chem
   use func_therm
   use chem
   use chem_var
   implicit none
   double precision,intent(in)   ::rho
   double precision,intent(in)   ::Y(2)
   double precision,intent(in)   ::E
   double precision,intent(inout)::T
   double precision,intent(inout)::n(ns)
   double precision,intent(out)  ::MWave
   double precision,intent(out)  ::kappa
   double precision,intent(out)  ::mu
   double precision,intent(out)  ::Yv(ns)
   double precision,intent(out)  ::vhi(ns)
   double precision,intent(out)  ::DHi(2)

   double precision,dimension(ne+2)::b0,b,bd,vpi
   double precision,dimension(ns)::vmurt,vert,Dlogn,logn
   integer,dimension(ns)::sec_num
   double precision,dimension(ne+2,ne+2)::Ad
   double precision tmp,logT,sn,vmurtn,DlogT,b0max,Dnmax,logrhoRu,tinyn,tinytinyn,logtinyn,logtinytinyn

   double precision dh,tmp2
   double precision vhe(ns),MWd(ne)

   integer sect
   double precision,dimension(ns)::muN
   double precision phi,denom,MWj,MWi

   integer i,j,k,counter
   logical flag,LUflag,REDUCEflag,debugflag
   integer nen,elist(ne+2),nelist(ne+2)

   double precision lambda1,lambda2,omega_var,logsn

   double precision,dimension(9)::vTmurt,vTcpr,vThrt
   double precision,dimension(4)::vTmu


   debugflag= .false.
   logrhoRu = log(rho*Ru/pst)

   !calc b0 of elements
   b0= Y(1)*b0f + Y(2)*b0o
   b0max = maxval(b0(1:ne),1)

   if(Y(1)<Y_eps) then
      REDUCEflag=.true.
      nen=neo
      elist =elisto
      nelist=nelisto
      n     =n *masko
      b0    =b0*maskbo
   else if(Y(2)<Y_eps) then
      REDUCEflag=.true.
      nen=nef
      elist =elistf
      nelist=nelistf
      n     =n *maskf
      b0    =b0*maskbf
   else
      REDUCEflag=.false.
   end if

   sn=sum(n(1:ns),1)
   tinyn        = TSIZE*sn
   logtinyn     = log(tinyn)
   tinytinyn    = TTSIZE*sn
   logtinytinyn = log(tinytinyn)

   !set vert
   logT=log(T)
   call calc_vThrt(T,logT,vThrt)
   b=0d0
   vert(1:ns)=0d0
   logn(1:ns)=logtinyn
   do j=1,ns
      call check_section_number(j,T,sec_num(j))
      if(n(j)>tinyn) then
         vert(j)  = dot_product(co(:,sec_num(j),j),vThrt)-1d0
         logn(j)  = log(n(j))
         do i=1,ne
            b(i) = b(i) + Ac(i,j)*n(j)
         end do
         b(ne+1)=b(ne+1)+vert(j)*n(j)
      end if
   end do

   !check convergence
   flag = .true.
   do i=1,ne
      if(b0(i) > eps .and. abs(b0(i)-b(i))>b0max*eps) flag = .false.
   end do
   if(abs(E-b(ne+1)*Ru*T) >eps*abs(E)) flag = .false.

   if(.not. flag) then
      counter=1
      !open(33,file="debug")
      do
         !set b0(ne+1)
         b0(ne+1)=E/(Ru*T)

         !set Ad, b and vmurtn
         b(1:ne)=0d0
         bd=b0
         bd(ne+1)=bd(ne+1)-b(ne+1)
         Ad=0d0
         call calc_vTmurt(T,logT,vTmurt)
         call calc_vTcpr(T,vTcpr)
         do k=1,ns
            vmurt(k)=dot_product(co(:,sec_num(k),k),vTmurt)+logrhoRu+logT+logn(k)

            if(n(k)>tinyn) then
               vmurtn=vmurt(k)*n(k)
               do i=1,ne
                  do j=i,ne
                     Ad(i,j)=Ad(i,j)+Ac(i,k)*Ac(j,k)*n(k)
                  end do
                  Ad(i,ne+1)=Ad(i,ne+1)+Ac(i,k)*vert(k)*n(k)
                  b(i)  = b(i)  + Ac(i,k)* n(k)
                  bd(i) = bd(i) + Ac(i,k)*(vmurtn-n(k))
               end do

               Ad(ne+1,ne+1)= Ad(ne+1,ne+1)+ n(k) *(vert(k)**2+dot_product(co(:,sec_num(k),k),vTcpr)-1d0)
               bd(ne+1)     = bd(ne+1)     + vmurtn*vert(k)
            end if
         end do

         do i=1,ne+1
            do j=i,ne+1
               Ad(j,i)=Ad(i,j)
            end do
         end do

         !calc vpi
         if(REDUCEflag) then
            do i=1,nen
               do j=1,nen
                  Ad(i,j)=Ad(elist(i),elist(j))
               end do
               bd(i)=bd(elist(i))
            end do
            call LU(Ad,bd,vpi,nen,LUflag)
            do i=nen,1,-1
               vpi(elist(i))=vpi(i)
            end do
            do i=1,ne+1-nen
               vpi(nelist(i))=-1d300
            end do
         else
            call LU(Ad,bd,vpi,ne+1,LUflag)
         end if

         if(LUflag) then
            debugflag=.true.
            !set new n
            n=n+initial_eps

            !set new n
            b(ne+1)=0d0
            Dnmax = 0d0
            logn(1:ns) = logtinytinyn
            do j=1,ns
               if(n(j)<tinytinyn) then
                   n(j)    = tinytinyn
               else
                   logn(j) = log(n(j))
               end if

               if(n(j)>tinyn) then
                   vert(j)  = dot_product(co(:,sec_num(j),j),vThrt)-1d0
                   b(ne+1) = b(ne+1)+vert(j)*n(j)
               end if
            end do
            sn = sum(n(1:ns),1)
         end if

         !calc Delta
         DlogT=vpi(ne+1)
         !write(33,*) counter,abs(DlogT)
         do j=1,ns
            tmp=0d0
            do i=1,ne
               tmp=tmp+Ac(i,j)*vpi(i)
            end do
            Dlogn(j)= -vmurt(j)+tmp+vert(j)*DlogT
         end do

         !calc omega_var
         lambda1=0d0
         lambda2=1d300
         logsn  =log(sn)
         do j=1,ns
            if(logn(j)-logsn>-18.420681d0) then
               lambda1=max(lambda1,abs(Dlogn(j)))
            else if(Dlogn(j)>0d0) then
               lambda2=min(lambda2,abs((-logn(j)+logsn-9.2103404d0)/Dlogn(j)))
            end if
         end do
         lambda1=2d0/max(5d0*abs(DlogT),abs(lambda1))
         omega_var=min(1d0,min(lambda1,lambda2))

         !set new T
         T=T*exp(omega_var*DlogT)
         logT=log(T)
         do j=1,ns
            call check_section_number(j,T,sec_num(j))
         end do

         !set new n
         b(ne+1)=0d0
         Dnmax = 0d0
         logn(1:ns) = logtinytinyn
         call calc_vThrt(T,logT,vThrt)
         do j=1,ns
            !set vert
            n(j)=n(j)*exp(omega_var*Dlogn(j))
            if(n(j)<tinytinyn) then
                n(j)    = tinytinyn
            else
                logn(j) = log(n(j))
            end if

            if(n(j)>tinyn) then
                vert(j)  = dot_product(co(:,sec_num(j),j),vThrt)-1d0
                b(ne+1) = b(ne+1)+vert(j)*n(j)
                tmp = abs(Dlogn(j))*n(j)
                if(Dnmax < tmp) Dnmax = tmp
            end if
         end do
         sn = sum(n(1:ns),1)

         !check convergence
         flag = .true.
         if(abs(Dnmax)>eps*sn) then
            flag = .false.
            !if(counter .eq. 501) print *,"flag1"
         end if
         do i=1,ne
            if(b0(i) > eps .and. abs(b0(i)-b(i))>b0max*eps)       then
               flag = .false.
               !if(counter .eq. 501) print *,"flag2",i
            end if
         end do
         if(abs(DlogT)>eps .and. abs(E-b(ne+1)*Ru*T) >eps*abs(E)) then
            flag = .false.
            !if(counter .eq. 501) print *,"flag3"
         end if
         if(counter>500) flag = .true.
         if(flag) exit

         counter=counter+1
      end do

      !close(33)
      if(counter>500) then
         !print *,E,T,Y,REDUCEflag,debugflag
         print *,"not converted. at cea"
      end if
   end if

   MWave=1d3/sn

   !calc kappa
   call calc_vTcpr(T,vTcpr)
   kappa=0d0
   do i=1,ns
      if(n(i)>tinyn) then
         kappa=kappa+n(i)*dot_product(co(:,sec_num(i),i),vTcpr)
      end if
   end do
   kappa=kappa/(kappa-sn)

   !calc vhi and DHi
   tmp=Ru*T
   vhe=0d0
   MWd=0d0
   do i=1,ns
      dh= 0d0
      do j=1,9
         dh =dh +co(j,sec_num(i),i)*vThrt(j)
      end do
      vhi( i)=tmp/MW(i)*dh*1d3
      tmp2   =tmp*n( i)*dh
      do j=1,ne
         vhe(j)=vhe(j)+YAc(j,i)*tmp2
         MWd(j)=MWd(j)+YAc(j,i)*n(i)
      end do
   end do

   DHi=0d0
   do i=1,ne
      tmp =1d0/(Y(1)*Yfe(i)+Y(2)*Yoe(i)+1d-300)&
          *(kappa*Ru*MWd(i)*T-(kappa-1d0)*vhe(i))
      DHi(1)=DHi(1)+Yfe(i)*tmp
      DHi(2)=DHi(2)+Yoe(i)*tmp
   end do

   !calc Yv
   Yv=n*MW*1d-3

   !calc mu
   call calc_vTmu(T,logT,vTmu)
   do i=1,nt
      call check_section_number_trans(i,T,sect)
      muN(i)=exp(dot_product(trans(:,sect,i),vTmu))
   end do

   mu = 0d0
   do i=1,nt
      denom = 0d0
      MWi=MW(tr2th(i))
      do j=1,nt
         MWj=MW(tr2th(j))
         phi = 0.25d0*(1d0+sqrt(muN(i)/muN(j))*(MWj/MWi)**0.25d0 )**2 &
                     *sqrt(2d0*MWj/(MWi+MWj))
         denom = denom + n(tr2th(j)) * phi
      end do
      mu = mu + n(tr2th(i)) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine cea!}}}
subroutine cea_hp(p,Y,H, T,n, MWave,kappa,mu,Yv,vhi,flag_outer)!{{{
   use const_chem
   use func_therm
   use chem
   use chem_var
   implicit none
   double precision,intent(in)   ::p
   double precision,intent(in)   ::Y(2)
   double precision,intent(in)   ::H
   double precision,intent(inout)::T
   double precision,intent(inout)::n(max_ns)
   double precision,intent(out)  ::MWave
   double precision,intent(out)  ::kappa
   double precision,intent(out)  ::mu
   double precision,intent(out)  ::Yv(max_ns)
   double precision,intent(out)  ::vhi(max_ns)
   logical         ,intent(out)  ::flag_outer

   double precision,dimension(ne+2)::b0,b,bd,vpi
   double precision,dimension(max_ns)::vmurt,vhrt,Dlogn,logn
   integer,dimension(max_ns)::sec_num
   double precision,dimension(ne+2,ne+2)::Ad
   double precision logP,tmp,logT,sn,logsn,Dlogsn,vmurtn,DlogT,b0max,Dnmax,logrhoRu,tinyn,tinytinyn,logtinyn,logtinytinyn

   double precision dh

   integer sect
   double precision,dimension(max_ns)::muN
   double precision phi,denom,MWj,MWi

   integer i,j,k,counter
   logical flag,LUflag,REDUCEflag
   integer nen,elist(ne+2),nelist(ne+2)

   double precision lambda1,lambda2,omega_var

   double precision,dimension(9)::vTmurt,vTcpr,vThrt
   double precision,dimension(4)::vTmu


   flag_outer=.true.
   logP = log(p/pst)

   !calc b0 of elements
   b0= Y(1)*b0f + Y(2)*b0o
   b0max = maxval(b0(1:ne),1)
   sn=sum(n(1:ns),1)
   b0(ne+1)=sn

   if(Y(1)<Y_eps) then
      REDUCEflag=.true.
      nen=neo+1
      elist =elisto
      nelist=nelisto
      n     =n*masko
   else if(Y(2)<Y_eps) then
      REDUCEflag=.true.
      nen=nef+1
      elist =elistf
      nelist=nelistf
      n     =n*maskf
   else
      REDUCEflag=.false.
   end if

   tinyn        = TSIZE*sn
   logtinyn     = log(tinyn)
   tinytinyn    = TTSIZE*sn
   logtinytinyn = log(tinytinyn)

   !set vhrt
   logT=log(T)
   call calc_vThrt(T,logT,vThrt)
   b=0d0
   b(ne+1) =sn
   vhrt(1:ns)=0d0
   logn(1:ns)=logtinyn
   do j=1,ns
      call check_section_number(j,T,sec_num(j))
      if(n(j)>tinyn) then
         vhrt(j)  = dot_product(co(:,sec_num(j),j),vThrt)
         logn(j)  = log(n(j))
         do i=1,ne
            b(i) = b(i) + Ac(i,j)*n(j)
         end do
         b(ne+2)=b(ne+2)+vhrt(j)*n(j)
      end if
   end do

   !check convergence
   flag = .true.
   do i=1,ne
      if(b0(i) > eps .and. abs(b0(i)-b(i))>b0max*eps) flag = .false.
   end do
   if(abs(H-b(ne+2)*Ru*T) >eps*abs(H)) flag = .false.

   if(.not. flag) then
      counter=1
      do
         !set b0(ne+2)
         b0(ne+2)=H/(Ru*T)

         !set Ad, b and vmurtn
         b(1:ne)=0d0
         bd=b0
         bd(ne+1)=b0(ne+1)-b(ne+1)
         bd(ne+2)=b0(ne+2)-b(ne+2)
         Ad=0d0
         call calc_vTmurt(T,logT,vTmurt)
         call calc_vTcpr(T,vTcpr)
         do k=1,ns
            vmurt(k)=dot_product(co(:,sec_num(k),k),vTmurt)+logn(k)-log(b0(ne+1))+logP

            if(n(k)>tinyn) then
               vmurtn=vmurt(k)*n(k)
               do i=1,ne
                  do j=i,ne
                     Ad(i,j)=Ad(i,j)+Ac(i,k)*Ac(j,k)*n(k)
                  end do
                  Ad(i,ne+1)= Ad(i,ne+1)+ Ac(i,k)*        n(k)
                  Ad(i,ne+2)= Ad(i,ne+2)+ Ac(i,k)*vhrt(k)*n(k)
                  b(i)      = b(i)      + Ac(i,k)        *n(k)
                  bd(i)     = bd(i)     + Ac(i,k)*(vmurtn-n(k))
               end do

               Ad(ne+1,ne+2)= Ad(ne+1,ne+2)+ n(k) * vhrt(k)
               Ad(ne+2,ne+2)= Ad(ne+2,ne+2)+ n(k) *(vhrt(k)**2+dot_product(co(:,sec_num(k),k),vTcpr))
               bd(ne+1)     = bd(ne+1)     + vmurtn
               bd(ne+2)     = bd(ne+2)     + vmurtn*vhrt(k)
            end if
         end do
         Ad(ne+1,ne+1)= b(ne+1)-b0(ne+1)

         do i=1,ne+2
            do j=i+1,ne+2
               Ad(j,i)=Ad(i,j)
            end do
         end do

         !calc vpi
         if(REDUCEflag) then
            do i=1,nen
               do j=1,nen
                  Ad(i,j)=Ad(elist(i),elist(j))
               end do
               bd(i)=bd(elist(i))
            end do
            call LU(Ad,bd,vpi,nen,LUflag)
            do i=nen,1,-1
               vpi(elist(i))=vpi(i)
            end do
            do i=1,ne+2-nen
               vpi(nelist(i))=-1d300
            end do
         else
            call LU(Ad,bd,vpi,ne+2,LUflag)
         end if

         if(LUflag) then
            !set new n
            n=n+initial_eps

            !set new n
            sn=0d0
            b(ne+2)=0d0
            logn(1:ns) = logtinytinyn
            do j=1,ns
               if(n(j)<tinytinyn) then
                   n(j)    = tinytinyn
               else
                   logn(j) = log(n(j))
               end if

               if(n(j)>tinyn) then
                   vhrt(j)  = dot_product(co(:,sec_num(j),j),vThrt)
                   b(ne+2)  = b(ne+2)+vhrt(j)*n(j)
                   sn       = sn     +        n(j)
               end if
            end do
            b(ne+1) = sn
            counter=counter+1
            cycle
         end if

         !calc Delta
         Dlogsn=vpi(ne+1)
         DlogT =vpi(ne+2)
         do j=1,ns
            tmp=0d0
            do i=1,ne
               tmp=tmp+Ac(i,j)*vpi(i)
            end do
            Dlogn(j)= -vmurt(j)+tmp+Dlogsn+vhrt(j)*DlogT
         end do

         !calc omega_var
         lambda1=0d0
         lambda2=1d300
         logsn  =log(b0(ne+1))
         do j=1,ns
            if(logn(j)-logsn>-18.420681d0) then
               lambda1=max(lambda1,abs(Dlogn(j)))
            else if(Dlogn(j)>0d0) then
               lambda2=min(lambda2,abs((-logn(j)+logsn-9.2103404d0)/(Dlogn(j)-Dlogsn)))
            end if
         end do
         lambda1=2d0/max(5d0*abs(DlogT),5d0*abs(Dlogsn),abs(lambda1))
         omega_var=min(1d0,min(lambda1,lambda2))

         !set new sn
         b0(ne+1)=b0(ne+1)*exp(omega_var*Dlogsn)

         !set new T
         T=T*exp(omega_var*DlogT)
         logT=log(T)
         do j=1,ns
            call check_section_number(j,T,sec_num(j))
         end do

         !set new n
         sn =0d0
         b(ne+2)=0d0
         Dnmax = 0d0
         logn(1:ns) = logtinytinyn
         call calc_vThrt(T,logT,vThrt)
         do j=1,ns
            !set vhrt
            n(j)=n(j)*exp(omega_var*Dlogn(j))
            if(n(j)<tinytinyn) then
                n(j)    = tinytinyn
            else
                logn(j) = log(n(j))
            end if

            if(n(j)>tinyn) then
                vhrt(j)  = dot_product(co(:,sec_num(j),j),vThrt)
                b(ne+2) = b(ne+2)+vhrt(j)*n(j)
                sn      = sn     +        n(j)
                tmp = abs(Dlogn(j))*n(j)
                if(Dnmax < tmp) Dnmax = tmp
            end if
         end do
         b(ne+1) = sn

         !check convergence
         flag = .true.
         if(abs(Dnmax)>eps*sn) flag = .false.
         do i=1,ne
            if(b0(i) > eps .and. abs(b0(i)-b(i))>b0max*eps)               flag = .false.
         end do
         if(abs(Dlogsn)>eps .or. abs(b0(ne+1)-b(ne+1))      >eps*abs(sn)) flag = .false.
         if(abs(DlogT) >eps .or. abs(H       -b(ne+2)*Ru*T) >eps*abs(H) ) flag = .false.
         if(counter>500) flag = .true.
         if(flag) exit

         counter=counter+1
      end do

      if(counter>500) then
         print *,"not converted. at cea_hp"
         flag_outer=.false.
      end if
   end if

   MWave=1d3/sn

   !calc kappa
   call calc_vTcpr(T,vTcpr)
   kappa=0d0
   do i=1,ns
      if(n(i)>tinyn) then
         kappa=kappa+n(i)*dot_product(co(:,sec_num(i),i),vTcpr)
      end if
   end do
   kappa=kappa/(kappa-sn)

   !calc vhi
   tmp=Ru*1d3*T
   do i=1,ns
      dh= 0d0
      do j=1,9
         dh =dh +co(j,sec_num(i),i)*vThrt(j)
      end do
      vhi(i)=tmp/MW(i)*dh
   end do
   tmp=Ru*T

   !calc Yv
   Yv=n*MW*1d-3

   !calc mu
   call calc_vTmu(T,logT,vTmu)
   do i=1,nt
      call check_section_number_trans(i,T,sect)
      muN(i)=exp(dot_product(trans(:,sect,i),vTmu))
   end do

   mu = 0d0
   do i=1,nt
      denom = 0d0
      MWi=MW(tr2th(i))
      do j=1,nt
         MWj=MW(tr2th(j))
         phi = 0.25d0*(1d0+sqrt(muN(i)/muN(j))*(MWj/MWi)**0.25d0 )**2 &
                     *sqrt(2d0*MWj/(MWi+MWj))
         denom = denom + n(tr2th(j)) * phi
      end do
      mu = mu + n(tr2th(i)) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine cea_hp!}}}
subroutine cea_tp(p,Y,T, n, H,MWave,kappa,mu,Yv,vhi,flag_outer)!{{{
   use const_chem
   use func_therm
   use chem
   use chem_var
   implicit none
   double precision,intent(in)   ::p
   double precision,intent(in)   ::Y(2)
   double precision,intent(in)   ::T
   double precision,intent(inout)::n(max_ns)
   double precision,intent(out)  ::H
   double precision,intent(out)  ::MWave
   double precision,intent(out)  ::kappa
   double precision,intent(out)  ::mu
   double precision,intent(out)  ::Yv(max_ns)
   double precision,intent(out)  ::vhi(max_ns)
   logical         ,intent(out)  ::flag_outer

   double precision,dimension(ne+2)::b0,b,bd,vpi
   double precision,dimension(max_ns)::vmurt,vhrt,Dlogn,logn
   integer,dimension(max_ns)::sec_num
   double precision,dimension(ne+2,ne+2)::Ad
   double precision logP,tmp,logT,sn,logsn,Dlogsn,vmurtn,DlogT,b0max,Dnmax,logrhoRu,tinyn,tinytinyn,logtinyn,logtinytinyn

   double precision dh

   integer sect
   double precision,dimension(max_ns)::muN
   double precision phi,denom,MWj,MWi

   integer i,j,k,counter
   logical flag,LUflag,REDUCEflag
   integer nen,elist(ne+2),nelist(ne+2)

   double precision lambda1,lambda2,omega_var

   double precision,dimension(9)::vTmurt,vTcpr,vThrt
   double precision,dimension(4)::vTmu


   flag_outer=.true.
   logP = log(p/pst)

   !calc b0 of elements
   b0= Y(1)*b0f + Y(2)*b0o
   b0max = maxval(b0(1:ne),1)
   sn=sum(n(1:ns),1)
   b0(ne+1)=sn

   if(Y(1)<Y_eps) then
      REDUCEflag=.true.
      nen=neo
      elist =elisto
      nelist=nelisto
      n     =n*masko
   else if(Y(2)<Y_eps) then
      REDUCEflag=.true.
      nen=nef
      elist =elistf
      nelist=nelistf
      n     =n*maskf
   else
      REDUCEflag=.false.
   end if

   tinyn        = TSIZE*sn
   logtinyn     = log(tinyn)
   tinytinyn    = TTSIZE*sn
   logtinytinyn = log(tinytinyn)

   !set vhrt
   logT=log(T)
   b=0d0
   b(ne+1) =sn
   vhrt(1:ns)=0d0
   logn(1:ns)=logtinyn
   do j=1,ns
      call check_section_number(j,T,sec_num(j))
      if(n(j)>tinyn) then
         logn(j)  = log(n(j))
         do i=1,ne
            b(i) = b(i) + Ac(i,j)*n(j)
         end do
      end if
   end do

   !check convergence
   flag = .true.
   do i=1,ne
      if(b0(i) > eps .and. abs(b0(i)-b(i))>b0max*eps) flag = .false.
   end do

   if(.not. flag) then
      counter=1
      do
         !set Ad, b and vmurtn
         b(1:ne)=0d0
         bd=b0
         bd(ne+1)=b0(ne+1)-b(ne+1)
         Ad=0d0
         call calc_vTmurt(T,logT,vTmurt)
         do k=1,ns
            vmurt(k)=dot_product(co(:,sec_num(k),k),vTmurt)+logn(k)-log(b0(ne+1))+logP

            if(n(k)>tinyn) then
               vmurtn=vmurt(k)*n(k)
               do i=1,ne
                  do j=i,ne
                     Ad(i,j)=Ad(i,j)+Ac(i,k)*Ac(j,k)*n(k)
                  end do
                  Ad(i,ne+1)= Ad(i,ne+1)+ Ac(i,k)*        n(k)
                  b(i)      = b(i)      + Ac(i,k)        *n(k)
                  bd(i)     = bd(i)     + Ac(i,k)*(vmurtn-n(k))
               end do

               bd(ne+1)     = bd(ne+1)     + vmurtn
            end if
         end do
         Ad(ne+1,ne+1)= b(ne+1)-b0(ne+1)

         do i=1,ne+1
            do j=i+1,ne+1
               Ad(j,i)=Ad(i,j)
            end do
         end do

         !calc vpi
         if(REDUCEflag) then
            do i=1,nen
               do j=1,nen
                  Ad(i,j)=Ad(elist(i),elist(j))
               end do
               bd(i)=bd(elist(i))
            end do
            call LU(Ad,bd,vpi,nen,LUflag)
            do i=nen,1,-1
               vpi(elist(i))=vpi(i)
            end do
            do i=1,ne+1-nen
               vpi(nelist(i))=-1d300
            end do
         else
            call LU(Ad,bd,vpi,ne+1,LUflag)
         end if

         if(LUflag) then
            !set new n
            n=n+initial_eps

            !set new n
            sn=0d0
            logn(1:ns) = logtinytinyn
            do j=1,ns
               if(n(j)<tinytinyn) then
                   n(j)    = tinytinyn
               else
                   logn(j) = log(n(j))
               end if

               if(n(j)>tinyn) then
                   vhrt(j)  = dot_product(co(:,sec_num(j),j),vThrt)
                   sn       = sn     +        n(j)
               end if
            end do
            b(ne+1) = sn
            counter=counter+1
            cycle
         end if

         !calc Delta
         Dlogsn=vpi(ne+1)
         do j=1,ns
            tmp=0d0
            do i=1,ne
               tmp=tmp+Ac(i,j)*vpi(i)
            end do
            Dlogn(j)= -vmurt(j)+tmp+Dlogsn
         end do

         !calc omega_var
         lambda1=0d0
         lambda2=1d300
         logsn  =log(b0(ne+1))
         do j=1,ns
            if(logn(j)-logsn>-18.420681d0) then
               lambda1=max(lambda1,abs(Dlogn(j)))
            else if(Dlogn(j)>0d0) then
               lambda2=min(lambda2,abs((-logn(j)+logsn-9.2103404d0)/(Dlogn(j)-Dlogsn)))
            end if
         end do
         lambda1=2d0/max(5d0*abs(Dlogsn),abs(lambda1))
         omega_var=min(1d0,min(lambda1,lambda2))

         !set new sn
         b0(ne+1)=b0(ne+1)*exp(omega_var*Dlogsn)

         !set new n
         sn =0d0
         Dnmax = 0d0
         logn(1:ns) = logtinytinyn
         do j=1,ns
            !set vhrt
            n(j)=n(j)*exp(omega_var*Dlogn(j))
            if(n(j)<tinytinyn) then
                n(j)    = tinytinyn
            else
                logn(j) = log(n(j))
            end if

            if(n(j)>tinyn) then
                sn      = sn     +        n(j)
                tmp = abs(Dlogn(j))*n(j)
                if(Dnmax < tmp) Dnmax = tmp
            end if
         end do
         b(ne+1) = sn

         !check convergence
         flag = .true.
         if(abs(Dnmax)>eps*sn) flag = .false.
         do i=1,ne
            if(b0(i) > eps .and. abs(b0(i)-b(i))>b0max*eps)               flag = .false.
         end do
         if(abs(Dlogsn)>eps .or. abs(b0(ne+1)-b(ne+1))      >eps*abs(sn)) flag = .false.
         if(counter>500) flag = .true.
         if(flag) exit

         counter=counter+1
      end do

      if(counter>500) then
         print *,"not converted. at cea_hp"
         flag_outer=.false.
      end if
   end if

   MWave=1d3/sn

   !calc kappa
   call calc_vTcpr(T,vTcpr)
   kappa=0d0
   do i=1,ns
      if(n(i)>tinyn) then
         kappa=kappa+n(i)*dot_product(co(:,sec_num(i),i),vTcpr)
      end if
   end do
   kappa=kappa/(kappa-sn)

   !calc H and vhi
   call calc_vThrt(T,logT,vThrt)
   H=0d0
   tmp=Ru*1d3*T
   do i=1,ns
      dh= 0d0
      do j=1,9
         dh =dh +co(j,sec_num(i),i)*vThrt(j)
      end do
      vhi(i)=tmp/MW(i)*dh
      H     =H+   n(i)*dh
   end do
   H=H*Ru*T

   !calc Yv
   Yv=n*MW*1d-3

   !calc mu
   call calc_vTmu(T,logT,vTmu)
   do i=1,nt
      call check_section_number_trans(i,T,sect)
      muN(i)=exp(dot_product(trans(:,sect,i),vTmu))
   end do

   mu = 0d0
   do i=1,nt
      denom = 0d0
      MWi=MW(tr2th(i))
      do j=1,nt
         MWj=MW(tr2th(j))
         phi = 0.25d0*(1d0+sqrt(muN(i)/muN(j))*(MWj/MWi)**0.25d0 )**2 &
                     *sqrt(2d0*MWj/(MWi+MWj))
         denom = denom + n(tr2th(j)) * phi
      end do
      mu = mu + n(tr2th(i)) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine cea_tp!}}}
subroutine calcH(T,Y,H)!{{{
   use const_chem
   use func_therm
   use chem
   use chem_var
   implicit none
   double precision,intent(in) ::T
   double precision,intent(in) ::Y(2)
   double precision,intent(out)::H

   double precision,dimension(max_ns)::n
   double precision,dimension(9)::vThrt
   integer sec_num,j

   n= Y(1)*nfini + Y(2)*noini
   call calc_vThrt(T,log(T),vThrt)
   H=0d0
   do j=1,ns
      call check_section_number(j,T,sec_num)
      H=H+dot_product(co(:,sec_num,j),vThrt)*n(j)
   end do
   H=H*Ru*T
end subroutine calcH!}}}
subroutine calcE_from_p(p,T,Y,rho,E)!{{{
   use const_chem
   use func_therm
   use chem
   use chem_var
   implicit none
   double precision,intent(in) ::p
   double precision,intent(in) ::T
   double precision,intent(in) ::Y(2)
   double precision,intent(out)::rho
   double precision,intent(out)::E

   double precision,dimension(max_ns)::n
   double precision,dimension(9)::vThrt
   integer sec_num,j

   n= Y(1)*nfini + Y(2)*noini
   rho=p/(sum(n(1:ns),1)*Ru*T)
   call calc_vThrt(T,log(T),vThrt)
   E=0d0
   do j=1,ns
      call check_section_number(j,T,sec_num)
      E=E+(dot_product(co(:,sec_num,j),vThrt)-1d0)*n(j)
   end do
   E=E*Ru*T
end subroutine calcE_from_p!}}}
subroutine calcE(T,Y,E)!{{{
   use const_chem
   use func_therm
   use chem
   use chem_var
   implicit none
   double precision,intent(in) ::T
   double precision,intent(in) ::Y(2)
   double precision,intent(out)::E

   double precision,dimension(max_ns)::n
   double precision,dimension(9)::vThrt
   integer sec_num,j

   n= Y(1)*nfini + Y(2)*noini
   call calc_vThrt(T,log(T),vThrt)
   E=0d0
   do j=1,ns
      call check_section_number(j,T,sec_num)
      E=E+(dot_product(co(:,sec_num,j),vThrt)-1d0)*n(j)
   end do
   E=E*Ru*T
end subroutine calcE!}}}
