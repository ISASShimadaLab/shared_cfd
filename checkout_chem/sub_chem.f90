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

   integer i,j,k,l,Nfound

   if(myid .eq. 0) then
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

   call MPI_Bcast(num_sctn,      ns,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(      MW,      ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  Trange,  ns*6*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(      co,  ns*6*9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(      Ac,   ne*ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
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
   call flame_sheet(Y,Ef,Tf, MWf,kappaf,muf,DHi,Yvf,vhif)
   rhof=pf/(Ru*1d3/MWf*Tf)
   call set_static_qw(pf,rhof,Tf,Ef,kappaf,muf,Y,qf,wf)

   Y(1)=0d0;Y(2)=1d0
   call calc_ini(po,To, no, Eo)
   call flame_sheet(Y,Eo,To, MWo,kappao,muo,DHi,Yvo,vhio)
   rhoo=po/(Ru*1d3/MWo*To)
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
   double precision,intent(in)::Y(2)
   double precision,intent(in)::E
   double precision,intent(inout)::T
   double precision,intent(out)::MWave
   double precision,intent(out)::kappa
   double precision,intent(out)::mu
   double precision,intent(out)::DHi(2)
   double precision,intent(out)::Yv(3)
   double precision,intent(out)::vhi(3)

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
      DHi_s(i)=1d0*(kappa-(kappa-1d0)*tmp)
      vhi_s(i)=1d0*tmp
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

!cea
subroutine init_pack_cea!{{{
   call read_cheminp
   call set_therm_data
   call set_trans_data
   call read_fo_composition

   call initialize_cea
end subroutine init_pack_cea!}}}
subroutine initialize_cea!{{{
   use chem
   use chem_var
   implicit none
   double precision sumo,sumf,sm
   double precision Y(2),DHi(2)
   integer i,j
   integer nshifto,nshiftf

   !set o/f ratio
   sumo=0d0
   sumf=0d0
   do j=1,ns
      sumo=sumo+MW(j)*no(j)
      sumf=sumf+MW(j)*nf(j)
   end do
   of = sumo/sumf

   !calc ini
   call calc_ini(pf,Tf, nf, Ef,MWf)
   call calc_ini(po,To, no, Eo,MWo)
 
   ! calc the amount of each element
   do j=1,ns
      do i=1,ne
         b0f(i)=b0f(i)+Ac(i,j)*nf(j)
         b0o(i)=b0o(i)+Ac(i,j)*no(j)
      end do
   end do

   nshiftf=0
   nshifto=0
   do i=1,ne
      if(b0f(i) .eq. 0d0) then
         nshiftf=nshiftf+1
         nelistf(nshiftf)=i
      else
         elistf(i-nshiftf)=i
      end if
      if(b0o(i) .eq. 0d0) then
         nshifto=nshifto+1
         nelisto(nshifto)=i
      else
         elisto(i-nshifto)=i
      end if
   end do
   elisto(ne+1-nshifto)=ne+1
   neo=ne+1-nshifto
   elistf(ne+1-nshiftf)=ne+1
   nef=ne+1-nshiftf

   !calc n,E
   Y(1)=1d0;Y(2)=0d0
   rhof=pf/(Ru*1d3/MWf*Tf)
   call flame_sheet(Y,Ef,Tf, MWf,kappaf,muf,DHi,Yvf,vhif)
   call set_static_qw(pf,rhof,Tf,Ef,kappaf,muf,Y,qf,wf)

   Y(1)=0d0;Y(2)=1d0
   rhoo=po/(Ru*1d3/MWo*To)
   call flame_sheet(Y,Eo,To, MWo,kappao,muo,DHi,Yvo,vhio)
   call set_static_qw(po,rhoo,To,Eo,kappao,muo,Y,qo,wo)
contains
   subroutine calc_ini(p,T, n, E,MWave)
      use func_therm
      implicit none
      double precision,intent(in) ::p
      double precision,intent(in) ::T
      double precision,intent(inout)::n(ns)
      double precision,intent(out)::E
      double precision,intent(out)::MWave
   
      integer j,k,nsc
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
subroutine cea(rho,Y,E, T,n, MWave,kappa,mu)!{{{
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

   double precision,dimension(ne+1)::b0,b,bd,vpi
   double precision,dimension(ns)::vmurt,vert,Dlogn,logn
   integer,dimension(ns)::sec_num
   double precision,dimension(ne+1,ne+1)::Ad
   double precision tmp,logT,sn,vmurtn,DlogT,b0max,Dnmax,logrhoRu,tinyn,logtinyn
   double precision,parameter::tinyratio=1d-3
   double precision,parameter::logtinyratio=-2.3026d0*3d0 !log(tinyratio)

   integer sect
   double precision,dimension(ns)::muN
   double precision phi,denom,MWj,MWi

   integer i,j,k,counter
   logical flag,LUflag,REDUCEflag
   integer nen,elist(ne+1),nelist(ne+1)

   double precision,dimension(9)::vTmurt,vTcpr,vThrt
   double precision,dimension(4)::vTmu


   logrhoRu = log(rho*Ru/pst)

   !calc b0 of elements
   b0= Y(1)*b0f + Y(2)*b0o
   b0max = maxval(b0(1:ne),1)

   if(Y(1)<1d-6) then
      REDUCEflag=.true.
      nen=neo
      elist =elisto
      nelist=nelisto
   else if(Y(2)<1d-6) then
      REDUCEflag=.true.
      nen=nef
      elist =elistf
      nelist=nelistf
   else
      REDUCEflag=.false.
   end if

   sn=sum(n(1:ns),1)
   tinyn =TSIZE*sn
   logtinyn = log(tinyn)

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
      if(b0(i) > 1d-6 .and. abs(b0(i)-b(i))>b0max*1d-6) flag = .false.
   end do
   if(abs(E-b(ne+1)*Ru*T)/abs(E) >eps) flag = .false.

   if(.not. flag) then
      counter=1
      do
         !set b0(ne+1)
         b0(ne+1)=E/(Ru*T)

         !set Ad, b and vmurtn
         b(1:ne)=0d0
         bd=b0
         bd(ne+1)=bd(ne+1)-b(ne+1)
         Ad(1:ne+1,1:ne+1)=0d0
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
            !set new n
            n=n+initial_eps

            !set new n
            b(ne+1)=0d0
            Dnmax = 0d0
            logn(1:ns) = logtinyn + logtinyratio
            do j=1,ns
               if(n(j)<tinyn*tinyratio) then
                   n(j)    = tinyn*tinyratio
               else
                   logn(j) = log(n(j))
               end if

               if(n(j)>tinyn) then
                   vert(j)  = dot_product(co(:,sec_num(j),j),vThrt)-1d0
                   b(ne+1) = b(ne+1)+vert(j)*n(j)
               end if
            end do
            sn = sum(n(1:ns),1)
            cycle
         end if

         !calc Delta
         DlogT=vpi(ne+1)
         do j=1,ns
            tmp=0d0
            do i=1,ne
               tmp=tmp+Ac(i,j)*vpi(i)
            end do
            Dlogn(j)= -vmurt(j)+tmp+vert(j)*DlogT
         end do

         !set new T
         T=T*(1d0+omega*DlogT)
         if(T < 0d0) stop "Negative Temperature."
         logT=log(T)
         do j=1,ns
            call check_section_number(j,T,sec_num(j))
         end do

         !set new n
         b(ne+1)=0d0
         Dnmax = 0d0
         logn(1:ns) = logtinyn + logtinyratio
         call calc_vThrt(T,logT,vThrt)
         do j=1,ns
            !set vert
            n(j)=n(j)*(1d0+omega*Dlogn(j))
            if(n(j)<tinyn*tinyratio) then
                n(j)    = tinyn*tinyratio
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
         if(abs(Dnmax)>1d-5*sn) flag = .false.
         do i=1,ne
            if(b0(i) > 1d-6 .and. abs(b0(i)-b(i))>b0max*1d-6)      flag = .false.
         end do
         if(abs(DlogT)>1d-5 .and. abs(E-b(ne+1)*Ru*T)/abs(E) >eps) flag = .false.
         if(counter>500) flag = .true.
         if(flag) exit

         counter=counter+1
      end do

      if(counter>500) print *,"not converted. at calc_therm"
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

!!for CH4/Air!!!!!!!!!!!!!!!
!subroutine calc_initial_prop(xi,T, n,E,MWini,b0)!CH4/Air{{{
!   use chem
!   use func_therm
!   implicit none
!   double precision,intent(in) ::xi
!   double precision,intent(in) ::T
!   double precision,intent(out)::n(ns)
!   double precision,intent(out)::E
!   double precision,intent(out)::MWini
!   double precision,intent(out)::b0(1:ne)
!
!   integer,parameter::iCH4 =7
!   integer,parameter::iN2  =144
!   integer,parameter::iO2  =157
!   double precision,parameter::etaO2=0.215d0
!   double precision,parameter::etaN2=1d0-etaO2
!
!   double precision logT
!   integer i
!   integer sCH4,sO2,sN2
!
!   !set mass fraction
!   !g-mol/kg
!   n=0d0
!   n(iCH4) =     xi /MW(iCH4)*1d3 !C3H6
!   n(iO2 ) =(1d0-xi)/(MW(iO2)*etaO2+MW(iN2)*etaN2)*etaO2*1d3 !O2
!   n(iN2 ) =(1d0-xi)/(MW(iO2)*etaO2+MW(iN2)*etaN2)*etaN2*1d3 !O2
!
!   !E
!   logT=log(T)
!   call check_section_number(iCH4,T,sCH4)
!   call check_section_number( iO2,T, sO2)
!   call check_section_number( iN2,T, sN2)
!   E=n(iCH4)*ert(iCH4,T,logT,sCH4)&
!    +n( iO2)*ert( iO2,T,logT, sO2)&
!    +n( iN2)*ert( iN2,T,logT, sN2)
!   E=E*Ru*T
!
!   !MWini
!   MWini=1d3/(n(iCH4)+n(iO2)+n(iN2))
!
!   !b0
!   do i=1,ne
!      b0(i)=Ac(i,iCH4)*n(iCH4)&
!           +Ac(i, iO2)*n( iO2)&
!           +Ac(i, iN2)*n( iN2)
!   end do
!
!   !modification for mathematical stability
!   n=n+initial_eps
!end subroutine calc_initial_prop!}}}
!subroutine calc_therm_flame_sheet_boundary(xi,T, E,MWave,kappa,mu)!CH4/Air{{{
!   use chem
!   use func_therm
!   double precision,intent(in)::xi
!   double precision,intent(in)::T
!   double precision,intent(out)::E
!   double precision,intent(out)::MWave
!   double precision,intent(out)::kappa
!   double precision,intent(out)::mu
!
!   double precision nCH4,nO2,nN2
!   double precision nCO2 ,nH2O
!
!   integer,parameter::iCH4 =7
!   integer,parameter::iN2  =144
!   integer,parameter::iO2  =157
!   integer,parameter::iCO2 =13
!   integer,parameter::iH2O =132
!
!   integer,parameter::itCH4 =2
!   integer,parameter::itN2  =19
!   integer,parameter::itO2  =24
!   integer,parameter::itCO2 =5
!   integer,parameter::itH2O =14
!
!   double precision,parameter::etaO2=0.215d0
!   double precision,parameter::etaN2=1d0-etaO2
!
!   integer,parameter:: ind(5) = (/ iCH4, iN2, iO2, iCO2, iH2O/)
!   integer,parameter::indt(5) = (/itCH4,itN2,itO2,itCO2,itH2O/)
!   integer            secN(5)
!   double precision   molN(5)
!   double precision    muN(5)
!
!   double precision MWf,MWo
!   double precision xi_st
!
!   integer sect
!   double precision cprnow,cvrnow,ERTnow,logT,sum_n,phi
!
!   integer i,j
!
!   !set MW
!   MWf=1d-3* MW(iCH4)
!   MWo=1d-3*(MW(iO2)*etaO2+MW(iN2)*etaN2)
!
!   !set xi_st
!   xi_st=MWf*etaO2/(MWf*etaO2+2d0*MWo)
!
!   !set n
!   nN2  = (1d0-xi)/MWo*etaN2
!   if(xi > xi_st) then !fuel rich
!      nO2  = 0d0
!      nCH4 = (xi*(2d0*MWo+MWf*etaO2)-MWf*etaO2)/(2d0*MWf*MWo)
!      nCO2 = 0.5d0*(1d0-xi)/MWo*etaO2
!   else !oxidizer rich
!      nO2  = (MWf*etaO2-xi*(2d0*MWo+MWf*etaO2))/(MWf*MWo)
!      nCH4 = 0d0
!      nCO2 = xi/MWf
!   end if
!   nH2O = 2d0*nCO2
!
!   molN(1)=nCH4
!   molN(2)=nN2
!   molN(3)=nO2
!   molN(4)=nCO2
!   molN(5)=nH2O
!
!   !calc T
!   logT=log(T)
!
!   ERTnow=0d0
!   cvrnow=0d0
!   do i=1,5
!      call check_section_number(ind(i),T,secN(i))
!      ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
!      cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
!   end do
!
!   sum_n=nCH4+nO2+nN2+nCO2+nH2O
!   MWave=1d3/sum_n
!
!   !calc kappa
!   cprnow= cvrnow+sum_n
!   kappa=cprnow/cvrnow
!   E=ERTnow*Ru*T
!
!   !calc mu
!   do i=1,5
!      call check_section_number_trans(indt(i),T,sect)
!      muN(i) = calc_mu(indt(i),T,logT,sect)
!   end do
!
!   mu = 0d0
!   do i=1,5
!      denom = 0d0
!      do j=1,5
!         phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MW(ind(j)) / MW(ind(i)) )**0.25d0 )**2 &
!                      * sqrt(2d0 *MW(ind(j))/(MW(ind(i))+MW(ind(j))))
!         denom = denom + molN(j) * phi
!      end do
!      mu = mu + molN(i) * muN(i) / denom
!   end do
!   mu = mu * 1d-7
!end subroutine calc_therm_flame_sheet_boundary!}}}
!subroutine calc_therm_flame_sheet_inert(xi,E, T, MWave,kappa,mu)!CH4/Air{{{
!   use chem
!   use func_therm
!   double precision,intent(in)::xi
!   double precision,intent(in)::E
!   double precision,intent(inout)::T
!   double precision,intent(out)::MWave
!   double precision,intent(out)::kappa
!   double precision,intent(out)::mu
!
!   double precision nCH4,nO2,nN2
!
!   integer,parameter::iCH4 =7
!   integer,parameter::iN2  =144
!   integer,parameter::iO2  =157
!
!   integer,parameter::itCH4 =2
!   integer,parameter::itN2  =19
!   integer,parameter::itO2  =24
!
!   double precision,parameter::etaO2=0.215d0
!   double precision,parameter::etaN2=1d0-etaO2
!
!   integer,parameter:: ind(3) = (/ iCH4, iN2, iO2/)
!   integer,parameter::indt(3) = (/itCH4,itN2,itO2/)
!   integer            secN(3)
!   double precision   molN(3)
!   double precision    muN(3)
!
!   double precision MWf,MWo
!
!   double precision T_old
!
!   integer sect
!   double precision cprnow,cvrnow,ERTnow,logT,sum_n,phi
!
!   integer,parameter::max_step = 200
!   integer i,j,k
!
!   !set MW
!   MWf=1d-3* MW(iCH4)
!   MWo=1d-3*(MW(iO2)*etaO2+MW(iN2)*etaN2)
!
!   !set n
!   nCH4 =      xi /MWf
!   nO2  = (1d0-xi)/MWo*etaO2
!   nN2  = (1d0-xi)/MWo*etaN2
!
!   molN(1)=nCH4
!   molN(2)=nN2
!   molN(3)=nO2
!
!   !calc T
!   logT=log(T)
!
!   ERTnow=0d0
!   cvrnow=0d0
!   do i=1,3
!      call check_section_number(ind(i),T,secN(i))
!      ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
!      cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
!   end do
!
!   k=0
!   T_old = T*0.5d0
!   do while(abs(E/Ru-ERTnow*T)>cvrnow*1d-10 .and. abs(T-T_old)>1d-8.and. k<max_step)
!      T_old =T
!      T     =T+omega*(E/Ru-ERTnow*T)/cvrnow
!      logT=log(T)
!
!      ERTnow=0d0
!      cvrnow=0d0
!      do i=1,3
!         call check_section_number(ind(i),T,secN(i))
!         ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
!         cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
!      end do
!      k=k+1
!   end do
!
!   if(k >= max_step) then
!      print *, "Not Converted at calc_therm_flame_sheet_inert",T
!      call exit(1)
!   end if
!
!
!   sum_n=nCH4+nO2+nN2
!   MWave=1d3/sum_n
!
!   !calc kappa
!   cprnow= cvrnow+sum_n
!   kappa=cprnow/cvrnow
!
!   !calc mu
!   do i=1,3
!      call check_section_number_trans(indt(i),T,sect)
!      muN(i) = calc_mu(indt(i),T,logT,sect)
!   end do
!
!   mu = 0d0
!   do i=1,3
!      denom = 0d0
!      do j=1,3
!         phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MW(ind(j)) / MW(ind(i)) )**0.25d0 )**2 &
!                      * sqrt(2d0 *MW(ind(j))/(MW(ind(i))+MW(ind(j))))
!         denom = denom + molN(j) * phi
!      end do
!      mu = mu + molN(i) * muN(i) / denom
!   end do
!   mu = mu * 1d-7
!end subroutine calc_therm_flame_sheet_inert!}}}
!subroutine calc_therm_flame_sheet_inert_boundary(xi,T, E,MWave,kappa,mu)!CH4/Air{{{
!   use chem
!   use func_therm
!   double precision,intent(in)::xi
!   double precision,intent(in)::T
!   double precision,intent(out)::E
!   double precision,intent(out)::MWave
!   double precision,intent(out)::kappa
!   double precision,intent(out)::mu
!
!   double precision nCH4,nO2,nN2
!
!   integer,parameter::iCH4 =7
!   integer,parameter::iN2  =144
!   integer,parameter::iO2  =157
!
!   integer,parameter::itCH4 =2
!   integer,parameter::itN2  =19
!   integer,parameter::itO2  =24
!
!   double precision,parameter::etaO2=0.215d0
!   double precision,parameter::etaN2=1d0-etaO2
!
!   integer,parameter:: ind(3) = (/ iCH4, iN2, iO2/)
!   integer,parameter::indt(3) = (/itCH4,itN2,itO2/)
!   integer            secN(3)
!   double precision   molN(3)
!   double precision    muN(3)
!
!   double precision MWf,MWo
!   double precision xi_st
!
!   integer sect
!   double precision cprnow,cvrnow,ERTnow,logT,sum_n,phi
!
!   integer i,j
!
!   !set MW
!   MWf=1d-3* MW(iCH4)
!   MWo=1d-3*(MW(iO2)*etaO2+MW(iN2)*etaN2)
!
!   !set n
!   nCH4 =      xi /MWf
!   nO2  = (1d0-xi)/MWo*etaO2
!   nN2  = (1d0-xi)/MWo*etaN2
!
!   molN(1)=nCH4
!   molN(2)=nN2
!   molN(3)=nO2
!
!   !calc T
!   logT=log(T)
!
!   ERTnow=0d0
!   cvrnow=0d0
!   do i=1,3
!      call check_section_number(ind(i),T,secN(i))
!      ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
!      cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
!   end do
!
!   sum_n=nCH4+nO2+nN2
!   MWave=1d3/sum_n
!
!   !calc kappa
!   cprnow= cvrnow+sum_n
!   kappa=cprnow/cvrnow
!   E=ERTnow*Ru*T
!
!   !calc mu
!   do i=1,3
!      call check_section_number_trans(indt(i),T,sect)
!      muN(i) = calc_mu(indt(i),T,logT,sect)
!   end do
!
!   mu = 0d0
!   do i=1,3
!      denom = 0d0
!      do j=1,3
!         phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MW(ind(j)) / MW(ind(i)) )**0.25d0 )**2 &
!                      * sqrt(2d0 *MW(ind(j))/(MW(ind(i))+MW(ind(j))))
!         denom = denom + molN(j) * phi
!      end do
!      mu = mu + molN(i) * muN(i) / denom
!   end do
!   mu = mu * 1d-7
!end subroutine calc_therm_flame_sheet_inert_boundary!}}}
!!end for CH4/Air!!!!!!!!!!!

