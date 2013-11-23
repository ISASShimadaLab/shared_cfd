! read data files
subroutine read_elmspc_to_use!{{{
   use mod_mpi
   use chem
   implicit none
   integer i

   if(myid .eq. 0) then
      open(22,file='elements_to_use.inp')
      do i=1,ne
         read(22,'(a2)') SYM_ELM(i)
      end do
      close(22)

      open(22,file='species_to_use.inp')
      do i=1,ns
         read(22,'(a18)') SYM_SPC(i)
      end do
      close(22)
   end if
   call MPI_Bcast( SYM_SPC,   ns*18,        MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast( SYM_ELM,    ne*2,        MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
end subroutine read_elmspc_to_use!}}}
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
         if(ind .eq. 0) stop "Bad format at control_chem.inp while reading Oxygen Compositions"
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

   !!! calculate other properties for oxygen and fuel
   call set_OF_ratio
   print *,of,(32d0*2d0)/16d0
   call calc_ini(po,To, no, Eo)
   call calc_ini(pf,Tf, nf, Ef)
   print *,po,To,Eo
   print *,pf,Tf,Ef
   !call rho_rel2abs(po,To, vrhoo, rhoo,Eo)
   !call calc_T(vrhoo,rhoo,Eo, To, kappao,MWo,DHit,vhio,muo)

   !call rho_rel2abs(pf,Tf, vrhof, rhof,Ef)
   !call calc_T(vrhof,rhof,Ef, Tf, kappaf,MWf,DHit,vhif,muf)
end subroutine read_fo_composition!}}}

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

!initialization
subroutine set_OF_ratio!{{{
   use chem
   use chem_var
   implicit none
   double precision sumo,sumf
   integer j

   sumo=0d0
   sumf=0d0
   do j=1,ns
      sumo=sumo+MW(j)*no(j)
      sumf=sumf+MW(j)*nf(j)
   end do
   of = sumo/sumf
end subroutine set_OF_ratio!}}}
subroutine calc_ini(p,T, n, E)!{{{
   use chem
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
      tmp=0d0
      do k=1,9
         tmp=tmp+vThrt(k)*co(k,nsc,j)
      end do
      print *,j,tmp*Ru*T,n(j)
      E=tmp*n(j)
   end do
   E=E*Ru*T

   n=n+initial_eps
end subroutine calc_ini!}}}

!main
subroutine calc_therm(rho,xi,E, T,n, MWave,kappa,mu)!{{{
   use const_chem
   use func_therm
   use chem
   use chem_var
   implicit none
   double precision,intent(in)   ::rho
   double precision,intent(in)   ::xi
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
   !double precision maxDlogn,logsn
   !double precision lambda,lambda1,lambda2

   integer ne_now

   integer sect
   double precision,dimension(ns)::muN
   double precision phi,denom,MWj,MWi

   integer i,j,k,counter
   logical flag,LUflag

   double precision,dimension(9)::vTmurt,vTcpr,vThrt


   ne_now = 4
   if(xi < 2.505d-6) ne_now=2

   logrhoRu = log(rho*Ru *1d-5)

   !calc b0 of elements
   do i=1,ne_now
      b0(i)= xi*b0f(i) + (1d0-xi)*b0o(i)
   end do
   b0max = maxval(b0(1:ne_now),1)

   sn=sum(n(1:ns),1)
   !logsn = log(sn)
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
         do i=1,ne_now
            b(i) = b(i) + Ac(i,j)*n(j)
         end do
         b(ne_now+1)=b(ne_now+1)+vert(j)*n(j)
      end if
   end do

   !check convergence
   flag = .true.
   do i=1,ne_now
      if(b0(i) > 1d-6 .and. abs(b0(i)-b(i))>b0max*1d-6) flag = .false.
   end do
   if(abs(E-b(ne_now+1)*Ru*T)/abs(E) >eps) flag = .false.

   if(.not. flag) then
      counter=1
      do
         !set b0(ne+1)
         b0(ne_now+1)=E/(Ru*T)

         !set b and vmurtn
         b(1:ne_now)=0d0
         bd=b0
         bd(ne_now+1)=bd(ne_now+1)-b(ne_now+1)
         Ad(1:ne_now+1,1:ne_now+1)=0d0
         call calc_vTmurt(T,logT,vTmurt)
         call calc_vTcpr(T,vTcpr)
         do k=1,ns
            vmurt(k)=dot_product(co(:,sec_num(k),k),vTmurt)+logrhoRu+logT+logn(k)

            if(n(k)>tinyn) then
               vmurtn=vmurt(k)*n(k)
               do i=1,ne_now
                  do j=i,ne_now
                     Ad(i,j)=Ad(i,j)+Ac(i,k)*Ac(j,k)*n(k)
                  end do
                  Ad(i,ne_now+1)=Ad(i,ne_now+1)+Ac(i,k)*vert(k)*n(k)
                  b(i)  = b(i)  + Ac(i,k)* n(k)
                  bd(i) = bd(i) + Ac(i,k)*(vmurtn-n(k))
               end do

               Ad(ne_now+1,ne_now+1)= Ad(ne_now+1,ne_now+1)+ n(k) *(vert(k)**2+dot_product(co(:,sec_num(k),k),vTcpr)-1d0)
               bd(ne_now+1)     = bd(ne_now+1)     + vmurtn*vert(k)
            end if
         end do

         do i=1,ne_now+1
            do j=i,ne_now+1
               Ad(j,i)=Ad(i,j)
            end do
         end do

         !print *,"Ad="
         !print '(3es9.1)',Ad(1:3,1)
         !print '(3es9.1)',Ad(1:3,2)
         !print '(3es9.1)',Ad(1:3,3)
 
         !calc vpi
         call LU(Ad,bd,vpi,ne_now+1,LUflag)
         if(LUflag) then
            !write(100,*) xi
            !print *,"xi:",xi
            !print *,"tiny n:",tinyn
            !print *,"n:"
            !do i=1,ns
            !   if(n(i)>=tinyn) print '(a15,x,es9.1)',species_name(i),n(i)
            !end do

            !set new n
            n=n+initial_eps

            !set new n
            b(ne_now+1)=0d0
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
                   b(ne_now+1) = b(ne_now+1)+vert(j)*n(j)
               end if
            end do
            sn = sum(n(1:ns),1)

            !call exit(0)
            cycle
         end if

         !calc Delta
         DlogT=vpi(ne_now+1)
         vpi(ne_now+1:ne)= -1d300
         !maxDlogn=0d0
         !lambda2=1d0
         do j=1,ns
            tmp=0d0
            do i=1,ne
               tmp=tmp+Ac(i,j)*vpi(i)
            end do
            Dlogn(j)= -vmurt(j)+tmp+vert(j)*DlogT

            !if(maxDlogn < abs(Dlogn(j))) maxDlogn = abs(Dlogn(j))
            !if(n(j) .eq. 0d0 .and. Dlogn(j) > 0d0) lambda2 = min(lambda2,(-logn(j)+logsn-9.21034d0)/Dlogn(j))
         end do

         !!calc lambda
         !lambda1 = 2d0/max(5d0*abs(DlogT),maxDlogn)
         !lambda  = min(1d0,lambda1,lambda2)

         !set new T
         T=T*(1d0+omega*DlogT)
         if(T < 0d0) then
            print *,"Negative Temperature. T=",T
            call exit(1)
         end if
         logT=log(T)
         do j=1,ns
            call check_section_number(j,T,sec_num(j))
         end do

         !set new n
         b(ne_now+1)=0d0
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
                b(ne_now+1) = b(ne_now+1)+vert(j)*n(j)
                tmp = abs(Dlogn(j))*n(j)
                if(Dnmax < tmp) Dnmax = tmp
            end if
         end do
         sn = sum(n(1:ns),1)
         !logsn = log(sn)

         !check convergence
         flag = .true.
         if(abs(Dnmax)>1d-5*sn) flag = .false.
         do i=1,ne_now
            if(b0(i) > 1d-6 .and. abs(b0(i)-b(i))>b0max*1d-6)      flag = .false.
         end do
         if(abs(DlogT)>1d-5 .and. abs(E-b(ne_now+1)*Ru*T)/abs(E) >eps) flag = .false.
         if(counter>500) flag = .true.
         if(flag) exit

         counter=counter+1
      end do

      !print '(a,i5)',"counter=",counter
      if(counter>500) then
         print *,"not converted. at calc_therm"
      end if
   end if

   !call calc_kappa(kappa)
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
   do i=1,nt
      call check_section_number_trans(i,T,sect)
      muN(i) = calc_mu(i,T,logT,sect)
   end do

   mu = 0d0
   do i=1,nt
      denom = 0d0
      MWi=MW(tr2th(i))
      do j=1,nt
         MWj=MW(tr2th(j))
         phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MWj / MWi)**0.25d0 )**2 &
                      * sqrt(2d0 *MWj/(MWi+MWj))
         denom = denom + n(tr2th(j)) * phi
      end do
      mu = mu + n(tr2th(i)) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine calc_therm!}}}

!!for C3H6/O2!!!!!!!!!!!!!!
!subroutine calc_initial_prop(xi,T, n,E,MWini)!C3H6/O2{{{
!   use chem
!   use func_therm
!   implicit none
!   double precision,intent(in) ::xi
!   double precision,intent(in) ::T
!   double precision,intent(out)::n(ns)
!   double precision,intent(out)::E
!   double precision,intent(out)::MWini
!
!   integer i
!
!   !set mass fraction
!   !g-mol/kg
!   n=0d0
!   n(40) =(    xi)/MW(40) *1d3 !C3H6
!   n(120)=(1d0-xi)/MW(120)*1d3 !O2
!
!   !E
!   E=n(40)*ert(40,T)+n(120)*ert(120,T)
!   E=E*Ru*T
!
!   !MWini
!   MWini=1d3/sum(n(1:ns),1)
!
!   !modification for mathematical stability
!   n=n+initial_eps
!end subroutine calc_initial_prop!}}}
!
!subroutine calc_therm_flame_sheet(xi,E, T, MWave,kappa)!C3H6/O2{{{
!   use chem
!   use func_therm
!   double precision,intent(in)::xi
!   double precision,intent(in)::E
!   double precision,intent(inout)::T
!   double precision,intent(out)::MWave
!   double precision,intent(out)::kappa
!
!   double precision nC3H6,nO2
!   double precision nCO2 ,nH2O
!
!   integer,parameter::iC3H6=48
!   integer,parameter::iO2  =157
!   integer,parameter::iCO2 =13
!   integer,parameter::iH2O =132
!
!   double precision MWf,MWo
!   double precision xi_st
!
!   double precision cprnow,cvrnow,ERTnow
!
!   !set MW
!   MWf=MW(iC3H6)
!   MWo=MW(iO2)
!
!   !set xi_st
!   xi_st=2d0*MWf/(2d0*MWf+9d0*MWo)
!
!   !set n
!   if(xi > xi_st) then !fuel rich
!      nC3H6= 1d3*(-2d0*MWf+(2d0*MWf+9d0*MWo)*xi)/(9d0*MWo*MWf)
!      nO2  = 1d3*0d0
!      nCO2 = 1d3*2d0/3d0*(1d0-xi)/MWo
!      nH2O = 1d3*nCO2
!   else !oxydizer rich
!      nC3H6= 1d3*0d0
!      nO2  = 1d3*(2d0*MWf-(2d0*MWf+9d0*MWo)*xi)/(2d0*MWo*MWf)
!      nCO2 = 1d3*3d0*xi/MWf
!      nH2O = 1d3*nCO2
!   end if
!
!   !calc T
!   ERTnow= ert(iC3H6,T)*nC3H6 &
!         + ert(iO2  ,T)*nO2   &
!         + ert(iCO2 ,T)*nCO2  &
!         + ert(iH2O ,T)*nH2O
!
!   cvrnow= cvr(iC3H6,T)*nC3H6 &
!         + cvr(iO2  ,T)*nO2   &
!         + cvr(iCO2 ,T)*nCO2  &
!         + cvr(iH2O ,T)*nH2O
!
!   do while(abs(E-ERTnow*Ru*T)>abs(E)*eps)
!      T     =T+omega*(E/Ru-ERTnow*T)/cvrnow
!
!      ERTnow= ert(iC3H6,T)*nC3H6 &
!            + ert(iO2  ,T)*nO2   &
!            + ert(iCO2 ,T)*nCO2  &
!            + ert(iH2O ,T)*nH2O
!
!      cvrnow= cvr(iC3H6,T)*nC3H6 &
!            + cvr(iO2  ,T)*nO2   &
!            + cvr(iCO2 ,T)*nCO2  &
!            + cvr(iH2O ,T)*nH2O
!   end do
!
!   MWave=1d3/(nC3H6+nO2+nCO2+nH2O)
!
!   !calc kappa
!   cprnow= cpr(iC3H6,T)*nC3H6 &
!         + cpr(iO2  ,T)*nO2   &
!         + cpr(iCO2 ,T)*nCO2  &
!         + cpr(iH2O ,T)*nH2O
!   kappa=cprnow/cvrnow
!end subroutine calc_therm_flame_sheet!}}}
!!end for C3H6/O2!!!!!!!!!!

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
!
!subroutine calc_therm_flame_sheet(xi,E, T, MWave,kappa,mu)!CH4/Air{{{
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
!   double precision T_old
!
!   integer sect
!   double precision cprnow,cvrnow,ERTnow,logT,sum_n,phi
!
!   integer i,j,k
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
!   else !oxydizer rich
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
!   k=0
!   T_old = T*0.5d0
!   do while(abs(E/Ru-ERTnow*T)>cvrnow*1d-10 .and. abs(T-T_old)>1d-8 .and. k<100)
!      T_old =T
!      T     =T+omega*(E/Ru-ERTnow*T)/cvrnow
!      logT=log(T)
!
!      ERTnow=0d0
!      cvrnow=0d0
!      do i=1,5
!         call check_section_number(ind(i),T,secN(i))
!         ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
!         cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
!      end do
!      k=k+1
!   end do
!
!   if(k >= 100) then
!      print *, "Not Converted at calc_therm_flame_sheet"
!      !call exit(1)
!   end if
!
!
!   sum_n=nCH4+nO2+nN2+nCO2+nH2O
!   MWave=1d3/sum_n
!
!   !calc kappa
!   cprnow= cvrnow+sum_n
!   kappa=cprnow/cvrnow
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
!end subroutine calc_therm_flame_sheet!}}}
!
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
!   else !oxydizer rich
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
!
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
!
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

