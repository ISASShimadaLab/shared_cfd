!read control.inp
subroutine read_fo_composition!{{{
   use chem
   use mod_mpi
   use chem_var
   implicit none
   character*100 buf,buf2,bufline(100)
   double precision amount
   integer ind

   integer i

   if(myid .eq. 0) then
      no =1d-20
      nf =1d-20

      open(8,file="control_chem.inp")
      read(8,'()')
      read(8,'(25x,es15.7)') po
      read(8,'(25x,es15.7)') To
      read(8,'()')
      !count
      numo=0
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit
         numo = numo + 1
         bufline(numo)=buf
      end do
      print *,po,To,numo
      allocate(SYM_OXID(numo),COMP_OXID(numo))
      !record
      do i=1,numo
         buf=adjustl(trim(bufline(i))) !delete front and back spaces
         ind=index(buf,' ')
         if(ind .eq. 0) stop "Bad format at control_chem.inp while reading Oxygen Compositions"
         read(buf(ind+1:),*) amount
         SYM_OXID(i)   = buf(1:ind-1)
         COMP_OXID(i) = amount
         print *,SYM_OXID(i),COMP_OXID(i)
      end do

      read(8,'(25x,es15.7)') pf
      read(8,'(25x,es15.7)') Tf
      read(8,'()')
      !count
      numf=0
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit
         numf = numf + 1
         bufline(numf)=buf
      end do
      print *,pf,Tf,numf
      allocate(SYM_FUEL(numf),COMP_FUEL(numf))
      !record
      do i=1,numf
         buf=adjustl(trim(bufline(i))) !delete front and back spaces
         ind=index(buf,' ')
         if(ind .eq. 0) stop "Bad format at control_chem.inp while reading Fuel Compositions"
         read(buf(ind+1:),*) amount
         SYM_FUEL(i)   = buf(1:ind-1)
         COMP_FUEL(i) = amount
         print *,SYM_FUEL(i),COMP_FUEL(i)
      end do
      close(8)
   end if

   !!!! MPI COMMUNICATIONS
   !call MPI_Bcast(po,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   !call MPI_Bcast(To,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   !call MPI_Bcast(vrhoo,   ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   !call MPI_Bcast(pf,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   !call MPI_Bcast(Tf,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   !call MPI_Bcast(vrhof,   ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)


   !!!! calculate other properties for oxygen and fuel
   !call rho_rel2abs(po,To, vrhoo, rhoo,Eo)
   !call calc_T(vrhoo,rhoo,Eo, To, kappao,MWo,DHit,vhio,muo)
   !vwo(1:nY) = vrhoo(1:nY)/rhoo

   !call rho_rel2abs(pf,Tf, vrhof, rhof,Ef)
   !call calc_T(vrhof,rhof,Ef, Tf, kappaf,MWf,DHit,vhif,muf)
   !vwf(1:nY) = vrhof(1:nY)/rhof
!contains
!   integer function search_species(str)
!      use chem
!      character(*),intent(in)::str
!      integer i
!   
!      do i=1,
!         if(trim(SYM_SPC(i)) .eq. trim(str)) then
!            search_species=i
!            return
!         end if
!      end do
!      search_species=-1
!   end function search_species
end subroutine read_fo_composition!}}}

subroutine set_therm_data!{{{
   use mod_mpi
   use const_chem
   use chem
   use chem_var
   implicit none
   character*18     sname
   character*300    comment
   character*6      optional_ID_code
   character*2      name_element(5)
   double precision num_element(5)
   integer          phase,num_T_exp
   double precision Texp(8),H0a

   logical to_record
   integer i,j,k,l

   !set elements name to use
   SYM_ELM(1)="O "
   SYM_ELM(2)="N "
   SYM_ELM(3)="C "
   SYM_ELM(4)="H "

   if(myid .eq. 0) then
      open(10,file='thermo.inp',status='old')
      do
         read(10,'(A)') comment
         if(comment(1:1) .ne. '!') exit
      end do
      read(10,'()')!ranges

      ns=1
      do_species1:do
            read(10,'(A18,A62)') sname, comment
            if(sname(1:3) .eq. 'END') exit

            read(10,'(I2,1x ,A6,1x ,5(A2,F6.2),1x ,I1,F13.7,F15.3)') &
               num_sctn(ns), optional_ID_code, (name_element(k), num_element(k), k=1,5),&
               phase, mw(ns), DH0(ns)

            if(existInReactant(sname)) then
               outer1:do k=1,5 
                  if(name_element(k) .eq. "  ") exit
                  do i=1,ne
                     if(SYM_ELM(i) .eq. name_element(k)) then
                        cycle outer1
                     end if
                  end do
               end do outer1
            end if

            do l=1, num_sctn(ns)*3
               read(10,'()')
            end do
      end do do_species1

      close(10)
      call exit(0)
      !!! END OF FETCH ATOMS


      !!! START RECORD SPECIES
      open(10,file='thermo.inp',status='old')

      !!! header
      do
         read(10,'(A)') comment
         if(comment(1:1) .ne. '!') exit
      end do
      read(10,'()')!ranges
      !!! end of header

      ns=1
      to_record = .true.
      do_species:do
            !first line --- name, comment
            read(10,'(A18,A62)') sname, comment
            if(sname(1:3) .eq. 'END') exit
            SYM_SPC(ns)=sname

            !second line --- other property
            !read other informations
            read(10,'(I2,1x ,A6,1x ,5(A2,F6.2),1x ,I1,F13.7,F15.3)') &
               num_sctn(ns), optional_ID_code, (name_element(k), num_element(k), k=1,5),&
               phase, mw(ns), DH0(ns)

            !for condensed phase
            if(phase .ne. 0) to_record = .false.

            !set number of elements
            Ac(1:ne,ns)=0d0
            outer:do k=1,5 
               if(name_element(k) .eq. "  ") exit

               do i=1,ne
                  if(SYM_ELM(i) .eq. name_element(k)) then
                     Ac(i,ns)=num_element(k)
                     cycle outer
                  end if
               end do
               to_record = .false.
            end do outer

            !other line ---coefficients of thermodynamical functions
            do l=1, num_sctn(ns)
               read(10,'(1x ,F10.3,1x ,F10.3, I1, 8F5.1, 2x ,F15.3)') &
                  Trange(1:2,l,ns), num_T_exp, Texp(1:8), H0a
               read(10,'(5D16.9)') co(1:5,l,ns)
               read(10,'(2D16.8,16x,2D16.8)') co(6:9,l,ns)

               if(num_T_exp .ne. 7) then
                  print *,"Error: Odd function at species:",SYM_SPC(j)
                  stop
               end if
            end do

            if(to_record) ns=ns+1
            to_record=.true.
      end do do_species
      ns=ns-1

      print *,"the number of species to use:",ns

      close(10)
   end if

   call MPI_Bcast(      ns,          1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(num_sctn,         ns,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(      MW,         ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(     DH0,         ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  Trange, max_ns*6*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(      co, max_ns*6*9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(      Ac,  ne*max_ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast( SYM_SPC,  max_ns*18,        MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast( SYM_ELM,      ne*18,        MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
contains
   logical function existInReactant(str)
      implicit none
      character(*),intent(in)::str
      integer i
      existInReactant = .true.
      do i=1,numo
         if(trim(SYM_OXID(i)) .eq. trim(str)) return
      end do
      do i=1,numf
         if(trim(SYM_FUEL(i)) .eq. trim(str)) return
      end do
      existInReactant = .false.
   end function existInReactant
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

   call MPI_Bcast(            nt,          1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(num_sctn_trans,     max_ns,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(         trans, 4*3*max_ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(  Trange_trans, 2*3*max_ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(         tr2th,     max_ns,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end subroutine set_trans_data!}}}

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
   double precision,intent(inout)::n(max_ns)
   double precision,intent(out)  ::MWave
   double precision,intent(out)  ::kappa
   double precision,intent(out)  ::mu

   double precision,dimension(ne+1)::b0,b,bd,vpi
   double precision,dimension(max_ns)::vmurt,vert,Dlogn,logn
   integer,dimension(max_ns)::sec_num
   double precision,dimension(ne+1,ne+1)::Ad
   double precision tmp,logT,sn,vmurtn,DlogT,b0max,Dnmax,logrhoRu,tinyn,logtinyn
   double precision,parameter::tinyratio=1d-3
   double precision,parameter::logtinyratio=-2.3026d0*3d0 !log(tinyratio)
   !double precision maxDlogn,logsn
   !double precision lambda,lambda1,lambda2

   integer ne_now

   integer sect
   double precision,dimension(max_ns)::muN
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

!subroutine calc_kappa_s(kappa)!{{{
!   use const_chem
!   use func_therm
!   use chem
!   implicit none
!   double precision,intent(out)::kappa
!
!   double precision      Ad(ne+2,  ne+2)
!   double precision Ad_temp(ne+1,ne+1)
!   double precision bd(ne+1)
!   double precision dpi(ne+2)
!   double precision temp 
!   integer i,j,k
!
!   double precision cpr_bar,cvr_bar,dLVdLP,dLVdLT
!
!   do j=1,ns
!      A(ne+1,j)  =1d0
!      A(ne+2,j)  =hrt(j)
!   end do
!
!   do i=1,ne+2
!      do j=i,ne+2
!         temp=0d0
!         do k=1,ns
!            temp=temp+A(i,k)*A(j,k)*n(k)
!         end do
!         Ad(i,j)=temp
!         Ad(j,i)=temp
!      end do
!   end do
!
!   temp=0d0
!   do k=1,ns
!      temp=temp+n(k)*cpr(k)
!   end do
!   Ad(ne+2,ne+2)=Ad(ne+2,ne+2)+temp
!
!
!   Ad_temp=Ad(1:ne+1,1:ne+1)
!   Ad_temp(ne+1,ne+1)=0d0
!   bd=-Ad(ne+2,1:ne+1)
!   call LU(Ad_temp,bd,dpi(1:ne+1),ne+1)
!   dLVdLT=1d0+dpi(ne+1)
!   !print *,"dLVdLT=",dLVdLT
!
!   dpi(ne+2)=1d0
!   cpr_bar=0d0
!   do k=1,ne+2
!      cpr_bar=cpr_bar+Ad(ne+2,k)*dpi(k)
!   end do
!   !print *,"cp_bar=",cpr_bar*Ru*1d-3
!
!   Ad_temp=Ad(1:ne+1,1:ne+1)
!   Ad_temp(ne+1,ne+1)=0d0
!   bd=Ad(ne+1,1:ne+1)
!   call LU(Ad_temp,bd,dpi(1:ne+1),ne+1)
!   dLVdLP=-1d0+dpi(ne+1)
!   !print *,"dLVdLP=",dLVdLP
!
!   kappa=1d0+sum(n(1:ns),1)*dLVdLT**2/dLVdLP/cpr_bar  !eq. (2.70)  kappa=Cv/Cp
!   kappa=1d0/kappa/(-dLVdLP)                          !eq. (2.73)
!end subroutine calc_kappa_s!}}}

integer function search_species(str)!{{{
   use chem
   character(*),intent(in)::str
   integer i

   do i=1,ns
      if(trim(SYM_ELM(i)) .eq. trim(str)) then
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

!!for C3H6/O2!!!!!!!!!!!!!!
!subroutine calc_initial_prop(xi,T, n,E,MWini)!C3H6/O2{{{
!   use chem
!   use func_therm
!   implicit none
!   double precision,intent(in) ::xi
!   double precision,intent(in) ::T
!   double precision,intent(out)::n(max_ns)
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

!for CH4/Air!!!!!!!!!!!!!!!
subroutine calc_initial_prop(xi,T, n,E,MWini,b0)!CH4/Air{{{
   use chem
   use func_therm
   implicit none
   double precision,intent(in) ::xi
   double precision,intent(in) ::T
   double precision,intent(out)::n(max_ns)
   double precision,intent(out)::E
   double precision,intent(out)::MWini
   double precision,intent(out)::b0(1:ne)

   integer,parameter::iCH4 =7
   integer,parameter::iN2  =144
   integer,parameter::iO2  =157
   double precision,parameter::etaO2=0.215d0
   double precision,parameter::etaN2=1d0-etaO2

   double precision logT
   integer i
   integer sCH4,sO2,sN2

   !set mass fraction
   !g-mol/kg
   n=0d0
   n(iCH4) =     xi /MW(iCH4)*1d3 !C3H6
   n(iO2 ) =(1d0-xi)/(MW(iO2)*etaO2+MW(iN2)*etaN2)*etaO2*1d3 !O2
   n(iN2 ) =(1d0-xi)/(MW(iO2)*etaO2+MW(iN2)*etaN2)*etaN2*1d3 !O2

   !E
   logT=log(T)
   call check_section_number(iCH4,T,sCH4)
   call check_section_number( iO2,T, sO2)
   call check_section_number( iN2,T, sN2)
   E=n(iCH4)*ert(iCH4,T,logT,sCH4)&
    +n( iO2)*ert( iO2,T,logT, sO2)&
    +n( iN2)*ert( iN2,T,logT, sN2)
   E=E*Ru*T

   !MWini
   MWini=1d3/(n(iCH4)+n(iO2)+n(iN2))

   !b0
   do i=1,ne
      b0(i)=Ac(i,iCH4)*n(iCH4)&
           +Ac(i, iO2)*n( iO2)&
           +Ac(i, iN2)*n( iN2)
   end do

   !modification for mathematical stability
   n=n+initial_eps
end subroutine calc_initial_prop!}}}

subroutine calc_therm_flame_sheet(xi,E, T, MWave,kappa,mu)!CH4/Air{{{
   use chem
   use func_therm
   double precision,intent(in)::xi
   double precision,intent(in)::E
   double precision,intent(inout)::T
   double precision,intent(out)::MWave
   double precision,intent(out)::kappa
   double precision,intent(out)::mu

   double precision nCH4,nO2,nN2
   double precision nCO2 ,nH2O

   integer,parameter::iCH4 =7
   integer,parameter::iN2  =144
   integer,parameter::iO2  =157
   integer,parameter::iCO2 =13
   integer,parameter::iH2O =132

   integer,parameter::itCH4 =2
   integer,parameter::itN2  =19
   integer,parameter::itO2  =24
   integer,parameter::itCO2 =5
   integer,parameter::itH2O =14

   double precision,parameter::etaO2=0.215d0
   double precision,parameter::etaN2=1d0-etaO2

   integer,parameter:: ind(5) = (/ iCH4, iN2, iO2, iCO2, iH2O/)
   integer,parameter::indt(5) = (/itCH4,itN2,itO2,itCO2,itH2O/)
   integer            secN(5)
   double precision   molN(5)
   double precision    muN(5)

   double precision MWf,MWo
   double precision xi_st

   double precision T_old

   integer sect
   double precision cprnow,cvrnow,ERTnow,logT,sum_n,phi

   integer i,j,k

   !set MW
   MWf=1d-3* MW(iCH4)
   MWo=1d-3*(MW(iO2)*etaO2+MW(iN2)*etaN2)

   !set xi_st
   xi_st=MWf*etaO2/(MWf*etaO2+2d0*MWo)

   !set n
   nN2  = (1d0-xi)/MWo*etaN2
   if(xi > xi_st) then !fuel rich
      nO2  = 0d0
      nCH4 = (xi*(2d0*MWo+MWf*etaO2)-MWf*etaO2)/(2d0*MWf*MWo)
      nCO2 = 0.5d0*(1d0-xi)/MWo*etaO2
   else !oxydizer rich
      nO2  = (MWf*etaO2-xi*(2d0*MWo+MWf*etaO2))/(MWf*MWo)
      nCH4 = 0d0
      nCO2 = xi/MWf
   end if
   nH2O = 2d0*nCO2

   molN(1)=nCH4
   molN(2)=nN2
   molN(3)=nO2
   molN(4)=nCO2
   molN(5)=nH2O

   !calc T
   logT=log(T)

   ERTnow=0d0
   cvrnow=0d0
   do i=1,5
      call check_section_number(ind(i),T,secN(i))
      ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
      cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
   end do

   k=0
   T_old = T*0.5d0
   do while(abs(E/Ru-ERTnow*T)>cvrnow*1d-10 .and. abs(T-T_old)>1d-8 .and. k<100)
      T_old =T
      T     =T+omega*(E/Ru-ERTnow*T)/cvrnow
      logT=log(T)

      ERTnow=0d0
      cvrnow=0d0
      do i=1,5
         call check_section_number(ind(i),T,secN(i))
         ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
         cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
      end do
      k=k+1
   end do

   if(k >= 100) then
      print *, "Not Converted at calc_therm_flame_sheet"
      !call exit(1)
   end if


   sum_n=nCH4+nO2+nN2+nCO2+nH2O
   MWave=1d3/sum_n

   !calc kappa
   cprnow= cvrnow+sum_n
   kappa=cprnow/cvrnow

   !calc mu
   do i=1,5
      call check_section_number_trans(indt(i),T,sect)
      muN(i) = calc_mu(indt(i),T,logT,sect)
   end do

   mu = 0d0
   do i=1,5
      denom = 0d0
      do j=1,5
         phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MW(ind(j)) / MW(ind(i)) )**0.25d0 )**2 &
                      * sqrt(2d0 *MW(ind(j))/(MW(ind(i))+MW(ind(j))))
         denom = denom + molN(j) * phi
      end do
      mu = mu + molN(i) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine calc_therm_flame_sheet!}}}

subroutine calc_therm_flame_sheet_boundary(xi,T, E,MWave,kappa,mu)!CH4/Air{{{
   use chem
   use func_therm
   double precision,intent(in)::xi
   double precision,intent(in)::T
   double precision,intent(out)::E
   double precision,intent(out)::MWave
   double precision,intent(out)::kappa
   double precision,intent(out)::mu

   double precision nCH4,nO2,nN2
   double precision nCO2 ,nH2O

   integer,parameter::iCH4 =7
   integer,parameter::iN2  =144
   integer,parameter::iO2  =157
   integer,parameter::iCO2 =13
   integer,parameter::iH2O =132

   integer,parameter::itCH4 =2
   integer,parameter::itN2  =19
   integer,parameter::itO2  =24
   integer,parameter::itCO2 =5
   integer,parameter::itH2O =14

   double precision,parameter::etaO2=0.215d0
   double precision,parameter::etaN2=1d0-etaO2

   integer,parameter:: ind(5) = (/ iCH4, iN2, iO2, iCO2, iH2O/)
   integer,parameter::indt(5) = (/itCH4,itN2,itO2,itCO2,itH2O/)
   integer            secN(5)
   double precision   molN(5)
   double precision    muN(5)

   double precision MWf,MWo
   double precision xi_st

   integer sect
   double precision cprnow,cvrnow,ERTnow,logT,sum_n,phi

   integer i,j

   !set MW
   MWf=1d-3* MW(iCH4)
   MWo=1d-3*(MW(iO2)*etaO2+MW(iN2)*etaN2)

   !set xi_st
   xi_st=MWf*etaO2/(MWf*etaO2+2d0*MWo)

   !set n
   nN2  = (1d0-xi)/MWo*etaN2
   if(xi > xi_st) then !fuel rich
      nO2  = 0d0
      nCH4 = (xi*(2d0*MWo+MWf*etaO2)-MWf*etaO2)/(2d0*MWf*MWo)
      nCO2 = 0.5d0*(1d0-xi)/MWo*etaO2
   else !oxydizer rich
      nO2  = (MWf*etaO2-xi*(2d0*MWo+MWf*etaO2))/(MWf*MWo)
      nCH4 = 0d0
      nCO2 = xi/MWf
   end if
   nH2O = 2d0*nCO2

   molN(1)=nCH4
   molN(2)=nN2
   molN(3)=nO2
   molN(4)=nCO2
   molN(5)=nH2O

   !calc T
   logT=log(T)

   ERTnow=0d0
   cvrnow=0d0
   do i=1,5
      call check_section_number(ind(i),T,secN(i))
      ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
      cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
   end do

   sum_n=nCH4+nO2+nN2+nCO2+nH2O
   MWave=1d3/sum_n

   !calc kappa
   cprnow= cvrnow+sum_n
   kappa=cprnow/cvrnow
   E=ERTnow*Ru*T

   !calc mu
   do i=1,5
      call check_section_number_trans(indt(i),T,sect)
      muN(i) = calc_mu(indt(i),T,logT,sect)
   end do

   mu = 0d0
   do i=1,5
      denom = 0d0
      do j=1,5
         phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MW(ind(j)) / MW(ind(i)) )**0.25d0 )**2 &
                      * sqrt(2d0 *MW(ind(j))/(MW(ind(i))+MW(ind(j))))
         denom = denom + molN(j) * phi
      end do
      mu = mu + molN(i) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine calc_therm_flame_sheet_boundary!}}}

subroutine calc_therm_flame_sheet_inert(xi,E, T, MWave,kappa,mu)!CH4/Air{{{
   use chem
   use func_therm
   double precision,intent(in)::xi
   double precision,intent(in)::E
   double precision,intent(inout)::T
   double precision,intent(out)::MWave
   double precision,intent(out)::kappa
   double precision,intent(out)::mu

   double precision nCH4,nO2,nN2

   integer,parameter::iCH4 =7
   integer,parameter::iN2  =144
   integer,parameter::iO2  =157

   integer,parameter::itCH4 =2
   integer,parameter::itN2  =19
   integer,parameter::itO2  =24

   double precision,parameter::etaO2=0.215d0
   double precision,parameter::etaN2=1d0-etaO2

   integer,parameter:: ind(3) = (/ iCH4, iN2, iO2/)
   integer,parameter::indt(3) = (/itCH4,itN2,itO2/)
   integer            secN(3)
   double precision   molN(3)
   double precision    muN(3)

   double precision MWf,MWo

   double precision T_old

   integer sect
   double precision cprnow,cvrnow,ERTnow,logT,sum_n,phi

   integer,parameter::max_step = 200
   integer i,j,k

   !set MW
   MWf=1d-3* MW(iCH4)
   MWo=1d-3*(MW(iO2)*etaO2+MW(iN2)*etaN2)

   !set n
   nCH4 =      xi /MWf
   nO2  = (1d0-xi)/MWo*etaO2
   nN2  = (1d0-xi)/MWo*etaN2

   molN(1)=nCH4
   molN(2)=nN2
   molN(3)=nO2

   !calc T
   logT=log(T)

   ERTnow=0d0
   cvrnow=0d0
   do i=1,3
      call check_section_number(ind(i),T,secN(i))
      ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
      cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
   end do

   k=0
   T_old = T*0.5d0
   do while(abs(E/Ru-ERTnow*T)>cvrnow*1d-10 .and. abs(T-T_old)>1d-8.and. k<max_step)
      T_old =T
      T     =T+omega*(E/Ru-ERTnow*T)/cvrnow
      logT=log(T)

      ERTnow=0d0
      cvrnow=0d0
      do i=1,3
         call check_section_number(ind(i),T,secN(i))
         ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
         cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
      end do
      k=k+1
   end do

   if(k >= max_step) then
      print *, "Not Converted at calc_therm_flame_sheet_inert",T
      call exit(1)
   end if


   sum_n=nCH4+nO2+nN2
   MWave=1d3/sum_n

   !calc kappa
   cprnow= cvrnow+sum_n
   kappa=cprnow/cvrnow

   !calc mu
   do i=1,3
      call check_section_number_trans(indt(i),T,sect)
      muN(i) = calc_mu(indt(i),T,logT,sect)
   end do

   mu = 0d0
   do i=1,3
      denom = 0d0
      do j=1,3
         phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MW(ind(j)) / MW(ind(i)) )**0.25d0 )**2 &
                      * sqrt(2d0 *MW(ind(j))/(MW(ind(i))+MW(ind(j))))
         denom = denom + molN(j) * phi
      end do
      mu = mu + molN(i) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine calc_therm_flame_sheet_inert!}}}

subroutine calc_therm_flame_sheet_inert_boundary(xi,T, E,MWave,kappa,mu)!CH4/Air{{{
   use chem
   use func_therm
   double precision,intent(in)::xi
   double precision,intent(in)::T
   double precision,intent(out)::E
   double precision,intent(out)::MWave
   double precision,intent(out)::kappa
   double precision,intent(out)::mu

   double precision nCH4,nO2,nN2

   integer,parameter::iCH4 =7
   integer,parameter::iN2  =144
   integer,parameter::iO2  =157

   integer,parameter::itCH4 =2
   integer,parameter::itN2  =19
   integer,parameter::itO2  =24

   double precision,parameter::etaO2=0.215d0
   double precision,parameter::etaN2=1d0-etaO2

   integer,parameter:: ind(3) = (/ iCH4, iN2, iO2/)
   integer,parameter::indt(3) = (/itCH4,itN2,itO2/)
   integer            secN(3)
   double precision   molN(3)
   double precision    muN(3)

   double precision MWf,MWo
   double precision xi_st

   integer sect
   double precision cprnow,cvrnow,ERTnow,logT,sum_n,phi

   integer i,j

   !set MW
   MWf=1d-3* MW(iCH4)
   MWo=1d-3*(MW(iO2)*etaO2+MW(iN2)*etaN2)

   !set n
   nCH4 =      xi /MWf
   nO2  = (1d0-xi)/MWo*etaO2
   nN2  = (1d0-xi)/MWo*etaN2

   molN(1)=nCH4
   molN(2)=nN2
   molN(3)=nO2

   !calc T
   logT=log(T)

   ERTnow=0d0
   cvrnow=0d0
   do i=1,3
      call check_section_number(ind(i),T,secN(i))
      ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
      cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
   end do

   sum_n=nCH4+nO2+nN2
   MWave=1d3/sum_n

   !calc kappa
   cprnow= cvrnow+sum_n
   kappa=cprnow/cvrnow
   E=ERTnow*Ru*T

   !calc mu
   do i=1,3
      call check_section_number_trans(indt(i),T,sect)
      muN(i) = calc_mu(indt(i),T,logT,sect)
   end do

   mu = 0d0
   do i=1,3
      denom = 0d0
      do j=1,3
         phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MW(ind(j)) / MW(ind(i)) )**0.25d0 )**2 &
                      * sqrt(2d0 *MW(ind(j))/(MW(ind(i))+MW(ind(j))))
         denom = denom + molN(j) * phi
      end do
      mu = mu + molN(i) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine calc_therm_flame_sheet_inert_boundary!}}}
!end for CH4/Air!!!!!!!!!!!

!for Air!!!!!!!!!!!!!!!!!!!
subroutine calc_therm_air(E, T, MWave,kappa,mu)!Air{{{
   use chem
   use func_therm
   double precision,intent(in)::E
   double precision,intent(inout)::T
   double precision,intent(out)::MWave
   double precision,intent(out)::kappa
   double precision,intent(out)::mu

   double precision nO2,nN2

   integer,parameter::iN2  =144
   integer,parameter::iO2  =157

   integer,parameter::itN2  =52
   integer,parameter::itO2  =59

   double precision,parameter::etaO2=0.215d0
   double precision,parameter::etaN2=1d0-etaO2

   integer,parameter:: ind(2) = (/ iN2, iO2/)
   integer,parameter::indt(2) = (/itN2,itO2/)
   integer            secN(2)
   double precision   molN(2)
   double precision    muN(2)

   double precision MWo

   double precision T_old

   integer sect
   double precision cprnow,cvrnow,ERTnow,logT,sum_n,phi

   integer i,j,k

   !set MW
   MWo=1d-3*(MW(iO2)*etaO2+MW(iN2)*etaN2)

   !set n
   nN2  = etaN2/MWo
   nO2  = etaO2/MWo

   molN(1)=nN2
   molN(2)=nO2

   !calc T
   logT=log(T)

   ERTnow=0d0
   cvrnow=0d0
   do i=1,2
      call check_section_number(ind(i),T,secN(i))
      ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
      cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
   end do

   k=0
   T_old = T*0.5d0
   do while(abs(E/Ru-ERTnow*T)>cvrnow*1d-10 .and. log(T-T_old)>1d-8.and. k<50)
      T_old =T
      T     =T+omega*(E/Ru-ERTnow*T)/cvrnow
      if(T<0d0) then
         print *,"Negative Temperature at calc_therm_air. T=",T,"T_old=",T_old
         call exit(1)
      end if
      logT=log(T)

      ERTnow=0d0
      cvrnow=0d0
      do i=1,2
         call check_section_number(ind(i),T,secN(i))
         ERTnow = ERTnow + ert(ind(i),T,logT,secN(i))*molN(i)
         cvrnow = cvrnow + cvr(ind(i),T,     secN(i))*molN(i)
      end do
      k=k+1
   end do

   if(k >= 50) then
      print *, "Not Converted at calc_therm_flame_sheet"
      call exit(1)
   end if


   sum_n=nO2+nN2
   MWave=1d3/sum_n

   !calc kappa
   cprnow= cvrnow+sum_n
   kappa=cprnow/cvrnow

   !calc mu
   do i=1,2
      call check_section_number_trans(indt(i),T,sect)
      muN(i) = calc_mu(indt(i),T,logT,sect)
   end do

   mu = 0d0
   do i=1,2
      denom = 0d0
      do j=1,2
         phi = 0.25d0 * (1d0+sqrt(muN(i)/muN(j)) * (MW(ind(j)) / MW(ind(i)) )**0.25d0 )**2 &
                      * sqrt(2d0 *MW(ind(j))/(MW(ind(i))+MW(ind(j))))
         denom = denom + molN(j) * phi
      end do
      mu = mu + molN(i) * muN(i) / denom
   end do
   mu = mu * 1d-7
end subroutine calc_therm_air!}}}
!end for Air!!!!!!!!!!!!!!!
