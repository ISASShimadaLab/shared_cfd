!init_package
subroutine init_therm!{{{
   call read_chemkin_parameter
   call set_reac_and_therm_data
   call set_trans_data
   call check_chem_and_flow
   call read_fo_composition
end subroutine init_therm!}}}

!read control.inp
subroutine read_chemkin_parameter!{{{
   use chem
   use mod_mpi
   use chem_var
   implicit none

   if(myid .eq. 0) then
      open(8,file="control_chem.inp")
      read(8,'()')
      read(8,'(25x,i10)')    ns_tocalc
      close(8)

      !!! adjust data !!!
      select case(ns_tocalc)
         case(:-1)
            ns_tocalc = ns
         case(0,ns+1:)
            print *,"Odd Number to calculation at flame sheet. : value = ",ns_tocalc
            stop
      end select
   end if

   !print *,"ns_tocalc = ",ns_tocalc

   !!! MPI COMMUNICATIONS
   call MPI_Bcast(ns_tocalc,1,          MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end subroutine read_chemkin_parameter!}}}
subroutine read_fo_composition!{{{
   use chem
   use mod_mpi
   use chem_var
   implicit none
   character*100 buf,buf2
   double precision amount
   integer ind

   double precision DHit(nY)
   integer i

   if(myid .eq. 0) then
      vrhoo =1d-20
      vrhof =1d-20

      open(8,file="control_chem.inp")
      read(8,'()')
      read(8,'()')
      read(8,'(25x,es15.7)') po
      read(8,'(25x,es15.7)') To
      read(8,'()')
      do
         read(8,'(a)') buf
         if(buf(1:3) .eq. 'end') exit
         buf=adjustl(trim(buf)) !delete front and back spaces
         ind=index(buf,' ')
         if(ind .eq. 0) stop "Bad format at control_chem.inp while reading Oxygen Compositions"
         read(buf(ind+1:),*) amount
         vrhoo(search_species(buf(1:ind-1))) = amount
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
         vrhof(search_species(buf(1:ind-1))) = amount
      end do
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


   !!! calculate other properties for oxygen and fuel
   call rho_rel2abs(po,To, vrhoo, rhoo,Eo)
   call calc_T(vrhoo,rhoo,Eo, To, kappao,MWo,DHit,vhio,muo)
   vwo(1:nY) = vrhoo(1:nY)/rhoo

   call rho_rel2abs(pf,Tf, vrhof, rhof,Ef)
   call calc_T(vrhof,rhof,Ef, Tf, kappaf,MWf,DHit,vhif,muf)
   vwf(1:nY) = vrhof(1:nY)/rhof
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
end subroutine read_fo_composition!}}}

!data consistency check between flow and chemistry
subroutine check_chem_and_flow!{{{
   use grbl_prmtr
   use const_chem
   implicit none
   if(nY .ne. ns) then
      print *,"Error. nY and ns is different. (nY,ns)=",nY,ns
      stop
   end if
end subroutine check_chem_and_flow!}}}

! flame sheet model
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
subroutine set_thermo_prop!{{{
   use mod_mpi
   use prmtr
   use variable
   use chem_var
   implicit none
   double precision ei,T,MW,kappa,xi
   integer i,j

   !$omp parallel do private(i,ei,T,MW,kappa,xi)
   do j=nys,nye
      do i=nxs,nxe
         !values necessary to calculate thermo values
         ei=q(nY+3,i,j)/w(1,i,j)-0.5d0*(w(2,i,j)**2+w(3,i,j)**2)
         T=w(4,i,j)/(w(1,i,j)*w(indxR,i,j))

         !calc therm
         call calc_T(     q(1:nY,i,j),w(1,i,j),ei, T, kappa,MW,DHi(:,i,j),vhi(:,i,j),w(indxMu,i,j))

         !calculate R_gas and w
         w(indxR,   i,j)=R_uni/MW
         w(     4,i,j)=w(1,i,j)*w(indxR,i,j)*T
         w(indxg, i,j)=kappa
         w(indxht,i,j)=(q(nY+3,i,j)+w(4,i,j))/w(1,i,j)
      end do
   end do
   !$omp end parallel do
end subroutine set_thermo_prop!}}}

! reaction
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
      n(j) = max(n(j),1d-20)
   end do

   if(T < 0d0) then
      print *,"negative temperature=",T
      call exit(0)
   end if

   if(tout < 0d0) then
      print *,"negative dt=",tout
      call exit(0)
   end if

   !set parameters
   n(ns+1)= T
   tt          = 0d0
   istate      = 1
   neq         = ns+1
   rtol        = 1d-10
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
         atol=1d-12
         num_recalc = num_recalc+1
         if(num_recalc>max_recalc) then
            print *,"Exceed max_recalc=",max_recalc," istate = ",istate
            call exit(1)
         end if
      else
         print *,"istate number =",istate
         call exit(1)
      end if
   end do

   !!set T
   !T=n(ns+1)

   !n to vrho
   do j=1,ns
      n(j) = max(n(j),1d-20)
      vrho(j)=n(j)*MWs(j)*1d3
   end do
end subroutine reaction!}}}
subroutine proceed_reaction!{{{
   use mod_mpi
   use variable
   implicit none
   integer i,j
   double precision T

   !$omp parallel do private(i,T)
   do j=nys,nye
      do i=nxs,nxe
         T=w(4,i,j)/w(1,i,j)/w(indxR,i,j)
         if(T>1500d0 .or. (w(5,i,j)<1d0-1d-11 .and. T>700d0)) then
            call reaction(T,q(1:nY,i,j),dt_grbl)
         !   write(myid+100,'(2(i4.4,x),2(es15.7,x),i1,x,i3.3)') i,j,x(i,j),r(i,j),1,myid
         !else
         !   write(myid+100,'(2(i4.4,x),2(es15.7,x),i1,x,i3.3)') i,j,x(i,j),r(i,j),0,myid
         end if
      end do
   end do
   !$omp end parallel do

   !close(myid+100)

   !call MPI_Barrier(MPI_COMM_WORLD)
   !call exit(0)
end subroutine proceed_reaction!}}}
subroutine distribute_reaction!{{{
   use mod_mpi
   use variable
   implicit none
   integer i,j,k,kk
   double precision T

   double precision tmp_send((nY+1)*bwmax*bwmax)
   integer data_send(2,bwmax*bwmax)
   integer tocalc

   integer num_send( 2,ng)
   integer num_receive(0:ng-1)
   integer list(2,1:ng)
   integer buf,buf2
   integer nsend,ndata
   integer bwrec,resrec
   integer res1,res2

   integer lists(2,1:ng)
   integer listr(2,1:ng)
   integer nums,numr

   double precision qrecv(nY+1,bwmax*bwmax)
   integer v_ireq(bwmax*bwmax)
   integer v_ierr(bwmax*bwmax)

   !!!! search S calculation
   nums=0
   tocalc=0
   do j=nys,nye
      do i=nxs,nxe
         T=w(4,i,j)/w(1,i,j)/w(indxR,i,j)
         if(T>1500d0 .or. (w(5,i,j)<1d0-1d-11 .and. T>700d0)) then
            tmp_send(tocalc*(nY+1)+1: tocalc*(nY+1)+nY)=q(1:nY,i,j)
            tmp_send(tocalc*(nY+1)+nY+1)=T
            tocalc=tocalc + 1
            data_send(1,tocalc)=i
            data_send(2,tocalc)=j
         end if
      end do
   end do
   !!!! end search S calculation

   if(myid .eq. 0) then
      ! receive each data
      nsend=0
      ndata=0
      do i=0,ng-1
         if(i .eq. 0) then
            buf=tocalc
         else
            call MPI_Recv(buf,1,MPI_INTEGER,i,0,MPI_COMM_WORLD, istatus, ierr)
         end if

         if(buf>0) then
            nsend=nsend+1
            ndata=ndata+buf
            num_send(1,nsend)=i
            num_send(2,nsend)=buf
         end if
      end do
      ! end receive each data

      !calc num_receive
      bwrec=ndata/ng
      resrec=ndata-ng*bwrec
      do i=0,ng-1
         if(i<resrec) then
            num_receive(i)=bwrec+1
         else
            num_receive(i)=bwrec
         end if
      end do

      !make send list
      j=0
      res2=0
      do i=1,nsend
         res1=num_send(2,i)
         k=0
         do
            k=k+1
            if(num_receive(j)-res2>=res1) then
               list(1,k)=j
               list(2,k)=res1

               res2=res2+res1
               if(num_receive(j) .eq. res2) then
                  res2=0
                  j=j+1
               end if
               exit
            else
               list(1,k)=j
               list(2,k)=num_receive(j)-res2

               res1=res1-list(2,k)
               res2=0
               j=j+1
            end if
         end do

         if(num_send(1,i) .eq. 0) then
            nums        =k
            lists(:,1:k)=list(:,1:k)
         else
            call MPI_Send(   k,  1,MPI_INTEGER,num_send(1,i),0,MPI_COMM_WORLD,ierr)
            call MPI_Send(list,2*k,MPI_INTEGER,num_send(1,i),0,MPI_COMM_WORLD,ierr)
         end if
      end do
      !end make send list

      !make receive list
      j=1
      res2=0
      do i=0,ng-1
         res1=num_receive(i)
         k=0
         do
            k=k+1
            if(num_send(2,j)-res2>=res1) then
               list(1,k)=num_send(1,j)
               list(2,k)=res1

               res2=res2+res1
               if(num_send(2,j) .eq. res2) then
                  res2=0
                  j=j+1
               end if
               exit
            else
               list(1,k)=num_send(1,j)
               list(2,k)=num_send(2,j)-res2

               res1=res1-list(2,k)
               res2=0
               j=j+1
            end if
         end do

         if(i .eq. 0) then
            numr        =k
            listr(:,1:k)=list(:,1:k)
         else
            call MPI_Send(k,     1,MPI_INTEGER,i,0,MPI_COMM_WORLD,ierr)
            call MPI_Send(list,2*k,MPI_INTEGER,i,0,MPI_COMM_WORLD,ierr)
         end if
      end do
      !end make receive list
   else
      call MPI_Send(tocalc,  1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
      if(tocalc .ne. 0) then
         call MPI_Recv( nums,     1,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
         call MPI_Recv(lists,2*nums,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
      end if
      call MPI_Recv( numr,     1,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
      call MPI_Recv(listr,2*numr,MPI_INTEGER,0,0,MPI_COMM_WORLD, istatus, ierr)
   end if

   !!!!!!!!!!!!!!!!!!!! SEND, CALC, RECV !!!!!!!!!!!!!!!!!
   buf=1
   do i=1,nums
      if(lists(1,i) .ne. myid) then
         call MPI_Isend(tmp_send(buf),lists(2,i)*(nY+1), MPI_DOUBLE_PRECISION,&
                         lists(1,i),0,MPI_COMM_WORLD,v_ireq(i),v_ierr(i))
         buf=buf+lists(2,i)*(nY+1)
      else
         buf2=0
         do j=1,numr
            if(listr(1,j) .ne. myid) then
               buf2=buf2+listr(2,j)
            else
               do k=1,listr(2,j)
                  qrecv(:,buf2+k)=tmp_send(buf:buf+nY)
                  buf=buf+nY+1
               end do
               exit
            end if
         end do
      end if
   end do

   buf=1
   do i=1,numr
      if(listr(1,i) .ne. myid) then
         call MPI_Recv(qrecv(1,buf),listr(2,i)*(nY+1), MPI_DOUBLE_PRECISION,&
                         listr(1,i),0,MPI_COMM_WORLD, istatus, ierr)
      end if
      buf=buf+listr(2,i)
   end do

   !$omp parallel do private(i)
   do i=1,buf-1
      call reaction(qrecv(nY+1,i),qrecv(1:nY,i),dt_grbl)
   end do
   !$omp end parallel do

   do i=1,nums
      if(lists(1,i) .ne. myid) call MPI_Wait(v_ireq(i),istatus,ierr)
   end do

   buf=1
   do i=1,numr
      if(listr(1,i) .ne. myid) then
         call MPI_Isend(qrecv(1,buf),listr(2,i)*(nY+1),MPI_DOUBLE_PRECISION,&
                         listr(1,i),0,MPI_COMM_WORLD,v_ireq(i),v_ierr(i))
      else
         buf2=0
         do j=1,nums
            if(lists(1,j) .ne. myid) then
               buf2=buf2+lists(2,j)
            else
               do k=1,lists(2,j)
                  q(1:nY,data_send(1,buf2+k),&
                         data_send(2,buf2+k))&
                         =qrecv(1:nY,buf-1+k)
               end do
               exit
            end if
         end do
      end if
      buf=buf+listr(2,i)
   end do

   buf=0
   do i=1,nums
      if(lists(1,i) .ne. myid) then
         call MPI_Recv(tmp_send,lists(2,i)*(nY+1), MPI_DOUBLE_PRECISION,&
                         lists(1,i),0,MPI_COMM_WORLD, istatus,ierr)
         do j=1,lists(2,i)
            q(1:nY,data_send(1,buf+j),data_send(2,buf+j))&
                    =tmp_send((j-1)*(nY+1)+1:(j-1)*(nY+1)+nY)
         end do
      end if
      buf=buf+lists(2,i)
   end do

   do i=1,numr
      if(listr(1,i) .ne. myid) call MPI_Wait(v_ireq(i),istatus,ierr)
   end do

   !print '(a,i3.2,a)',"ID:",myid," OK"
   !call MPI_Barrier(MPI_COMM_WORLD,ierr)
   !call exit(0)
end subroutine distribute_reaction!}}}

!initial condition
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

!boundary condition
subroutine calc_boundary(p,T,deg, wt,vhi)!{{{
   use grbl_prmtr
   use prmtr
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
   wt(indxR) =R_uni/MWave
end subroutine calc_boundary!}}}
