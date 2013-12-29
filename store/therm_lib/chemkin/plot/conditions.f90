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
   character*2 chartype
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
   read(line,*) chartype

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
end subroutine set_Tref!}}}
subroutine out_Tref!{{{
   integer i,j,k

   !start calculations
   if(chartype .eq. "P ") then
      do k=1,numYf
         do j=1,numT
            call print_head(j,k)
            write(55,*) "#T=",list_T(j)," Yf=",list_Yf(k)
            do i=1,numP
               write(55,'(100es15.7)') list_T(j),list_P(i),list_Yf(k),Tref(:,j,i,k)
            end do
            close(55)
         end do
      end do
   else if(chartype .eq. "T ") then
      do k=1,numYf
         do j=1,numP
            call print_head(j,k)
            write(55,*) "#p=",list_p(j)," Yf=",list_Yf(k)
            do i=1,numT
               write(55,'(100es15.7)') list_T(i),list_P(j),list_Yf(k),Tref(:,i,j,k)
            end do
            close(55)
         end do
      end do
   else
      do k=1,numP
         do j=1,numT
            call print_head(j,k)
            write(55,*) "#T=",list_T(j)," P=",list_p(k)
            do i=1,numYf
               write(55,'(100es15.7)') list_T(j),list_P(k),list_Yf(i),Tref(:,j,k,i)
            end do
            close(55)
         end do
      end do
   end if
contains
subroutine print_head(j,k)
   integer j,k
   character*50 line
   write(line,'("plot.",2(i3.3,"."),"dat")') j,k
   open(55,file=trim(line))
   write(55,'(a)') "# Tini = Initial Temperature, MF = mass fraction, &
             &oxid = oxidizer, prod = product, &
             &br = before reaction, Tign = Ignition time, &
             &Teq = Equilibrium Temperature"
   write(55,'(a1,a14,100a15)') "#","Tini(K)","pressure(Pa)","br fuel MF",&
                  "Tign(s)","Teq(K)"
end subroutine print_head
end subroutine out_Tref!}}}
subroutine out_plt!{{{
   integer j,k
   double precision p,Yf,Tini

   open( 55,file="plot.plt")
   write(55,'(a)') "#!/usr/bin/gnuplot"
   write(55,'(a)') "set terminal postscript enhanced color"
   if(chartype .eq. "P ") then
      write(55,'(a)') "set logscale y"
      write(55,'(a)') "set ylabel 'Time(sec)'"
      call sub("Pressure(Pa)","2",&
               "Mole Fraction of fuel", numYf,list_Yf,'f6.3',&
               "Initial Temperature(K)",numT, list_T, 'f5.0',&
               "Ignition Time(s) at ","Tign","4")
      write(55,'(a)') "reset"
      write(55,'(a)') "set ylabel 'Temperature(K)'"
      call sub("Pressure(Pa)","2",&
               "Mole Fraction of fuel", numYf,list_Yf,'f6.3',&
               "Initial Temperature(K)",numT, list_T, 'f5.0',&
               "Equilibrium Temperature(K) at ","Teq","5")
   else if(chartype .eq. "T ") then
      write(55,'(a)') "set logscale y"
      write(55,'(a)') "set ylabel 'Time(sec)'"
      call sub("Initial Temperature(K)","1",&
               "Mole Fraction of fuel",numYf,list_Yf,'f6.3',&
               "Pressure(Pa)",numP,list_p,'es9.1',&
               "Ignition Time(s) at ","Tign","4")
      write(55,'(a)') "reset"
      write(55,'(a)') "set ylabel 'Temperature(K)'"
      call sub("Initial Temperature(K)","1",&
               "Mole Fraction of fuel",numYf,list_Yf,'f6.3',&
               "Pressure(Pa)",numP,list_p,'es9.1',&
               "Equilibrium Temperature(K) at ","Teq","5")
   else
      write(55,'(a)') "set logscale y"
      write(55,'(a)') "set ylabel 'Time(sec)'"
      call sub("Mole Fraction of fuel","3",&
               "Pressure(Pa)",numP,list_p,'es9.1',&
               "Initial Temperature(K)",numT, list_T, 'f5.0',&
               "Ignition Time(s) at ","Tign","4")
      write(55,'(a)') "reset"
      write(55,'(a)') "set ylabel 'Temperature(K)'"
      call sub("Mole Fraction of fuel","3",&
               "Pressure(Pa)",numP,list_p,'es9.1',&
               "Initial Temperature(K)",numT, list_T, 'f5.0',&
               "Equilibrium Temperature(K) at ","Teq","5")
   end if
   close(55)
contains
subroutine sub(stringa,chara,stringb,numb,listb,formb,stringc,numc,listc,formc,title,stringd,chard)
   implicit none
   character(*)                 ,intent(in)::stringa
   character                    ,intent(in)::chara
   character(*)                 ,intent(in)::stringb
   integer                      ,intent(in)::numb
   double precision,dimension(:),intent(in)::listb
   character(*)                 ,intent(in)::formb
   character(*)                 ,intent(in)::stringc
   integer                      ,intent(in)::numc
   double precision,dimension(:),intent(in)::listc
   character(*)                 ,intent(in)::formc
   character(*)                 ,intent(in)::title
   character(*)                 ,intent(in)::stringd
   character                    ,intent(in)::chard

   integer k,j
   write(55,'(a)') "set xlabel '"//trim(stringa)//"'"
   if(numc>numb) then
      do k=1,numb
         write(55,'(a,i3.3,a)') 'set output "'//trim(stringd)//'.',k,'.eps"'
         write(55,'(a,'//trim(formb)//',a)') 'set title  "'//trim(title)//trim(stringb)//'= ',listb(k),'"'
         write(55,'(a)') 'plot \'
         do j=1,numc-1
            write(55,'(a,i3.3,a,i3.3,a,'//trim(formc)//',a)') &
                    '"plot.',j,'.',k,'.dat" u '//chara//':'//trim(chard)//'  w l title "'//trim(stringc)//'=',listc(j),'",\'
         end do
         write(55,'(a,i3.3,a,i3.3,a,'//trim(formc)//',a)') &
                 '"plot.',j,'.',k,'.dat" u '//chara//':'//trim(chard)//'  w l title "'//trim(stringc)//'=',listc(j),'"'
      end do
   else
      do j=1,numc
         write(55,'(a,i3.3,a)') 'set output "'//trim(stringd)//'.',j,'.eps"'
         write(55,'(a,'//trim(formc)//',a)') 'set title  "'//trim(title)//trim(stringc)//' = ',listc(j),'"'
         write(55,'(a)') 'plot \'
         do k=1,numb-1
            write(55,'(a,i3.3,a,i3.3,a,'//trim(formb)//',a)') &
                    '"plot.',j,'.',k,'.dat" u '//chara//':'//trim(chard)//'  w l title "'//trim(stringb)//'=',listb(k),'",\'
         end do
         write(55,'(a,i3.3,a,i3.3,a,'//trim(formb)//',a)') &
                 '"plot.',j,'.',k,'.dat" u '//chara//':'//trim(chard)//'  w l title "'//trim(stringb)//'=',listb(k),'"'
      end do
   end if
end subroutine sub
end subroutine out_plt!}}}
end module conditions

