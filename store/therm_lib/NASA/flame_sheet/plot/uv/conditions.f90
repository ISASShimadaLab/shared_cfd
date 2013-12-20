module conditions
   use chem
   use chem_var
   implicit none
   double precision,allocatable,dimension(:)::list_T
   double precision,allocatable,dimension(:)::list_Yf
   integer numT
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

   chartype=trim(next_word())
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
subroutine calcTeq!{{{
   double precision p,H,T,Tini, MWave,kappa,mu,rho,E
   double precision Y(2),DHi(2)
   double precision,dimension(3)::n,vhi,Yv
   integer i,j

   !start calculations
   T=To
   if(chartype .eq. "T ") then
      do j=1,numYf
         Y(1)=list_Yf(j)
         Y(2)=1d0-Y(1)
         call print_head(j)
         write(55,*) "#Yf=",list_Yf(j)
         do i=1,numT
            Tini=list_T(i)
            call calcE(Tini,Y,E)

            call flame_sheet(Y,E, T, MWave,kappa,mu,DHi,Yv,vhi)
            write(55,'(100es15.7)') Tini,Y,T,E,MWave,kappa,mu
         end do
         close(55)
      end do
   else
      do j=1,numT
         Tini=list_T(j)
         call print_head(j)
         write(55,*) "#T=",list_T(j)
         do i=1,numYf
            Y(1)=list_Yf(i)
            Y(2)=1d0-Y(1)
            call calcE(Tini,Y,E)

            call flame_sheet(Y,E, T, MWave,kappa,mu,DHi,Yv,vhi)
            write(55,'(100es15.7)') Tini,Y,T,E,MWave,kappa,mu
         end do
         close(55)
      end do
   end if
contains
subroutine print_head(j)
   integer j
   character*50 line
   write(line,'("plot.",i3.3,".dat")') j
   open(55,file=trim(line))
   write(55,'(a)') "# Tini = Initial Temperature, &
             &MF = mass fraction, oxid = oxidizer, prod = product, &
             &br = before reaction, &
             &MW = Molecular Weight, kappa = specific heat ratio, &
             &vis = viscosity coefficient"
   write(55,'(a1,a14,100a15)') "#","Tini(K)","br fuel MF",&
                  "br oxid MF","Temperature(K)","Energy(J/kg)",&
                  "MW (g/mol)","kappa","vis(Pa*s)"
end subroutine print_head
end subroutine calcTeq!}}}
subroutine out_plt!{{{
   integer j,k
   double precision Yf,Tini

   open( 55,file="plot.plt")
   write(55,'(a)') "#!/usr/bin/gnuplot"
   write(55,'(a)') "set terminal postscript enhanced color"
   write(55,'(a)') "set ylabel 'Temperature(K)'"
   if(chartype .eq. "T ") then
      call sub("Initial Temperature(K)","1",&
               "Mole Fraction of fuel",numYf,list_Yf,'f6.3')
   else
      call sub("Mole Fraction of fuel","2",&
               "Initial Temperature(K)",numT, list_T, 'f5.0')
   end if
   close(55)
contains
subroutine sub(stringa,chara,stringc,numc,listc,formc)
   implicit none
   character(*)                 ,intent(in)::stringa
   character                    ,intent(in)::chara
   character(*)                 ,intent(in)::stringc
   integer                      ,intent(in)::numc
   double precision,dimension(:),intent(in)::listc
   character(*)                 ,intent(in)::formc

   integer j
   write(55,'(a)') "set xlabel '"//trim(stringa)//"'"
   write(55,'(a,i3.3,a)') 'set output "out.eps"'
   write(55,'(a)') 'set title "Equilibrium Temperature"'
   write(55,'(a)') 'plot \'
   do j=1,numc-1
      write(55,'(a,i3.3,a,'//trim(formc)//',a)') &
              '"plot.',j,'.dat" u '//chara//':4  w l title "'//trim(stringc)//'=',listc(j),'",\'
   end do
   write(55,'(a,i3.3,a,'//trim(formc)//',a)') &
           '"plot.',j,'.dat" u '//chara//':4  w l title "'//trim(stringc)//'=',listc(j),'"'
end subroutine sub
end subroutine out_plt!}}}
end module conditions
