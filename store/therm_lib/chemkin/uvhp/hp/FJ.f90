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
   double precision sn,logsn
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
   double precision,dimension(ns)::vmu0rt,vhrt,vcpr
   double precision vlogM,Dr,r
   double precision vlogk(2)
   double precision dh,cp,sn,logsn
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
      temp2=0d0
      temp3=0d0
      do i=1,7
         tmp=coeff(i,sec(j),j)
         temp1=temp1+tmp*cmurt(i)
         temp2=temp2+tmp*chrt( i)
         temp3=temp3+tmp*ccpr( i)
      end do
      vmu0rt(j)=temp1
      vhrt(  j)=temp2
      vcpr(  j)=temp3

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

   dh=0d0; cp=0d0
   do j=1,ns
      dh=dh+vhrt(j)*dndt(j)
      cp=cp+vcpr(j)*n(   j)
   end do
   dndt(ns+1)=-dh*T/cp
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
   double precision,dimension(ns)::vmu0rt,vhrt,vcpr
   double precision vlogM,Dr,r
   double precision vlogk(2)
   double precision dh,cp,sn,logsn
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
      temp2=0d0
      temp3=0d0
      do i=1,7
         tmp=coeff(i,sec(j),j)
         temp1=temp1+tmp*cmurt(i)
         temp2=temp2+tmp*chrt( i)
         temp3=temp3+tmp*ccpr( i)
      end do
      vmu0rt(j)=temp1
      vhrt(  j)=temp2
      vcpr(  j)=temp3

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

   dh=0d0; cp=0d0
   do j=1,ns
      dh=dh+vhrt(j)*(dndtp(j)+dndtm(j))
      cp=cp+vcpr(j)*n(   j)
   end do
   dndtp(ns+1)=max(-dh*T/cp,0d0)
   dndtm(ns+1)=min(-dh*T/cp,0d0)
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
   double precision,dimension(ns)::vmu0rt,vdmu0rt,vhrt,vcpr,vdcprdT
   double precision numu0rt,nudmu0rt,vlogM,vM,dlnkdlnM
   double precision vlogk(2),dlnkdT(2),r
   double precision Dr
   double precision dh,cp,dcp,logsn,sn
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
      temp3= 0d0
      temp4= 0d0
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
      vhrt(   j)=temp3
      vcpr(   j)=temp4
      vdcprdT(j)=temp5

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
   dh=0d0; cp=0d0; dcp=0d0
   temp1=0d0;temp2=0d0
   do j=1,ns
      temp1=temp1+   vhrt(j)*dndn(j,ns+1)
      dh   =dh   +   vhrt(j)*dndt(j)
      dcp  =dcp  +   vcpr(j)*dndt(j)
      cp   =cp   +   vcpr(j)*   n(j)
      temp2=temp2+vdcprdT(j)*   n(j)
   end do
   dndt(ns+1)=-dh*T/cp
   dndn(ns+1,ns+1)=-(T*temp1+dcp+dndt(ns+1)*temp2)/cp

   !dTdn
   do j=1,ns
      temp1=0d0
      do k=1,ns
         temp1=temp1+vhrt(k)*dndn(k,j)
      end do
      dndn(ns+1,j)=-(T*temp1+dndt(ns+1)*vcpr(j))/cp
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

