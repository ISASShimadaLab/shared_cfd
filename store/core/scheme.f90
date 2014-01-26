subroutine set_w!{{{
   use mod_mpi
   use grbl_prmtr
   use variable
   implicit none
   integer i,j,k,plane
   double precision tmp

   !w1=rho, w2=u, w3=v, w4=p, w_(dimq+1)=gamma, w_(dimq+2)=ht(J/kg)
   do plane = nps,npe
      !$omp parallel do private(i,k,tmp)
      do j=nys(plane),nye(plane)
         do i=nxs(plane),nxe(plane)
            tmp =0d0
            do k=1,nY
               tmp=tmp+q(k,i,j,plane)
            end do
            w(1,i,j,plane)= tmp
            tmp=1d0/tmp
            w(2,i,j,plane)= q(nY+1,i,j,plane)*tmp
            w(3,i,j,plane)= q(nY+2,i,j,plane)*tmp
            do k=1,nY
               w(k+4,i,j,plane)=q(k,i,j,plane)*tmp
            end do
         end do
      end do
      !$omp end parallel do
   end do
end subroutine set_w!}}}

subroutine set_TG!{{{
   !************************************!
   ! set fH from uH which have already  !
   ! determined at set_HV               !
   !************************************!
   use mod_mpi
   use grbl_prmtr
   use variable
   implicit none
   integer i,j,plane

   !local variables
   double precision,dimension(dimw)::wHl
   double precision,dimension(dimw)::wHr
   double precision,dimension(dimq)::tg
   double precision,dimension(2)::vn

   double precision ul,vl,rhol,gmmal,unl,phil(dimq),al,Ml,betal,pl
   double precision ur,vr,rhor,gmmar,unr,phir(dimq),ar,Mr,betar,pr

   double precision N(dimq)
   double precision Drho,Dp,a_bar,un_bar_abs,p_bar
   double precision g,gp,gm,chi,M_tilde,dm,un_bar_abs_p,un_bar_abs_m
   double precision p_tilde

   integer k

   do plane = nps,npe
      ! i-direction{{{
      !$omp parallel do private(i,&
      !$omp                     wHl,wHr,vn,tg,&
      !$omp                     ul,vl,rhol,gmmal,unl,phil,al,Ml,betal,pl,&
      !$omp                     ur,vr,rhor,gmmar,unr,phir,ar,Mr,betar,pr,&
      !$omp                     N,&
      !$omp                     Drho,Dp,a_bar,un_bar_abs,p_bar,&
      !$omp                     g,gp,gm,chi,M_tilde,dm,un_bar_abs_p,un_bar_abs_m,&
      !$omp                     p_tilde,&
      !$omp                     k)
      do j=nys(plane),nye(plane)
         do i=nxs(plane)-1,nxe(plane)
            wHl=wHli(1:dimw,i,j,plane)
            wHr=wHri(1:dimw,i,j,plane)
            vn =vni(1:2,i,j,plane)
            !calc wlr
            rhol =whl(1)
            ul   =whl(2)
            vl   =whl(3)
            pl   =whl(4)
            gmmal=whl(indxg)
            do k=1,nY
               phil(k)=whl(4+k)
            end do
            phil(nY+1)=ul
            phil(nY+2)=vl
            phil(nY+3)=whl(indxht) !enthalpy
            al=sqrt(gmmal*pl/rhol)
            unl=vn(1)*ul+vn(2)*vl

            rhor =whr(1)
            ur   =whr(2)
            vr   =whr(3)
            pr   =whr(4)
            gmmar=whr(indxg)
            do k=1,nY
               phir(k)=whr(4+k)
            end do
            phir(nY+1)=ur
            phir(nY+2)=vr
            phir(nY+3)=whr(indxht) !enthalpy
            ar=sqrt(gmmar*pr/rhor)
            unr=vn(1)*ur+vn(2)*vr

                     !print '("phir")'
                     !print '(es15.7)',phir
                     !print '("phil")'
                     !print '(es15.7)',phil

            !set average and difference
            dp  =pr-pl
            a_bar =0.5d0*(ar +al )
            p_bar =0.5d0*(pr +pl )

            ml=unl/a_bar
            mr=unr/a_bar

            !set dm
            m_tilde=min(1d0,1d0/a_bar*sqrt(0.5d0*(ul**2+vl**2+ur**2+vr**2)))
            chi=(1d0-m_tilde)**2

            gp=-max(min(ml,0d0),-1d0)
            gm= min(max(mr,0d0), 1d0)
            g=gp*gm

            un_bar_abs  =(abs(unr)*rhor+abs(unl)*rhol)/(rhor+rhol)
            un_bar_abs_p=(1d0-g)*un_bar_abs+g*abs(unl)
            un_bar_abs_m=(1d0-g)*un_bar_abs+g*abs(unr)

            dm=0.5d0*(rhol*(unl+un_bar_abs_p)+rhor*(unr-un_bar_abs_m)-chi/a_bar*dp)

            !set p_tilde
            n(:)=0d0
            n(nY+1)=vn(1)
            n(nY+2)=vn(2)

            if(abs(mr)<1d0) then
               betar=0.25d0*(2d0+mr)*(mr-1d0)**2
            else
               betar=0.5d0*(1d0-sign(1d0,mr))
            end if

            if(abs(ml)<1d0) then
               betal=0.25d0*(2d0-ml)*(ml+1d0)**2
            else
               betal=0.5d0*(1d0+sign(1d0,ml))
            end if

            p_tilde=p_bar+0.5d0*(betal-betar)*(pl-pr)+(1d0-chi)*(betar+betal-1d0)*p_bar

            !calculate tg
            tg=0.5d0*(dm+abs(dm))*phil&
              +0.5d0*(dm-abs(dm))*phir&
              +p_tilde*n
            TGi(1:dimq,i,j,plane)=tg
         end do
      end do
      !$omp end parallel do
      !}}}

      ! j-direction{{{
      !$omp parallel do private(i,&
      !$omp                     ul,vl,rhol,gmmal,unl,phil,al,Ml,betal,pl,&
      !$omp                     wHl,wHr,vn,tg,&
      !$omp                     ur,vr,rhor,gmmar,unr,phir,ar,Mr,betar,pr,&
      !$omp                     N,&
      !$omp                     Drho,Dp,a_bar,un_bar_abs,p_bar,&
      !$omp                     g,gp,gm,chi,M_tilde,dm,un_bar_abs_p,un_bar_abs_m,&
      !$omp                     p_tilde,&
      !$omp                     k)
      do j=nys(plane)-1,nye(plane)
         do i=nxs(plane),nxe(plane)
            wHl=wHlj(1:dimw,i,j,plane)
            wHr=wHrj(1:dimw,i,j,plane)
            vn =vnj(1:2,i,j,plane)
            !calc wlr
            rhol =whl(1)
            ul   =whl(2)
            vl   =whl(3)
            pl   =whl(4)
            gmmal=whl(indxg)
            do k=1,nY
               phil(k)=whl(4+k)
            end do
            phil(nY+1)=ul
            phil(nY+2)=vl
            phil(nY+3)=whl(indxht) !enthalpy
            al=sqrt(gmmal*pl/rhol)
            unl=vn(1)*ul+vn(2)*vl

            rhor =whr(1)
            ur   =whr(2)
            vr   =whr(3)
            pr   =whr(4)
            gmmar=whr(indxg)
            do k=1,nY
               phir(k)=whr(4+k)
            end do
            phir(nY+1)=ur
            phir(nY+2)=vr
            phir(nY+3)=whr(indxht) !enthalpy
            ar=sqrt(gmmar*pr/rhor)
            unr=vn(1)*ur+vn(2)*vr

                     !print '("phir")'
                     !print '(es15.7)',phir
                     !print '("phil")'
                     !print '(es15.7)',phil

            !set average and difference
            dp  =pr-pl
            a_bar =0.5d0*(ar +al )
            p_bar =0.5d0*(pr +pl )

            ml=unl/a_bar
            mr=unr/a_bar

            !set dm
            m_tilde=min(1d0,1d0/a_bar*sqrt(0.5d0*(ul**2+vl**2+ur**2+vr**2)))
            chi=(1d0-m_tilde)**2

            gp=-max(min(ml,0d0),-1d0)
            gm= min(max(mr,0d0), 1d0)
            g=gp*gm

            un_bar_abs  =(abs(unr)*rhor+abs(unl)*rhol)/(rhor+rhol)
            un_bar_abs_p=(1d0-g)*un_bar_abs+g*abs(unl)
            un_bar_abs_m=(1d0-g)*un_bar_abs+g*abs(unr)

            dm=0.5d0*(rhol*(unl+un_bar_abs_p)+rhor*(unr-un_bar_abs_m)-chi/a_bar*dp)

            !set p_tilde
            n(:)=0d0
            n(nY+1)=vn(1)
            n(nY+2)=vn(2)

            if(abs(mr)<1d0) then
               betar=0.25d0*(2d0+mr)*(mr-1d0)**2
            else
               betar=0.5d0*(1d0-sign(1d0,mr))
            end if

            if(abs(ml)<1d0) then
               betal=0.25d0*(2d0-ml)*(ml+1d0)**2
            else
               betal=0.5d0*(1d0+sign(1d0,ml))
            end if

            p_tilde=p_bar+0.5d0*(betal-betar)*(pl-pr)+(1d0-chi)*(betar+betal-1d0)*p_bar

            !calculate tg
            tg=0.5d0*(dm+abs(dm))*phil&
              +0.5d0*(dm-abs(dm))*phir&
              +p_tilde*n
            TGj(1:dimq,i,j,plane)=tg
         end do
      end do
      !$omp end parallel do
      !}}}
   end do
end subroutine set_TG!}}}


