subroutine set_geojac!{{{
   use grbl_prmtr
   use mod_mpi
   use variable
   implicit none
   double precision dxdxi,dydxi,dxdeta,dydeta
   double precision det
   integer i,j

   ! set geojaci{{{
   do j=nys,nye
      do i=nxs-1,nxe
         dxdxi  = x(i+1,j)-x(i,j)
         dydxi  = r(i+1,j)-r(i,j)

         dxdeta = (x(i+1,j+1)+x(i,j+1)&
                  -x(i+1,j-1)-x(i,j-1))*0.25d0
         dydeta = (r(i+1,j+1)+r(i,j+1)&
                  -r(i+1,j-1)-r(i,j-1))*0.25d0

         det= dxdxi*dydeta - dydxi*dxdeta

         geojaci(1,1,i,j) = dydeta/det!d( xi)/d(x)
         geojaci(1,2,i,j) =-dxdeta/det!d( xi)/d(y)
         geojaci(2,1,i,j) =-dydxi /det!d(eta)/d(x)
         geojaci(2,2,i,j) = dxdxi /det!d(eta)/d(y)
      end do
   end do
   !}}}

   ! set geojacj{{{
   do j=nys-1,nye
      do i=nxs,nxe
         dxdxi = (x(i+1,j+1)+x(i+1,j)&
                 -x(i-1,j+1)-x(i-1,j))*0.25d0
         dydxi = (r(i+1,j+1)+r(i+1,j)&
                 -r(i-1,j+1)-r(i-1,j))*0.25d0

         dxdeta= x(i,j+1)-x(i,j)
         dydeta= r(i,j+1)-r(i,j)


         det= dxdxi*dydeta - dydxi*dxdeta

         geojacj(1,1,i,j) = dydeta/det!d( xi)/d(x)
         geojacj(1,2,i,j) =-dxdeta/det!d( xi)/d(y)
         geojacj(2,1,i,j) =-dydxi /det!d(eta)/d(x)
         geojacj(2,2,i,j) = dxdxi /det!d(eta)/d(y)
      end do
   end do
   !}}}
end subroutine set_geojac!}}}

subroutine set_TGv!{{{
   use grbl_prmtr
   use mod_mpi
   use variable
   implicit none
   double precision mu,k_heat
   double precision,dimension(dimq)::Ev,Fv
   double precision,dimension(dimw)::dwdxi,dwdeta,dwdx,dwdr
   double precision,dimension(0:ni+1,0:nj+1)::Tdeg
   double precision dudxi,dudeta
   double precision dudx, dudr
   double precision dvdxi,dvdeta
   double precision dvdx, dvdr
   double precision dTdxi,dTdeta
   double precision dTdx, dTdr
   double precision rho,u,v,T
   double precision Rgas,gamm
   double precision t_xx,t_xr,t_rr
   double precision qx,qr
   double precision divu
   double precision lhi(nY),hidwdx,hidwdr

   double precision D
   double precision,parameter::Sc=1d0
   double precision,parameter::Pr=1d0
   integer i,j,k

   !set Tdeg{{{
   !$omp parallel do default(none) shared(j,w,Tdeg,nxs,nxe,nys,nye) private(i)
   do j=nys-1,nye+1
      do i=nxs-1,nxe+1
         Tdeg(i,j)= w(4,i,j)/w(1,i,j)/w(indxR,i,j)
      end do
   end do
   !$omp end parallel do
   !}}}

   !set TGvi!{{{
   !$omp parallel do default(shared),&
   !$omp             private(i,k,&
   !$omp                     rho,u,v,gamm,T,Rgas,mu,&
   !$omp                     D,k_heat,&
   !$omp                     dudxi,dudeta,dudx,dudr,&
   !$omp                     dvdxi,dvdeta,dvdx,dvdr,&
   !$omp                     dTdxi,dTdeta,dTdx,dTdr,&
   !$omp                     dwdxi,dwdeta,dwdx,dwdr,&
   !$omp                     lhi,hidwdx,hidwdr,&
   !$omp                     divu,t_xx,t_xr,t_rr,qx,qr,&
   !$omp                     Ev,Fv)
   do j=nys,nye
      do i=nxs-1,nxe
         rho = (w(1    ,i,j)+w(1    ,i+1,j))*0.5d0
         u   = (w(2    ,i,j)+w(2    ,i+1,j))*0.5d0
         v   = (w(3    ,i,j)+w(3    ,i+1,j))*0.5d0
         gamm= (w(indxg,i,j)+w(indxg,i+1,j))*0.5d0
         T   =    (Tdeg(i,j)+   Tdeg(i+1,j))*0.5d0 !total energy
         Rgas=   (w(indxR,i,j)+  w(indxR,i+1,j))*0.5d0
         mu  =  (w(indxMu,i,j)+ w(indxMu,i+1,j))*0.5d0
         do k=1,nY
            lhi(k)=(vhi(k,i,j) +vhi(  k,i+1,j))*0.5d0
         end do

         D     =mu/rho/Sc
         k_heat=mu*gamm/(gamm-1d0)*Rgas/Pr               !k from mu by Prandtl number

         dudxi  =  w(2,i+1,j) -w(2,i,j)
         dvdxi  =  w(3,i+1,j) -w(3,i,j)
         dTdxi  = Tdeg(i+1,j)-Tdeg(i,j)
         do k=1,nY
            dwdxi(4+k)  =  w(4+k,i+1,j) -w(4+k,i,j)
         end do

         dudeta =  (w(2,i+1,j+1) +w(2,i,j+1)&
                   -w(2,i+1,j-1) -w(2,i,j-1))*0.25d0
         dvdeta =  (w(3,i+1,j+1) +w(3,i,j+1)&
                   -w(3,i+1,j-1) -w(3,i,j-1))*0.25d0
         dTdeta = (Tdeg(i+1,j+1)+Tdeg(i,j+1)&
                  -Tdeg(i+1,j-1)-Tdeg(i,j-1))*0.25d0
         do k=1,nY
            dwdeta(4+k) =  (w(4+k,i+1,j+1) +w(4+k,i,j+1)&
                           -w(4+k,i+1,j-1) -w(4+k,i,j-1))*0.25d0
         end do

         dudx=dudxi*geojaci(1,1,i,j)+dudeta*geojaci(2,1,i,j)
         dvdx=dvdxi*geojaci(1,1,i,j)+dvdeta*geojaci(2,1,i,j)
         dTdx=dTdxi*geojaci(1,1,i,j)+dTdeta*geojaci(2,1,i,j)
         dudr=dudxi*geojaci(1,2,i,j)+dudeta*geojaci(2,2,i,j)
         dvdr=dvdxi*geojaci(1,2,i,j)+dvdeta*geojaci(2,2,i,j)
         dTdr=dTdxi*geojaci(1,2,i,j)+dTdeta*geojaci(2,2,i,j)
         hidwdx=0d0;hidwdr=0d0
         do k=1,nY
            dwdx(4+k)=dwdxi(4+k)*geojaci(1,1,i,j)+dwdeta(4+k)*geojaci(2,1,i,j)
            hidwdx=hidwdx+dwdx(4+k)*lhi(k)
            dwdr(4+k)=dwdxi(4+k)*geojaci(1,2,i,j)+dwdeta(4+k)*geojaci(2,2,i,j)
            hidwdr=hidwdr+dwdr(4+k)*lhi(k)
         end do

         !!!2 dimentional-plane
         divu=dudx+dvdr

         t_xx=-2d0/3d0*mu*divu+2d0*mu*dudx
         t_rr=-2d0/3d0*mu*divu+2d0*mu*dvdr
         t_xr= mu*(dvdx+dudr)
         qx  = -k_heat*dTdx-rho*D*hidwdx
         qr  = -k_heat*dTdr-rho*D*hidwdr

         do k=1,nY
            Ev(k)=rho*D*dwdx(4+k)
         end do
         Ev(nY+1)=t_xx
         Ev(nY+2)=t_xr
         Ev(nY+3)=u*t_xx+v*t_xr-qx

         do k=1,nY
            Fv(k)=rho*D*dwdr(4+k)
         end do
         Fv(nY+1)=t_xr
         Fv(nY+2)=t_rr
         Fv(nY+3)=u*t_xr+v*t_rr-qr

         TGvi(:,i,j)=Ev*vni(1,i,j)+Fv*vni(2,i,j)
      end do
   end do
   !$omp end parallel do
   !}}}

   !set TGvj!{{{
   !$omp parallel do default(shared),&
   !$omp             private(i,k,&
   !$omp                     rho,u,v,gamm,T,Rgas,mu,&
   !$omp                     D,k_heat,&
   !$omp                     dudxi,dudeta,dudx,dudr,&
   !$omp                     dvdxi,dvdeta,dvdx,dvdr,&
   !$omp                     dTdxi,dTdeta,dTdx,dTdr,&
   !$omp                     dwdxi,dwdeta,dwdx,dwdr,&
   !$omp                     lhi,hidwdx,hidwdr,&
   !$omp                     divu,t_xx,t_xr,t_rr,qx,qr,&
   !$omp                     Ev,Fv)
   do j=nys-1,nye
      do i=nxs,nxe
         rho = (w(1    ,i,j)+w(1    ,i,j+1))*0.5d0
         u   = (w(2    ,i,j)+w(2    ,i,j+1))*0.5d0
         v   = (w(3    ,i,j)+w(3    ,i,j+1))*0.5d0
         gamm= (w(indxg,i,j)+w(indxg,i,j+1))*0.5d0
         T   =    (Tdeg(i,j)+   Tdeg(i,j+1))*0.5d0
         Rgas=   (w(indxR,i,j)+  w(indxR,i,j+1))*0.5d0
         mu  =  (w(indxMu,i,j)+ w(indxMu,i,j+1))*0.5d0
         do k=1,nY
            lhi(k)=(vhi(k,i,j) +vhi(  k,i,j+1))*0.5d0
         end do

         D     =mu/rho/Sc
         k_heat=mu*gamm/(gamm-1d0)*Rgas/Pr               !k from mu by Prandtl number

         dudxi  =  (w(2,i+1,j+1) +w(2,i+1,j)&
                   -w(2,i-1,j+1) -w(2,i-1,j))*0.25d0
         dvdxi  =  (w(3,i+1,j+1) +w(3,i+1,j)&
                   -w(3,i-1,j+1) -w(3,i-1,j))*0.25d0
         dTdxi  = (Tdeg(i+1,j+1)+Tdeg(i+1,j)&
                  -Tdeg(i-1,j+1)-Tdeg(i-1,j))*0.25d0
         do k=1,nY
            dwdxi(4+k)  =  (w(4+k,i+1,j+1) +w(4+k,i+1,j)&
                           -w(4+k,i-1,j+1) -w(4+k,i-1,j))*0.25d0
         end do

         dudeta =  w(2,i,j+1) -w(2,i,j)
         dvdeta =  w(3,i,j+1) -w(3,i,j)
         dTdeta = Tdeg(i,j+1)-Tdeg(i,j)
         do k=1,nY
            dwdeta(4+k) =  w(4+k,i,j+1) -w(4+k,i,j)
         end do

         dudx=dudxi*geojacj(1,1,i,j)+dudeta*geojacj(2,1,i,j)
         dvdx=dvdxi*geojacj(1,1,i,j)+dvdeta*geojacj(2,1,i,j)
         dTdx=dTdxi*geojacj(1,1,i,j)+dTdeta*geojacj(2,1,i,j)
         dudr=dudxi*geojacj(1,2,i,j)+dudeta*geojacj(2,2,i,j)
         dvdr=dvdxi*geojacj(1,2,i,j)+dvdeta*geojacj(2,2,i,j)
         dTdr=dTdxi*geojacj(1,2,i,j)+dTdeta*geojacj(2,2,i,j)
         hidwdx=0d0;hidwdr=0d0
         do k=1,nY
            dwdx(4+k)=dwdxi(4+k)*geojacj(1,1,i,j)+dwdeta(4+k)*geojacj(2,1,i,j)
            hidwdx=hidwdx+dwdx(4+k)*lhi(k)
            dwdr(4+k)=dwdxi(4+k)*geojacj(1,2,i,j)+dwdeta(4+k)*geojacj(2,2,i,j)
            hidwdr=hidwdr+dwdr(4+k)*lhi(k)
         end do

         !!2 dimentional-plane
         divu=dudx+dvdr

         t_xx=-2d0/3d0*mu*divu+2d0*mu*dudx
         t_rr=-2d0/3d0*mu*divu+2d0*mu*dvdr
         t_xr= mu*(dvdx+dudr)
         qx  = -k_heat*dTdx-rho*D*hidwdx
         qr  = -k_heat*dTdr-rho*D*hidwdr

         do k=1,nY
            Ev(k)=rho*D*dwdx(4+k)
         end do
         Ev(nY+1)=t_xx
         Ev(nY+2)=t_xr
         Ev(nY+3)=u*t_xx+v*t_xr-qx

         do k=1,nY
            Fv(k)=rho*D*dwdr(4+k)
         end do
         Fv(nY+1)=t_xr
         Fv(nY+2)=t_rr
         Fv(nY+3)=u*t_xr+v*t_rr-qr

         TGvj(:,i,j)=Ev*vnj(1,i,j)+Fv*vnj(2,i,j)
      end do
   end do
   !$omp end parallel do
   !}}}
end subroutine set_TGv!}}}

