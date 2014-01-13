subroutine set_geojac!{{{
   use grbl_prmtr
   use mod_mpi
   use variable
   implicit none
   double precision dxdxi,dydxi,dxdeta,dydeta
   double precision det
   integer i,j,plane

   ! set geojaci
   do plane=nps,npe
      do j=nys(plane),nye(plane)
         do i=nxs(plane)-1,nxe(plane)
            dxdxi  = x(i+1,j,plane)-x(i,j,plane)
            dydxi  = r(i+1,j,plane)-r(i,j,plane)

            dxdeta = (x(i+1,j+1,plane)+x(i,j+1,plane)&
                     -x(i+1,j-1,plane)-x(i,j-1,plane))*0.25d0
            dydeta = (r(i+1,j+1,plane)+r(i,j+1,plane)&
                     -r(i+1,j-1,plane)-r(i,j-1,plane))*0.25d0

            det= dxdxi*dydeta - dydxi*dxdeta

            geojaci(1,1,i,j,plane) = dydeta/det!d( xi)/d(x)
            geojaci(1,2,i,j,plane) =-dxdeta/det!d( xi)/d(y)
            geojaci(2,1,i,j,plane) =-dydxi /det!d(eta)/d(x)
            geojaci(2,2,i,j,plane) = dxdxi /det!d(eta)/d(y)
         end do
      end do
   end do

   ! set geojacj
   do plane=nps,npe
      do j=nys(plane)-1,nye(plane)
         do i=nxs(plane),nxe(plane)
            dxdxi = (x(i+1,j+1,plane)+x(i+1,j,plane)&
                    -x(i-1,j+1,plane)-x(i-1,j,plane))*0.25d0
            dydxi = (r(i+1,j+1,plane)+r(i+1,j,plane)&
                    -r(i-1,j+1,plane)-r(i-1,j,plane))*0.25d0

            dxdeta= x(i,j+1,plane)-x(i,j,plane)
            dydeta= r(i,j+1,plane)-r(i,j,plane)


            det= dxdxi*dydeta - dydxi*dxdeta

            geojacj(1,1,i,j,plane) = dydeta/det!d( xi)/d(x)
            geojacj(1,2,i,j,plane) =-dxdeta/det!d( xi)/d(y)
            geojacj(2,1,i,j,plane) =-dydxi /det!d(eta)/d(x)
            geojacj(2,2,i,j,plane) = dxdxi /det!d(eta)/d(y)
         end do
      end do
   end do
end subroutine set_geojac!}}}

subroutine set_TGv!{{{
   use grbl_prmtr
   use mod_mpi
   use variable
   implicit none
   double precision mu,k_heat
   double precision,dimension(dimq)::Ev,Fv
   double precision,dimension(dimw)::dwdxi,dwdeta,dwdx,dwdr
   double precision,dimension(nV)::dYvdxi,dYvdeta
   double precision,dimension(0:nimax+1,0:njmax+1)::Tdeg
   double precision dudxi,dudeta
   double precision dudx, dudr
   double precision dvdxi,dvdeta
   double precision dvdx, dvdr
   double precision dTdxi,dTdeta
   double precision dTdx, dTdr
   double precision dYvdx,dYvdr
   double precision rho,u,v,T
   double precision Rgas,gamm
   double precision t_xx,t_xr,t_rr
   double precision qx,qr
   double precision divu
   double precision lhi(nV),hidYvdx,hidYvdr

   double precision D
   double precision,parameter::Sc=1d0
   double precision,parameter::Pr=1d0
   integer i,j,k,plane

   do plane=nps,npe
      !set Tdeg{{{
      !$omp parallel do private(i)
      do j=nys(plane)-1,nye(plane)+1
         do i=nxs(plane)-1,nxe(plane)+1
            Tdeg(i,j)= w(4,i,j,plane)/w(1,i,j,plane)/w(indxR,i,j,plane)
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
      !$omp                     dYvdxi,dYvdeta,dYvdx,dYvdr,&
      !$omp                     lhi,hidYvdx,hidYvdr,&
      !$omp                     divu,t_xx,t_xr,t_rr,qx,qr,&
      !$omp                     Ev,Fv)
      do j=nys(plane),nye(plane)
         do i=nxs(plane)-1,nxe(plane)
            rho = (w(1    ,i,j,plane)+w(1    ,i+1,j,plane))*0.5d0
            u   = (w(2    ,i,j,plane)+w(2    ,i+1,j,plane))*0.5d0
            v   = (w(3    ,i,j,plane)+w(3    ,i+1,j,plane))*0.5d0
            gamm= (w(indxg,i,j,plane)+w(indxg,i+1,j,plane))*0.5d0
            T   =    (Tdeg(i,j)      +   Tdeg(i+1,j)      )*0.5d0 !total energy
            Rgas=   (w(indxR,i,j,plane)+  w(indxR,i+1,j,plane))*0.5d0
            mu  =  (w(indxMu,i,j,plane)+ w(indxMu,i+1,j,plane))*0.5d0
            do k=1,nV
               lhi(k)=(vhi(k,i,j,plane) +vhi(  k,i+1,j,plane))*0.5d0
            end do

            D     =mu/rho/Sc
            k_heat=mu*gamm/(gamm-1d0)*Rgas/Pr               !k from mu by Prandtl number

            dudxi  =  w(2,i+1,j,plane) -w(2,i,j,plane)
            dvdxi  =  w(3,i+1,j,plane) -w(3,i,j,plane)
            dTdxi  = Tdeg(i+1,j)      -Tdeg(i,j)
            do k=1,nY
               dwdxi(4+k)  =  w(4+k,i+1,j,plane) -w(4+k,i,j,plane)
            end do
            do k=1,nV
               dYvdxi(k)  =  Yv(k,i+1,j,plane) -Yv(k,i,j,plane)
            end do

            dudeta =  (w(2,i+1,j+1,plane) +w(2,i,j+1,plane)&
                      -w(2,i+1,j-1,plane) -w(2,i,j-1,plane))*0.25d0
            dvdeta =  (w(3,i+1,j+1,plane) +w(3,i,j+1,plane)&
                      -w(3,i+1,j-1,plane) -w(3,i,j-1,plane))*0.25d0
            dTdeta = (Tdeg(i+1,j+1)      +Tdeg(i,j+1)     &
                     -Tdeg(i+1,j-1)      -Tdeg(i,j-1)     )*0.25d0
            do k=1,nY
               dwdeta(4+k) =  (w(4+k,i+1,j+1,plane) +w(4+k,i,j+1,plane)&
                              -w(4+k,i+1,j-1,plane) -w(4+k,i,j-1,plane))*0.25d0
            end do
            do k=1,nV
               dYvdeta(k) =  (Yv(k,i+1,j+1,plane) +Yv(k,i,j+1,plane)&
                             -Yv(k,i+1,j-1,plane) -Yv(k,i,j-1,plane))*0.25d0
            end do

            dudx=dudxi*geojaci(1,1,i,j,plane)+dudeta*geojaci(2,1,i,j,plane)
            dvdx=dvdxi*geojaci(1,1,i,j,plane)+dvdeta*geojaci(2,1,i,j,plane)
            dTdx=dTdxi*geojaci(1,1,i,j,plane)+dTdeta*geojaci(2,1,i,j,plane)
            dudr=dudxi*geojaci(1,2,i,j,plane)+dudeta*geojaci(2,2,i,j,plane)
            dvdr=dvdxi*geojaci(1,2,i,j,plane)+dvdeta*geojaci(2,2,i,j,plane)
            dTdr=dTdxi*geojaci(1,2,i,j,plane)+dTdeta*geojaci(2,2,i,j,plane)
            do k=1,nY
               dwdx(4+k)=dwdxi(4+k)*geojaci(1,1,i,j,plane)+dwdeta(4+k)*geojaci(2,1,i,j,plane)
               dwdr(4+k)=dwdxi(4+k)*geojaci(1,2,i,j,plane)+dwdeta(4+k)*geojaci(2,2,i,j,plane)
            end do
            hidYvdx=0d0;hidYvdr=0d0
            do k=1,nV
               dYvdx=dYvdxi(k)*geojaci(1,1,i,j,plane)+dYvdeta(k)*geojaci(2,1,i,j,plane)
               hidYvdx=hidYvdx+dYvdx*lhi(k)
               dYvdr=dYvdxi(k)*geojaci(1,2,i,j,plane)+dYvdeta(k)*geojaci(2,2,i,j,plane)
               hidYvdr=hidYvdr+dYvdr*lhi(k)
            end do

            !!!2 dimentional-plane
            divu=dudx+dvdr

            t_xx=-2d0/3d0*mu*divu+2d0*mu*dudx
            t_rr=-2d0/3d0*mu*divu+2d0*mu*dvdr
            t_xr= mu*(dvdx+dudr)
            qx  = -k_heat*dTdx-rho*D*hidYvdx
            qr  = -k_heat*dTdr-rho*D*hidYvdr

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

            TGvi(:,i,j,plane)=Ev*vni(1,i,j,plane)+Fv*vni(2,i,j,plane)
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
      !$omp                     dYvdxi,dYvdeta,dYvdx,dYvdr,&
      !$omp                     lhi,hidYvdx,hidYvdr,&
      !$omp                     divu,t_xx,t_xr,t_rr,qx,qr,&
      !$omp                     Ev,Fv)
      do j=nys(plane)-1,nye(plane)
         do i=nxs(plane),nxe(plane)
            rho = (w(1    ,i,j,plane)+w(1    ,i,j+1,plane))*0.5d0
            u   = (w(2    ,i,j,plane)+w(2    ,i,j+1,plane))*0.5d0
            v   = (w(3    ,i,j,plane)+w(3    ,i,j+1,plane))*0.5d0
            gamm= (w(indxg,i,j,plane)+w(indxg,i,j+1,plane))*0.5d0
            T   =    (Tdeg(i,j)      +   Tdeg(i,j+1)      )*0.5d0
            Rgas=   (w(indxR,i,j,plane)+  w(indxR,i,j+1,plane))*0.5d0
            mu  =  (w(indxMu,i,j,plane)+ w(indxMu,i,j+1,plane))*0.5d0
            do k=1,nV
               lhi(k)=(vhi(k,i,j,plane) +vhi(  k,i,j+1,plane))*0.5d0
            end do

            D     =mu/rho/Sc
            k_heat=mu*gamm/(gamm-1d0)*Rgas/Pr               !k from mu by Prandtl number

            dudxi  =  (w(2,i+1,j+1,plane) +w(2,i+1,j,plane)&
                      -w(2,i-1,j+1,plane) -w(2,i-1,j,plane))*0.25d0
            dvdxi  =  (w(3,i+1,j+1,plane) +w(3,i+1,j,plane)&
                      -w(3,i-1,j+1,plane) -w(3,i-1,j,plane))*0.25d0
            dTdxi  = (Tdeg(i+1,j+1)      +Tdeg(i+1,j)     &
                     -Tdeg(i-1,j+1)      -Tdeg(i-1,j)     )*0.25d0
            do k=1,nY
               dwdxi(4+k)  =  (w(4+k,i+1,j+1,plane) +w(4+k,i+1,j,plane)&
                              -w(4+k,i-1,j+1,plane) -w(4+k,i-1,j,plane))*0.25d0
            end do
            do k=1,nV
               dYvdxi(k)  =  (Yv(k,i+1,j+1,plane) +Yv(k,i+1,j,plane)&
                             -Yv(k,i-1,j+1,plane) -Yv(k,i-1,j,plane))*0.25d0
            end do

            dudeta =  w(2,i,j+1,plane) -w(2,i,j,plane)
            dvdeta =  w(3,i,j+1,plane) -w(3,i,j,plane)
            dTdeta = Tdeg(i,j+1)      -Tdeg(i,j)
            do k=1,nY
               dwdeta(4+k) =  w(4+k,i,j+1,plane) -w(4+k,i,j,plane)
            end do
            do k=1,nV
               dYvdeta(k) =  Yv(k,i,j+1,plane) -Yv(k,i,j,plane)
            end do

            dudx=dudxi*geojacj(1,1,i,j,plane)+dudeta*geojacj(2,1,i,j,plane)
            dvdx=dvdxi*geojacj(1,1,i,j,plane)+dvdeta*geojacj(2,1,i,j,plane)
            dTdx=dTdxi*geojacj(1,1,i,j,plane)+dTdeta*geojacj(2,1,i,j,plane)
            dudr=dudxi*geojacj(1,2,i,j,plane)+dudeta*geojacj(2,2,i,j,plane)
            dvdr=dvdxi*geojacj(1,2,i,j,plane)+dvdeta*geojacj(2,2,i,j,plane)
            dTdr=dTdxi*geojacj(1,2,i,j,plane)+dTdeta*geojacj(2,2,i,j,plane)
            do k=1,nY
               dwdx(4+k)=dwdxi(4+k)*geojacj(1,1,i,j,plane)+dwdeta(4+k)*geojacj(2,1,i,j,plane)
               dwdr(4+k)=dwdxi(4+k)*geojacj(1,2,i,j,plane)+dwdeta(4+k)*geojacj(2,2,i,j,plane)
            end do
            hidYvdx=0d0;hidYvdr=0d0
            do k=1,nV
               dYvdx=dYvdxi(k)*geojacj(1,1,i,j,plane)+dYvdeta(k)*geojacj(2,1,i,j,plane)
               hidYvdx=hidYvdx+dYvdx*lhi(k)
               dYvdr=dYvdxi(k)*geojacj(1,2,i,j,plane)+dYvdeta(k)*geojacj(2,2,i,j,plane)
               hidYvdr=hidYvdr+dYvdr*lhi(k)
            end do

            !!2 dimentional-plane
            divu=dudx+dvdr

            t_xx=-2d0/3d0*mu*divu+2d0*mu*dudx
            t_rr=-2d0/3d0*mu*divu+2d0*mu*dvdr
            t_xr= mu*(dvdx+dudr)
            qx  = -k_heat*dTdx-rho*D*hidYvdx
            qr  = -k_heat*dTdr-rho*D*hidYvdr

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

            TGvj(:,i,j,plane)=Ev*vnj(1,i,j,plane)+Fv*vnj(2,i,j,plane)
         end do
      end do
      !$omp end parallel do
      !}}}
   end do
end subroutine set_TGv!}}}

subroutine init_Svq
   use variable
   implicit none
   Svq(:,:,:,:)=0d0
end subroutine init_Svq
