subroutine set_thermo_prop(step)
   use mod_mpi
   use variable
   use prmtr
   use chem_var
   implicit none
   integer,intent(in)::step
   double precision ei,T,MWave,kappa,mu
   integer i,j,plane

   double precision Tfs,Tcea
   double precision tYf_cfh_ok,tYo_cfh_ok
   double precision tYf_cfh_ng,tYo_cfh_ng
   double precision Nall,Nselect,tNall,tNselect

   if(mod(step,Dstep_cfh).eq. 0d0) then
      !!!!!!!!!!!!!!!!!!!!! CHECK THRESHOLDS !!!!!!!!!!!!!!!!!!!!!
      !initialize cfh parameters
      tYf_cfh_ng= 1d0
      tYo_cfh_ng= 1d0
      tYf_cfh_ok= 0d0
      tYo_cfh_ok= 0d0
      tNall   =0d0
      tNselect=0d0
      do plane = nps,npe
         !$omp parallel do private(i,ei,T,MWave,kappa,mu,Tfs,Tcea) reduction(+:tNall,tNselect) &
         !$omp                            reduction(min:tYf_cfh_ng,tYo_cfh_ng) &
         !$omp                            reduction(max:tYf_cfh_ok,tYo_cfh_ok)
         do j=nys(plane),nye(plane)
            do i=nxs(plane),nxe(plane)
               w(5,i,j,plane)=min(max(0d0,w(5,i,j,plane)),1d0)
               w(6,i,j,plane)=min(max(0d0,w(6,i,j,plane)),1d0)
               ei = q(nY+3,i,j,plane)/w(1,i,j,plane)-0.5d0*(w(2,i,j,plane)**2+w(3,i,j,plane)**2)
               T  = w(4,i,j,plane)/(w(1,i,j,plane)*w(indxR,i,j,plane))

               call flame_sheet_cfh(w(5:6,i,j,plane),ei,&
                        T,&
                        MWave,kappa,mu,DHi(:,i,j,plane),Yv(:,i,j,plane),vhi(:,i,j,plane))
               Tfs =T
               call cea(w(1,i,j,plane),w(5:6,i,j,plane),ei,&
                        T,n_save(:,i,j,plane),&
                        MWave,kappa,mu,Yv(:,i,j,plane),vhi(:,i,j,plane),DHi(:,i,j,plane))
               Tcea=T

               DHi(:,   i,j,plane) = 0d0
               w(4,     i,j,plane) = w(1,i,j,plane)*(R_uni/MWave)*T
               w(indxg, i,j,plane) = kappa
               w(indxht,i,j,plane) = (q(nY+3,i,j,plane)+w(4,i,j,plane))/w(1,i,j,plane)
               w(indxMu,i,j,plane) = mu
               w(indxR, i,j,plane) = R_uni/MWave

               !calc cfh parameters
               tNall=tNall+1d0
               if(abs(Tfs-Tcea)<Tdiff_cfh) then
                  tNselect=tNselect+1d0
                  if(w(5,i,j,plane)<0.5d0) then
                     tYf_cfh_ok=max(tYf_cfh_ok,w(5,i,j,plane))
                  else         
                     tYo_cfh_ok=max(tYo_cfh_ok,w(6,i,j,plane))
                  end if
               else
                  if(w(5,i,j,plane)<0.5d0) then
                     tYf_cfh_ng=min(tYf_cfh_ng,w(5,i,j,plane))
                  else 
                     tYo_cfh_ng=min(tYo_cfh_ng,w(6,i,j,plane))
                  end if
               end if
            end do
         end do
         !$omp end parallel do
      end do
     
      !Compute cfl parameters
      call MPI_Reduce(tNselect, Nselect,1,MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tNall,    Nall,   1,MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(tYf_cfh_ok,Yf_cfh,1,MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(tYo_cfh_ok,Yo_cfh,1,MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_WORLD,ierr)
      tYf_cfh_ng=min(tYf_cfh_ng,Yf_cfh)
      tYo_cfh_ng=min(tYo_cfh_ng,Yo_cfh)
      call MPI_Allreduce(tYf_cfh_ng,Yf_cfh,1,MPI_DOUBLE_PRECISION, MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(tYo_cfh_ng,Yo_cfh,1,MPI_DOUBLE_PRECISION, MPI_MIN,MPI_COMM_WORLD,ierr)
      if(myid .eq. 0) then
         print '(a,i12,x,a,f8.2,x,2(a,es9.1,x))',"step=",step,&
                                                  "fs/all(%)=",Nselect/Nall*1d2,&
                                                  "Yf_cfh=",Yf_cfh,&
                                                  "Yo_cfh=",Yo_cfh
      end if
   else
      !!!!!!!!!!!!!!!!!!!!! NORMAL OPERATION !!!!!!!!!!!!!!!!!!!!!
      do plane = nps,npe
         !$omp parallel do private(i,ei,T,MWave,kappa,mu)
         do j=nys(plane),nye(plane)
            do i=nxs(plane),nxe(plane)
               w(5,i,j,plane)=min(max(0d0,w(5,i,j,plane)),1d0)
               w(6,i,j,plane)=min(max(0d0,w(6,i,j,plane)),1d0)
               ei = q(nY+3,i,j,plane)/w(1,i,j,plane)-0.5d0*(w(2,i,j,plane)**2+w(3,i,j,plane)**2)
               T  = w(4,i,j,plane)/(w(1,i,j,plane)*w(indxR,i,j,plane))

               if(w(5,i,j,plane)<Yf_cfh .or. w(6,i,j,plane)<Yo_cfh) then
                  call flame_sheet_cfh(w(5:6,i,j,plane),ei,&
                        T,&
                        MWave,kappa,mu,DHi(:,i,j,plane),Yv(:,i,j,plane),vhi(:,i,j,plane))
               else
                 call cea(w(1,i,j,plane),w(5:6,i,j,plane),ei,&
                       T,n_save(:,i,j,plane),&
                       MWave,kappa,mu,Yv(:,i,j,plane),vhi(:,i,j,plane),DHi(:,i,j,plane))
               end if

               DHi(:,   i,j,plane) = 0d0
               w(4,     i,j,plane) = w(1,i,j,plane)*(R_uni/MWave)*T
               w(indxg, i,j,plane) = kappa
               w(indxht,i,j,plane) = (q(nY+3,i,j,plane)+w(4,i,j,plane))/w(1,i,j,plane)
               w(indxMu,i,j,plane) = mu
               w(indxR, i,j,plane) = R_uni/MWave
            end do
         end do
         !$omp end parallel do
      end do
   end if
end subroutine set_thermo_prop

subroutine YPT2w(Yf,p,T,wt,vhit)
   use grbl_prmtr
   use chem
   use prmtr
   use chem_var
   implicit none
   double precision,intent(in) ::Yf
   double precision,intent(in) ::p
   double precision,intent(in) ::T
   double precision,intent(out)::wt(dimw)
   double precision,intent(out)::vhit(nY)

   double precision Y(2),H
   double precision MWave,kappa,mu,rho
   double precision,dimension(max_ns)::n,Yv
   logical flag_outer


   Y(1)=Yf
   Y(2)=1d0-Yf

   n=sum(no(1:ns)+nf(1:ns),1)/ns
   call cea_tp(p,Y,T, n, H,MWave,kappa,mu,Yv,vhit,flag_outer)
   if(.not.flag_outer) stop "Not converged at YPT2w"

   rho = P/(R_uni/MWave*T)
   wt(1)=rho
   wt(2)=0d0
   wt(3)=0d0
   wt(4)=p
   wt(5)=Y(1)
   wt(6)=Y(2)
   wt(indxg ) = kappa
   wt(indxht) = H
   wt(indxR ) = R_uni/MWave
   wt(indxMu) = mu
end subroutine YPT2w

