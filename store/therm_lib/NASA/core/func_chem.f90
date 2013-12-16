module func_therm
   use chem
contains
double precision function sr(species,T)!{{{
   implicit none
   integer,intent(in)::species
   double precision,intent(in)::T
   
   integer num_section
   integer l

   !consult section number{{{
   num_section=-1
   if (Trange(1,1,species) > T) then
      print *,"Error: T is out of range",T
      call exit(0)
   end if

   do l=1,num_sctn(species)
      if(Trange(2,l,species)>T) then
         num_section=l
         exit
      end if
   end do

   if(num_section .eq. -1) then
      print *,"Error: T is out of range",T
      call exit(0)
   end if
   !}}}

   sr=-co(1,num_section,species)*T**(-2)/2d0&
      -co(2,num_section,species)*T**(-1)&
      +co(3,num_section,species)*log(T)&
      +co(4,num_section,species)*T   &
      +co(5,num_section,species)*T**2/2d0&
      +co(6,num_section,species)*T**3/3d0&
      +co(7,num_section,species)*T**4/4d0&
      +co(9,num_section,species)
end function sr!}}}
subroutine check_section_number(species,T,num_section)!{{{
   implicit none
   integer,         intent(in):: species
   double precision,intent(in):: T
   integer,         intent(out)::num_section

   integer l

   num_section=-1
   if (Trange(1,1,species) > T) then
      num_section=1
   else
   do l=1,num_sctn(species)
      if(Trange(2,l,species)>T) then
         num_section=l
         exit
      end if
   end do
   end if

   if(num_section .eq. -1) num_section=num_sctn(species)
end subroutine check_section_number!}}}

subroutine calc_vThrt(T,logT,vT)!{{{
   implicit none
   double precision,intent(in)::T
   double precision,intent(in)::logT
   double precision,intent(out)::vT(9)

   double precision TT

   TT = 1d0/T**2

   vT(1) = -TT
   TT=TT*T
   vT(2) =  TT*logT
   vT(8) =  TT
   TT=TT*T
   vT(3) =  TT
   TT=TT*T
   vT(4) =  TT/2d0
   TT=TT*T
   vT(5) =  TT/3d0
   TT=TT*T
   vT(6) =  TT/4d0
   TT=TT*T
   vT(7) =  TT/5d0

   vT(9) =  0d0
 end subroutine calc_vThrt!}}}
subroutine calc_vTcpr(T,vT)!{{{
   implicit none
   double precision,intent(in)::T
   double precision,intent(out)::vT(9)

   double precision TT

   TT=1d0/T**2

   vT(1)=TT
   TT=TT*T
   vT(2)=TT
   TT=TT*T
   vT(3)=TT
   TT=TT*T
   vT(4)=TT
   TT=TT*T
   vT(5)=TT
   TT=TT*T
   vT(6)=TT
   TT=TT*T
   vT(7)=TT

   vT(8)=0d0
   vT(9)=0d0
end subroutine calc_vTcpr!}}}
subroutine calc_vTmurt(T,logT,vT)!{{{
   implicit none
   double precision,intent(in)::T
   double precision,intent(in)::logT
   double precision,intent(out)::vT(9)

   double precision TT

   TT = 1d0/T**2

   vT(1) = -TT/2d0
   TT=TT*T
   vT(2) =  TT*(1d0+logT)
   vT(8) =  TT
   TT=TT*T
   vT(3) =  TT*(1d0-logT)
   vT(9) = -TT
   TT=TT*T
   vT(4) = -TT/2d0
   TT=TT*T
   vT(5) = -TT/6d0
   TT=TT*T
   vT(6) = -TT/12d0
   TT=TT*T
   vT(7) = -TT/20d0
end subroutine calc_vTmurt!}}}
subroutine calc_vTmu(T,logT,vT)!{{{
   implicit none
   double precision,intent(in)::T
   double precision,intent(in)::logT
   double precision,intent(out)::vT(4)

   double precision TT


   TT = 1d0/T**2
   vT(1) =  logT
   vT(2) =  TT*T
   vT(3) =  TT
   vT(4) =  1d0
end subroutine calc_vTmu!}}}

double precision function calc_mu(species,T,logT,num_section)!{{{
   implicit none
   integer,intent(in)::species
   double precision,intent(in)::T
   double precision,intent(in)::logT
   integer,intent(in)::num_section
  
   calc_mu= trans(1,num_section,species)*logT&
          +(trans(2,num_section,species)*T&
          + trans(3,num_section,species))/(T**2)&
          + trans(4,num_section,species)

   calc_mu=exp(calc_mu)
end function calc_mu!}}}
subroutine check_section_number_trans(species,T,num_section)!{{{
   implicit none
   integer,         intent(in):: species
   double precision,intent(in):: T
   integer,         intent(out)::num_section

   integer l

   num_section=-1
   if (Trange_trans(1,1,species) > T) num_section = 1

   do l=1,num_sctn_trans(species)
      if(Trange_trans(2,l,species)>T) then
         num_section=l
         exit
      end if
   end do

   if(num_section .eq. -1) num_section=num_sctn_trans(species)
end subroutine check_section_number_trans!}}}
end module func_therm
