module const_chem
   implicit none
   double precision,parameter::Ru=8.3144621d0   !universal gas constant
   double precision,parameter::pst=1d5          !standard state pressure
   double precision,parameter::omega=0.1d0      !relaxation factor
   double precision,parameter::eps=1d-6         !epsilon for convergence
   double precision,parameter::n_eps=1d-20      !epsilon for mol/kg
   double precision,parameter::initial_eps=1d-5  !epsilon for initial mole/mass fraction
   double precision,parameter::TSIZE=1d-11
   integer,parameter::ne=4     !the number of elements
   integer,parameter::max_ns=600  !the maximum number of species
   integer          ns
   integer          nt
end module const_chem

module chem
   use const_chem
   implicit none
   integer         ,dimension(max_ns)::num_sctn
   double precision,dimension(max_ns)::MW,DH0
   double precision,dimension(2,6,max_ns)::Trange
   double precision,dimension(9,6,max_ns)::co

   character*18,dimension(max_ns)::SYM_SPC
   character*2 ,dimension(ne)::    SYM_ELM

   !for trans
   double precision trans(4,3,max_ns),Trange_trans(2,3,max_ns)
   character*18     species_name_trans(max_ns)
   integer         ,dimension(max_ns)::num_sctn_trans
   integer         ,dimension(max_ns)::tr2th

   double precision,dimension(ne,max_ns)::Ac
end module chem

module chem_var
   use grbl_prmtr
   use const_chem
   implicit none
   double precision      n_save(1:max_ns,1:ni,1:nj)

   integer            numf
   double precision   rhof
   double precision     pf
   double precision     Tf
   double precision     Ef
   double precision    MWf
   double precision kappaf
   double precision    muf
   double precision    b0f(1:ne)
   double precision     nf(1:max_ns)
   character*18,    dimension(:),allocatable::SYM_FUEL
   double precision,dimension(:),allocatable::COMP_FUEL

   integer            numo
   double precision   rhoo
   double precision     po
   double precision     To
   double precision     Eo
   double precision    MWo
   double precision kappao
   double precision    muo
   double precision    b0o(1:ne)
   double precision     no(1:max_ns)
   character*18,    dimension(:),allocatable::SYM_OXID
   double precision,dimension(:),allocatable::COMP_OXID
end module chem_var
