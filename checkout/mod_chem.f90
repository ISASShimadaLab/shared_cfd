module const_chem
   implicit none
   integer,parameter::ne=3 !the number of elements
   integer,parameter::ns=4 !the number of species

   double precision,parameter::Ru         = 8.3144621d0 !universal gas constant
   double precision,parameter::pst        = 1d5         !standard state pressure
   double precision,parameter::omega      = 0.1d0       !relaxation factor
   double precision,parameter::eps        = 1d-6        !epsilon for convergence
   double precision,parameter::n_eps      = 1d-20       !epsilon for mol/kg
   double precision,parameter::initial_eps= 1d-5        !epsilon for initial mole/mass fraction
   double precision,parameter::TSIZE      = 1d-11
   integer          nt
end module const_chem

module chem
   use const_chem
   implicit none
   integer         ,dimension(ns)::num_sctn
   double precision,dimension(ns)::MW
   double precision,dimension(2,6,ns)::Trange
   double precision,dimension(9,6,ns)::co

   character*2 ,dimension(ne)::SYM_ELM
   character*18,dimension(ns)::SYM_SPC

   !for trans
   double precision trans(4,3,ns),Trange_trans(2,3,ns)
   character*18     species_name_trans(ns)
   integer         ,dimension(ns)::num_sctn_trans
   integer         ,dimension(ns)::tr2th

   double precision,dimension(ne,ns)::Ac
end module chem

module chem_var
   use grbl_prmtr
   use const_chem
   implicit none
   double precision n_save(1:ns,1:ni,1:nj)

   double precision rhof,pf,Tf,Ef,MWf,kappaf,muf
   double precision b0f(1:ne)
   double precision nf(1:ns)

   double precision rhoo,po,To,Eo,MWo,kappao,muo
   double precision b0o(1:ne)
   double precision no(1:ns)

   double precision np(1:ns)
   double precision of
end module chem_var
