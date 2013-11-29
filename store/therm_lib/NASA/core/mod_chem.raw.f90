module const_chem
   implicit none
   integer,parameter::ne     = NE !the number of elements
   integer,parameter::max_ns = NS !the maximum number of species

   double precision,parameter::Ru         = 8.3144621d0 !universal gas constant
   double precision,parameter::pst        = 1d5         !standard state pressure
   double precision,parameter::omega      = 0.5d0       !relaxation factor
   double precision,parameter::eps        = 1d-8        !epsilon for convergence
   double precision,parameter::n_eps      = 1d-20       !epsilon for mol/kg
   double precision,parameter::initial_eps= 1d-5        !epsilon for initial mole/mass fraction
   double precision,parameter::TSIZE      = 1d-11
   integer          ns
   integer          nt
end module const_chem

module chem
   use const_chem
   implicit none
   integer         ,dimension(max_ns)::num_sctn
   double precision,dimension(max_ns)::MW
   double precision,dimension(2,6,max_ns)::Trange
   double precision,dimension(9,6,max_ns)::co

   character*2 ,dimension(ne)::SYM_ELM
   character*18,dimension(max_ns)::SYM_SPC

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
   double precision n_save(max_ns,nimax,njmax,Nplane)

   double precision qf(dimq),wf(dimw)
   double precision rhof,pf,Tf,Ef,MWf,kappaf,muf
   double precision Yvf(nV),vhif(nV)
   double precision b0f(ne)
   double precision nf(max_ns)

   double precision qo(dimq),wo(dimw)
   double precision rhoo,po,To,Eo,MWo,kappao,muo
   double precision Yvo(nV),vhio(nV)
   double precision b0o(ne)
   double precision no(max_ns)

   !for flame sheet model
   double precision np(max_ns)
   double precision of
end module chem_var
