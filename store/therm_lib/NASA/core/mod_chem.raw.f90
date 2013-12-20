module const_chem
   implicit none
   integer,parameter::ne     = NE !the number of elements
   integer,parameter::max_ns = NS !the maximum number of species

   double precision,parameter::Ru         = 8.3144621d0 !universal gas constant
   double precision,parameter::pst        = 1d5         !standard state pressure
   double precision,parameter::omega      = 0.5d0       !relaxation factor
   double precision,parameter::eps        = 1d-9        !epsilon for convergence
   double precision,parameter::initial_eps= 1d-5        !epsilon to stabilize calculation
   double precision,parameter::Y_eps      = 1d-6        !epsilon for reduction determination
   double precision,parameter::TSIZE      = 1d-11       !epsilon for calc E or H
   double precision,parameter::TTSIZE     = 1d-14       !epsilon for calc n
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
   double precision rhof,pf,Tf,Ef,Hf,MWf,kappaf,muf
   double precision Yvf(nV),vhif(nV)
   double precision b0f(ne+2)
   double precision nf(max_ns)
   double precision nfini(max_ns)
   integer          nef
   integer          elistf(ne+2)
   integer          nelistf(ne+2)
   double precision maskf(max_ns)
   double precision maskbf(ne+2)

   double precision qo(dimq),wo(dimw)
   double precision rhoo,po,To,Eo,Ho,MWo,kappao,muo
   double precision Yvo(nV),vhio(nV)
   double precision b0o(ne+2)
   double precision no(max_ns)
   double precision noini(max_ns)
   integer          neo
   integer          elisto(ne+2)
   integer          nelisto(ne+2)
   double precision masko(max_ns)
   double precision maskbo(ne+2)

   !for flame sheet model
   double precision np(max_ns)
   double precision of
end module chem_var
