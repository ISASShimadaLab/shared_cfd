module const_chem
   implicit none
   double precision,parameter::erg2cal=1d0/4.184e7
   double precision,parameter::cal2erg=1d0/erg2cal
   double precision,parameter::b2Pa=0.1d0
   double precision,parameter::Pa2b=1d0/b2Pa

   double precision,parameter::Ru=8.3144621d7        !universal gas constant in cgs unit
   double precision,parameter::invRc=1d0/(1.987e0)
                                            !inverse of universal gas constant in mol/cal
   double precision,parameter::pst=1.01325d6             !standard state pressure(dyn/cm^2)
   double precision,parameter::logPR=- 4.407418526701446!log(Pst/Ru)

   double precision,parameter::n_eps = 1d-40            !epsilon for small n
   double precision rtol_user,atol_user
  
   !!!!!!!!!!!!!!!!! CHANGE BELOW IF YOU CHANGE CHEM.INP !!!!!!!!!!!!!!!!!!!!!
   logical,parameter::isColdFlow = ISCOLDFLOW
   integer,parameter::max_ne=NE !num of elements
   integer,parameter::max_ns=NS !num of species
   integer,parameter::max_nr=NR !num of reactions
   integer ne,ns,nr,num_reac

   integer,parameter::LRW=22+9*(max_ns+1)+2*(max_ns+1)**2
   integer,parameter::LIW=30+  (max_ns+1)

   integer,parameter::max_recalc=50
   !!!!!!!!!!!!!!!!! CHANGE ABOVE IF YOU CHANGE CHEM.INP !!!!!!!!!!!!!!!!!!!!!
end module const_chem

module chem
   use const_chem
   implicit none
   integer nt,ns_tocalc

   character(3), dimension(max_ne)::SYM_ELM
   character(18),dimension(max_ns)::SYM_SPC
   character(80),dimension(max_nr)::SYM_RCT

   double precision    ES(max_ne,max_ns)

   double precision Tthre(3,  max_ns)
   double precision coeff(7,2,max_ns)
   double precision   MWs(    max_ns)
   double precision invMW(    max_ns)

   integer          Rstate(max_nr,2) !(,1)...0:bi-direction;1:mono-direction
                                     !(,2)...0:none;1:rev;2:low;3:troe;4:+M
   integer          IndNu(3,2,max_nr)
   integer          NumNu(  2,max_nr)
   integer          snu(      max_nr)

   logical          duplicate(max_nr)

   logical          exist_M(    max_nr)
   integer          IndM(max_ns,max_nr)
   integer          NumM(       max_nr)
   double precision Men( max_ns,max_nr)

   double precision    ABE(3,max_nr)
   double precision   cABE(3,max_nr)
   double precision   TROE(4,max_nr)

   character(18),dimension(max_ns)::species_name_trans
   integer,dimension(max_ns)::num_sctn_trans,tr2th
   double precision Trange_trans(2,5,max_ns)
   double precision trans(       4,5,max_ns)
end module chem

module chem_var
   use const_chem
   use grbl_prmtr
   implicit none

   double precision     of

   double precision     pf
   double precision   rhof
   double precision     Tf
   double precision     Ef
   double precision  vrhof(max_ns)
   double precision    vwf(max_ns)
   double precision   vhif(max_ns)
   double precision    MWf
   double precision kappaf
   double precision    muf
   double precision     qf(dimq)
   double precision     wf(dimw)

   double precision     po
   double precision   rhoo
   double precision     To
   double precision     Eo
   double precision  vrhoo(max_ns)
   double precision    vwo(max_ns)
   double precision   vhio(max_ns)
   double precision    MWo
   double precision kappao
   double precision    muo
   double precision     qo(dimq)
   double precision     wo(dimw)
end module chem_var
