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
  
   !!!!!!!!!!!!!!!!! CHANGE BELOW IF YOU CHANGE CHEM.INP !!!!!!!!!!!!!!!!!!!!!
   integer,parameter::ne=NumOfElements  !num of elements
   integer,parameter::ns=NumOfSpecies   !num of species
   integer,parameter::nr=NumOfReactions !num of reactions

   integer,parameter::LRW=22+9*(ns+1)+2*(ns+1)**2
   integer,parameter::LIW=30+  (ns+1)

   integer,parameter::max_recalc=50
   !!!!!!!!!!!!!!!!! CHANGE ABOVE IF YOU CHANGE CHEM.INP !!!!!!!!!!!!!!!!!!!!!
end module const_chem

module chem
   use const_chem
   implicit none
   integer nt,ns_tocalc

   character(3), dimension(ne)::SYM_ELM
   character(18),dimension(ns)::SYM_SPC
   character(80),dimension(nr)::SYM_RCT

   double precision    ES(ne,ns)

   double precision Tthre(3,  ns)
   double precision coeff(7,2,ns)
   double precision   MWs(    ns)
   double precision invMW(    ns)

   integer          Rstate(nr,2) !(,1)...0:bi-direction;1:mono-direction
                                      !(,2)...0:none;1:rev;2:low;3:troe;4:+M
   integer          IndNu(3,2,nr)
   integer          NumNu(  2,nr)
   integer          snu(      nr)

   logical          exist_M(nr)
   integer          IndM(ns,nr)
   integer          NumM(   nr)
   double precision Men( ns,nr)

   double precision    ABE(3,nr)
   double precision   cABE(3,nr)
   double precision   TROE(4,nr)

   character(18),dimension(ns)::species_name_trans
   integer,dimension(ns)::num_sctn_trans,tr2th
   double precision Trange_trans(2,5,ns)
   double precision trans(4,5,ns)
end module chem

module chem_var
   use const_chem
   implicit none

   double precision     pf
   double precision   rhof
   double precision     Tf
   double precision     Ef
   double precision  vrhof(ns)
   double precision    vwf(ns)
   double precision   vhif(ns)
   double precision    MWf
   double precision kappaf
   double precision    muf

   double precision     po
   double precision   rhoo
   double precision     To
   double precision     Eo
   double precision  vrhoo(ns)
   double precision    vwo(ns)
   double precision   vhio(ns)
   double precision    MWo
   double precision kappao
   double precision    muo
end module chem_var
