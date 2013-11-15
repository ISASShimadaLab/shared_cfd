module grbl_prmtr
   implicit none
   !!!! BELOW ARE USER DEFINISIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer,parameter::Nplane = NPLANE
   integer,parameter::ni(Nplane)=(/NumI/)  !the number of split of i-direction
   integer,parameter::nj(Nplane)=(/NumJ/)  !the number of split of j-direction
   integer,parameter::nimax=NIMAX
   integer,parameter::njmax=NJMAX
   integer,parameter::nY=NumY
   character*100,parameter::file_coordinate="grid/GridFileName"
   !!!! ABOVE ARE USER DEFINISIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer,parameter::dimq  =nY+3
   integer,parameter::dimw  =dimq+5
   integer,parameter::indxg =dimq+2
   integer,parameter::indxht=dimq+3
   integer,parameter::indxR =dimq+4
   integer,parameter::indxMu=dimq+5
end module grbl_prmtr

