module grbl_prmtr
   implicit none
   !!!! BELOW ARE USER DEFINISIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer,parameter::Nplane = 1
   integer,parameter::ni(Nplane)=(/1/)  !the number of split of i-direction
   integer,parameter::nj(Nplane)=(/1/)  !the number of split of j-direction
   integer,parameter::nimax=1
   integer,parameter::njmax=1
   integer,parameter::nY=2
   integer,parameter::nV=158
   character*100,parameter::file_coordinate="grid/dummy"
   !!!! ABOVE ARE USER DEFINISIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer,parameter::dimq  =nY+3
   integer,parameter::dimw  =dimq+5
   integer,parameter::indxg =dimq+2
   integer,parameter::indxht=dimq+3
   integer,parameter::indxR =dimq+4
   integer,parameter::indxMu=dimq+5
end module grbl_prmtr

