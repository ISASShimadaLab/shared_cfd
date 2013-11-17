module grbl_prmtr
   implicit none
   !!!! BELOW ARE USER DEFINISIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer,parameter::Nplane = 1
   integer,parameter::ni(Nplane)=(/223/)  !the number of split of i-direction
   integer,parameter::nj(Nplane)=(/99/)  !the number of split of j-direction
   integer,parameter::nimax=223
   integer,parameter::njmax=99
   integer,parameter::nY=1
   character*100,parameter::file_coordinate="grid/ideal_bl.x"
   !!!! ABOVE ARE USER DEFINISIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer,parameter::dimq  =nY+3
   integer,parameter::dimw  =dimq+5
   integer,parameter::indxg =dimq+2
   integer,parameter::indxht=dimq+3
   integer,parameter::indxR =dimq+4
   integer,parameter::indxMu=dimq+5
end module grbl_prmtr

