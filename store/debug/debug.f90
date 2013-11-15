subroutine debug_ABpm
   use grbl_prmtr
   use half_value
   use variable
   use var_lusgs
   implicit none
   double precision q_org(1:dimq,-1:ni+2,-1:nj+2)
   double precision,dimension(1:dimq)::q_tmp,Eu,El,Fu,Fl
   double precision,dimension(1:dimq,1:dimq)::A_diff,B_diff
   double precision delta
   integer i,j,k,l,m

   q_org = q

   i = 88
   j = 50
   do l=1,dimq
      q_tmp = q_org(:,i,j)
      q_tmp(l) = q_tmp(l)*1.01d0
      delta = q_tmp(l)

      q(:,i  ,j  ) = q_tmp
      q(:,i+1,j  ) = q_tmp
      q(:,i  ,j+1) = q_tmp
      call set_w
      call set_thermo_prop
      call set_HV(w,wHli,wHri,wHlj,wHrj)
      call set_TG

      Eu(:)=TGi(:,i,j)
      Fu(:)=TGj(:,i,j)

      q_tmp = q_org(:,i,j)
      q_tmp(l) = q_tmp(l)*0.99d0
      delta = delta-q_tmp(l)

      q(:,i  ,j  ) = q_tmp
      q(:,i+1,j  ) = q_tmp
      q(:,i  ,j+1) = q_tmp
      call set_w
      call set_thermo_prop
      call set_HV(w,wHli,wHri,wHlj,wHrj)
      call set_TG

      El(:)=TGi(:,i,j)
      !print '(2es9.2)',(q_tmp(m),El(m),Fl(m),m=1,dimq)
      Fl(:)=TGj(:,i,j)
      do m=1,dimq
         !print '(2i3,2es9.1)',m,l,Ap(m,l,i,j)+Am(m,l,i,j),(Eu(m)-El(m))/delta
         print '(2i3,2es9.1)',m,l,Bp(m,l,i,j)+Bm(m,l,i,j),(Fu(m)-Fl(m))/delta
      end do
      print '()'
   end do
   call exit(0)
end subroutine debug_ABpm
