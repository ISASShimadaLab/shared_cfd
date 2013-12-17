subroutine LU(A, b, x, n,flag)
   use const_chem
   implicit none
   double precision,intent(inout)::A(ne+2, ne+2)
   double precision,intent(inout)::b(ne+2)
   double precision,intent(out)::x(ne+2)
   integer,intent(in)::n
   logical,intent(out)::flag
   
   double precision c(n)
   integer i, j, k
   double precision dtemp,p
   integer ip

   !print *,"A="
   !print '(3es9.1)',A(1:3,1)
   !print '(3es9.1)',A(1:3,2)
   !print '(3es9.1)',A(1:3,3)
   flag = .false.
 
   !--- LU decomposition ----------------------
   do k=1, n
     !pibotting
     p=abs(A(k,k)); ip=k

     do i=k+1,n
        dtemp =abs(A(i,k))
        if(p<dtemp) then
           p =dtemp
           ip=i
        end if
     end do

     if(p<1d-30) then
        !print *,"Error! Singular Matrix! p=",p,"k=",k
        flag = .true.
        return
     end if

     if(ip .ne. k) then
        do i=1,n
           dtemp  =A( k,i)
           A(k,i) =A(ip,i)
           A(ip,i)=dtemp
        end do
        dtemp = b(k)
        b(k)  = b(ip)
        b(ip) = dtemp
     end if

     dtemp = 1.0d0 / A(k, k)
     do i=k+1, n
       A(i, k) = A(i, k)*dtemp
     enddo
     do j=k+1, n
       dtemp = A(k, j)
       do i=k+1, n
         A(i, j) = A(i, j) - dtemp * A(i, k)
       enddo
     enddo
   enddo
   !-------------------------------------------
   
   !--- Forward substitution ------------------ 
   do k=1, n
     c(k) = b(k)
     do j=1, k-1
       c(k) = c(k) - A(k, j)*c(j)
     enddo
   enddo
   !--------------------------------------------
   
   !--- Backward substitution ------------------  
   x(n) = c(n) / A(n, n)
   do k=n-1, 1, -1
     x(k) = c(k)
     do j=k+1, n
       x(k) = x(k) - A(k, j)*x(j)
     enddo
     x(k) = x(k) / A(k, k)
   enddo
   !-------------------------------------------
end subroutine LU
