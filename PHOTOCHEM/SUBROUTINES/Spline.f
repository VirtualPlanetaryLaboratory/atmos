      SUBROUTINE spline(x,y,n,yp1,ypn,y2) 
!from Numerical recipes in F77

      INTEGER n,NMAX 
      REAL*8 yp1,ypn,x(n),y(n),y2(n) 
      PARAMETER (NMAX=500)
!Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with x1 < x2 < ... < xN, and given values yp1 and ypn for the first derivative of the inter- polating function at points 1 and n, respectively, this routine returns an array y2(1:n) of length n which contains the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set the corresponding boundary condition for a natural spline, with zero second derivative on that boundary.
!Parameter: NMAX is the largest anticipated value of n. INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      
      if (yp1.gt..99e30) then !The lower boundary condition is set either to be “natural” 
         y2(1)=0.
         u(1)=0. 
      else                    !or else to have a specified first derivative.
         y2(1)=-0.5
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      endif

!This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors.
      do i=2,n-1 
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1)) 
         p=sig*y2(i-1)+2. 
         y2(i)=(sig-1.)/p 
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) / 
     $        (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo

      if (ypn.gt..99e30) then  !The upper boundary condition is set either to be “natural”
         qn=0.
         un=0. 
      else                     !or else to have a specified first derivative.
         qn=0.5
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))) 
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k) 
      enddo

      return 
      END


      SUBROUTINE splint(xa,ya,y2a,n,x,y) 
      INTEGER n 
      REAL*8 x,y,xa(n),y2a(n),ya(n)
!Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xa's in order), and given the array y2a(1:n), which is the output from spline above, and given a value of x, this routine returns a cubic-spline interpolated value y.
      INTEGER k,khi,klo 
      REAL*8 a,b,h 
      klo=1       !We will find the right place in the table by means of bisection.
      khi=n
!This is optimal if sequential calls to this routine are at random values of x. If sequential calls are in order, and closely spaced, one would do better to store previous values of klo and khi and test if they remain appropriate on the next call

 1    if (khi-klo.gt.1) then 
         k=(khi+klo)/2
         if(xa(k).gt.x)then 
            khi=k
         else 
            klo=k
         endif 
      goto 1
      endif   !klo and khi now bracket the input value of x.
      
      h=xa(khi)-xa(klo) 
      if (h.eq.0.) then
         print *, "bad xa input in splint" !The xa’s must be distinct. 
         STOP
      endif   
      a=(xa(khi)-x)/h    !Cubic spline polynomial is now evaluated. 
      b=(x-xa(klo))/h 
      y=a*ya(klo)+b*ya(khi) + 
     $  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

      return 
      END
