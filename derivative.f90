!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Differentiation schemes 
!     forward difference, central difference, and Pade schemes
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Aug. 27, 2015
!-----------------------------------------------------------------------------!

program differentiation_test
implicit none
integer::j,n,np
real*8 ::h,xp,fp
real*8,allocatable::x(:),u(:),upe(:)
real*8,allocatable::up1(:),up2(:),up3(:),up4(:)

!number of grid points
n = 8

!step size between x:[2,6]
h = 4.0d0/dfloat(n)

allocate(x(0:n),u(0:n),upe(0:n))

do j=0,n
!domain
x(j)   = 2.0d0 + dfloat(j)*h
!function
u(j)   = dsin(x(j))/(x(j)**3)
!exact solution for derivative
upe(j) = (x(j)*dcos(x(j)) - 3.0d0*dsin(x(j)))/(x(j)**4)
end do

!numerical solutions for the derivative 
allocate(up1(0:n),up2(0:n),up3(0:n),up4(0:n))

!up1 = using first-order forward difference
call f1(u,up1,h,n)

!up2 = using second-order central difference
call c2(u,up2,h,n)

!up3 = using fourth-order central difference
call c4(u,up3,h,n)

!up4 = using fourth-order Pade scheme
call pade4(u,up4,h,n)

! Writing data to text file
open(11, file="numerical.plt")
write(11,*)'variables ="x","f1","c2","c4","pade4","exact"'
do j=0,n
write(11,100) x(j),up1(j),up2(j),up3(j),up4(j),upe(j)
end do
close(11)

100 format(6ES18.6)

! Writing exact solution using np points
np = 2000 !use 2000 points to plot curve
open(12, file="exact_curve.plt")
write(12,*)'variables ="x","upe"'
	do j=0,np
		xp = 2.0d0 + dfloat(j)*4.0d0/(dfloat(np))
		fp = (xp*dcos(xp) - 3.0d0*dsin(xp))/(xp**4)
		write(12,*) xp,fp
	end do
close(12)

     
end


!-----------------------------------------------------------------------------!
!f1: 	1st-order forward difference scheme for the first order derivative(up)	
!-----------------------------------------------------------------------------!
subroutine f1(u,up,h,n)
implicit none
integer :: n,i
real*8  :: h
real*8, dimension (0:n)	:: u,up

do i=0,n-1
up(i) = (u(i+1)-u(i))/h
end do

!sided difference for i=n
i=n
up(i) = (u(i)-u(i-1))/h

return
end


!-----------------------------------------------------------------------------!
!c2: 	2nd-order central difference scheme for the first order derivative(up)	
!-----------------------------------------------------------------------------!
subroutine c2(u,up,h,n)
implicit none
integer :: n,i
real*8  :: h
real*8, dimension (0:n)	:: u,up

do i=1,n-1
up(i) = (u(i+1)-u(i-1))/(2.0d0*h)
end do

!sided difference for i=0
i=0
!up(i) = (u(i+1)-u(i))/h
up(i) = (-3.0d0*u(i) + 4.0d0*u(i+1) - u(i+2))/(2.0d0*h)

!sided difference for i=n
i=n
!up(i) = (u(i)-u(i-1))/h
up(i) = (-3.0d0*u(i) + 4.0d0*u(i-1) - u(i-2))/(-2.0d0*h)

return
end

!-----------------------------------------------------------------------------!
!c4: 	4th-order central difference scheme for the first order derivative(up)	
!-----------------------------------------------------------------------------!
subroutine c4(u,up,h,n)
implicit none
integer :: n,i
real*8  :: h
real*8, dimension (0:n)	:: u,up

do i=2,n-2
up(i) = (-u(i+2)+8.0d0*u(i+1)-8.0d0*u(i-1)+u(i-2))/(12.0d0*h)
end do

!sided difference for i=0 (2nd-order)
i=0
up(i) = (-3.0d0*u(i) + 4.0d0*u(i+1) - u(i+2))/(2.0d0*h)

!sided difference for i=1 (3rd-order)
i=1
up(i) = (-2.0d0*u(i-1) - 3.0d0*u(i) + 6.0d0*u(i+1) - u(i+2))/(6.0d0*h)

!sided difference for i=n (2nd-order)
i=n
up(i) = (-3.0d0*u(i) + 4.0d0*u(i-1) - u(i-2))/(-2.0d0*h)

!sided difference for i=n-1 (3rd-order)
i=n-1
up(i) = (-2.0d0*u(i+1) - 3.0d0*u(i) + 6.0d0*u(i-1) - u(i-2))/(-6.0d0*h)


return
end




!-----------------------------------------------------------------------------!
!pade4: 4th-order compact Pade scheme for the first order derivative(up)
!		3-4-3
!       3rd-order b.c.
!       4th-order interior		
!-----------------------------------------------------------------------------!
subroutine pade4(u,up,h,n)
implicit none
integer :: n,i
real*8  :: h
real*8, dimension (0:n)	:: u,a,b,c,r,up

i=0
b(i) = 1.0d0
c(i) = 2.0d0
r(i) = (-5.0d0*u(i) + 4.0d0*u(i+1) + u(i+2))/(2.0d0*h)

do i=1,n-1
a(i) = 1.0d0
b(i) = 4.0d0
c(i) = 1.0d0
r(i) = (3.0d0/h)*(u(i+1)-u(i-1))
end do

i=n
a(i) = 2.0d0
b(i) = 1.0d0
r(i) = (-5.0d0*u(i) + 4.0d0*u(i-1) + u(i-2))/(-2.0d0*h)

call tdma(a,b,c,r,up,0,n)

return
end


!------------------------------------------------------------------!
!Tridiagonal matrix algorithm (TDMA)
!Thomas algorithm
!solution tridiagonal systems
!a: lower diagonal
!b: main diagonal
!c: upper diagonal
!r: source vector
!x: solution vector
!   for indices s(start) to e(end)
!   i: s,s+1,s+2, ....,i,....,e 
!
!Note: a(s) and c(e) are dummy coefficients, not used.
!------------------------------------------------------------------!

subroutine tdma(a,b,c,r,x,s,e)
implicit none
integer s,e,i
real*8, dimension(s:e) ::a,b,c,r,x    

! forward elimination phase
do i=s+1,e
b(i) = b(i) - a(i)/b(i-1)*c(i-1)
r(i) = r(i) - a(i)/b(i-1)*r(i-1)
end do
! backward substitution phase 
x(e) = r(e)/b(e)
do i=e-1,s,-1
x(i) = (r(i)-c(i)*x(i+1))/b(i)
end do

return
end
