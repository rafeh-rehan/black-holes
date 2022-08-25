program shadows


use tools

!!!!!!!!!
! 
! 	This code is designed to draw shadows of spherically symmetric 
! black holes. Given a metric and its derivatives, the geodesic equation 
! can be solved symplectically in order to trace the path of light from 
! a source location to an observer. 
! 
! 	The Schwarzschild metric is used as the base case for spherically 
! symmetric black holes and their asymptotic solutions. A second 
! spherically symmetric metric found in https://arxiv.org/pdf/0806.0679.pdf
! is also included. This second metric can be changed to study other 
! spacetimes in comparison with the Schwarzschild geometry.
! 
! An implicit Runge-Kutta method for solving sets of differential equations
! is used to solve the geodesic equations (reference: https://arxiv.org/abs/1704.04114)
! 
! Compile as:
! gfortran -O3 -fdefault-real-8 tools.f90 shadows.f90
! 
!!!!!!!!! 

implicit real (a-h,o-z)
real  :: y(8)  ! State vector


pi = acos(-1.0)
dt = 0.0004    ! timestep for Gauss-Legendre solver


! Black hole mass
bh_m  = 1.0 

! Model parameters for metric in https://arxiv.org/pdf/0806.0679.pdf
a   = -3e-3    ! a in [-3e-3, 3e-3]
n   = 1        ! n = 1,2,3,4
eta = n*bh_m/a ! eta = n*M/a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose spacetime for black hole shadow     ! 
! 0 for Schwarzschild, 1 for other spacetime ! 
nn=0



! Shoot photons towards a black hole with increasing impact parameters
jjend = 200     ! Number of iterations of the impact parameter
xjend  = jjend

do jj = 1, jjend
xi = jj

! Increase impact param until photon doesnt plunge 
call remap(1.0,xjend,0.92,1.01,xi,u)
! May need to change  ^    ^  range of impact parameters considered

! Initial photon coordinates (circle parameterization)
ecks = -100.0
why  = u*3.0*sqrt(3.0) !!! Impact parameter
zed  = 0.0
call cart2sph(ecks, why, zed, are, thta, p)

r_0   = are     ! Observer's radial position
thta0 = pi/2.0  ! Observer's inclination


!   Pack initial conditions into state vector
    y(1) = 0.0       ! t
    y(3) = are       ! r - units of GM = c = 1
    y(5) = p         ! \phi
    y(7) = thta      ! \theta

!   Null geodesic conditions. Change ll=1 case if using another metric
!   dt/ds for ll=0 schwarz
    if (ll.eq.0) y(2) = sqrt(sin(y(7))**2*cos(y(5))**2*(1-2.0/y(3))**(-2) &
           + cos(y(7))**2*cos(y(5))**2*(1-2.0/y(3))**(-1) &
           + sin(y(5))**2*(1-2.0/y(3))**(-1) )                        
    
!   dt/ds for ll=1 metric in https://arxiv.org/pdf/0806.0679.pdf
    if (ll.eq.1) y(2) = sqrt( ( grr(y(3))*(sin(y(7))*cos(y(5)))**2 &
           + y(3)**2*( (cos(y(7))*cos(y(5))/y(3))**2 & 
           + (sin(y(5))/(y(3)))**2 ) )/(-gtt(y(3))) ) 

    
!   Initial velocity in the +x direction only
    y(4) = sin(y(7))*cos(y(5))         ! dr/ds
    y(6) = -sin(y(5))/(y(3)*sin(y(7))) ! d\phi/ds
    y(8) = cos(y(7))*cos(y(5))/y(3)    ! d\theta/ds

    E0 = energy(y)      
    A0 = am(y)


! Evolution loop for solving geodesic equations
jk = 0
do l = 0, 640000 ! Number of steps for solver

!       Write x,y,z in spherical coords
	r1 = y(3)*sin(y(7))*cos(y(5))
	r2 = y(3)*sin(y(7))*sin(y(5))
	r3 = y(3)*cos(y(7))


!       Sky Coordinates for an observer's sky 
!       - https://arxiv.org/abs/gr-qc/0308023
	rr = 0.0
	! if light hits the observer at r = r0
	if ( (abs(y(3)-r_0) .le. 0.05) .and. (r1 .gt. 3.0)) then
	a = -y(3)**2 * sin(thta0) * y(6)/y(4)
	b =  y(3)**2 * y(8)/y(4)
	
	
!       Sketch circle with radius a**2 + b**2, due to radial symmetry 
	rr = sqrt(a**2+b**2)
	do kk = 1,100
	xk = kk
	call remap(1.0,100.0,0.0,2.0*pi,xk,ss)
	write (*,'(4g24.16)') rr*cos(ss), rr*sin(ss)
	end do

	exit
	endif 
 
! Which equations to use: nn=0 for Schwarz, nn=1 for other 
if (nn.eq.1) then
! Stitch solver to Schwarzschild metric - avoids undefined terms in other metric
if ( y(3) .lt. abs(eta) ) then ! radial conditions
ll = 1 ! Other metric
else 
ll = 0 ! Schwarzschild metric
end if
end if

if (nn.eq.0) ll=0 
call gl10(y, dt, ll)  
end do

if (rr .ne. 0.0) exit
end do




contains

! Energy of particle
function energy(vec)
real energy, vec(8)

if (ll.eq.0) energy = (1-2.0/vec(3)) * vec(2)
if (ll.eq.1) energy = -gtt(vec(3)) * vec(2)
end function energy

! Angular momentum of particle
function am(vec)
real am, vec(8)

am = (vec(3)*sin(vec(7)))**2 *vec(6)
end function am


!!! Metric equations and derivatives. Change these to study other
!!! spherically symmetric spactimes


function gtt(r)
gtt = -(1.0+3.0*bh_m/eta-2.0*bh_m/r-(1.0+6.0*bh_m/eta)*(r/eta) & 
        + (r/eta)**2*(1.0+(1.0+6.0*bh_m/eta)*log(a*(1.0+eta/r))) )
end function

function drgtt(r)
drgtt = - 2.0*bh_m/r**2 + (1.0+6.0*bh_m/eta)*(1.0+1.0/(1.0+eta/r))/eta & 
        - 2.0*(r/eta**2)*(1.0+(1.0+6.0*bh_m/eta)*log(a*(1.0+eta/r)))
end function


function grr(r)
grr = -1.0/gtt(r)
end function

function drgrr(r)
drgrr = (1.0/gtt(r))**2*drgtt(r)
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Geodesic equations of motion for metric given in https://arxiv.org/pdf/0806.0679.pdf
! ds^2 = -exp(\nu(r))dt^2 + exp(\lambda(r))dr^2 + r^2(d\theta^2 + sin\theta d\phi^2)

subroutine evalf(y, dydx)
real        :: y(8), dydx(8) 

!!!     Derivatives of state vector elements from geodesic 
!!!     equation. Change as desired with metric

!       dt/ds & d^2t/ds^2
        dydx(1) = y(2)
        dydx(2) = -1.0/gtt(y(3))*drgtt(y(3))*y(4)*y(2)
        
!       dr/ds & d^2r/ds^2
        dydx(3) = y(4)
        dydx(4) = -(-0.5/grr(y(3))*drgtt(y(3))*(y(2))**2 & !time
                  +0.5/grr(y(3))*drgrr(y(3))*y(4)**2 &     !radial
                  -1.0/grr(y(3))*y(3)*y(8)**2 &            !theta
                  -1.0/grr(y(3))*y(3)*sin(y(7))**2 *(A0/(y(3)**2))**2) !phi
        
        
!       d\phi/ds & d^2phi/ds^2
        dydx(5) = y(6)
        dydx(6) = - ( 2.0/y(3) * y(4)*y(6) &
        + 2.0*cos(y(7))/sin(y(7)) * y(8)*y(6) )

        
!       d\theta/ds & d^2\theta/ds^2
        dydx(7) = y(8)
        dydx(8) = - ( 2.0/y(3) * y(4)*y(8) & 
        - sin(y(7))*cos(y(7)) * y(6)**2 )


end subroutine evalf


!!! State vector derivatives for Schwarzschild metric

subroutine evalf2(y, dydx)
        real y(8), dydx(8)
        
        
!       dt/ds & d^2t/ds^2
        dydx(1) = y(2)
        dydx(2) = - 2.0/(y(3)*(y(3)-2.0)) * y(2) * y(4)
        

!       dr/ds & d^2r/ds^2
        dydx(3) = y(4)
        dydx(4) = - ( 1.0/(y(3)*(2.0-y(3))) * y(4)**2 &
        + 1.0*(y(3)-2.0)/(y(3)**3) * y(2)**2 &
        + (2.0-y(3))*(y(8)**2+(sin(y(7)))**2*y(6)**2) )
        
!       d\phi/ds & d^2phi/ds^2
        dydx(5) = y(6)
        dydx(6) = - ( 2.0/y(3) * y(4)*y(6) &
        + 2.0*cos(y(7))/sin(y(7)) * y(8)*y(6) )
     
!       d\theta/ds & d^2\theta/ds^2
        dydx(7) = y(8)
        dydx(8) = - ( 2.0/y(3)*y(4)*y(8) & 
        - sin(y(7))*cos(y(7))*y(6)**2 ) 

end subroutine


! Implicit Runge-Kutta method for solving sets of differential equations.
! From https://arxiv.org/abs/1704.04114.

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt, ll)
integer, parameter :: s = 5, n = 8	      
real y(n), g(n,s), dt; integer i, k, ll


! Butcher table for 10th order Gauss-Legendre method
real, parameter :: a(s,s) = reshape([ &
0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
0.5923172126404727187856601017997934066Q-1], [s,s])
real, parameter ::   b(s) = [ &
1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
1.1846344252809454375713202035995868132Q-1]

! Iterate trial steps
g = 0.0; do k = 1,16
g = matmul(g,a)
do i = 1,s

if (ll .eq. 0) then
call evalf2(y + g(:,i)*dt, g(:,i)) ! Schwarzschild
else if (ll .eq. 1) then
call evalf(y + g(:,i)*dt, g(:,i))  ! Other metric
end if 

end do
end do

! Update the solution
y = y + matmul(g,b)*dt
end subroutine gl10


end program
