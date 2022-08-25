program structure

use tools
use roots
use integrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This code is designed to calculate the percent differences
! of the peaks of flux, temperature, and luminosity distributions
! between two spherically symmetric geometries. 
! 
! The Schwarzschild geometry is taken to be the reference case. 
! Also provided are two conformally related geometries corresponding to 
! scalar tensor theories of modified gravity. These conformally related 
! geometries can be changed to study other spacetimes in comparison to the 
! Schwarzschild solution. 
! 
! Compile as:
! gfortran tools.f90 roots.f90 integrate.f90 structure.f90
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit double precision (a-h,o-z)
! # of intervals for integration (via simpsons 3/8 rule)
integer, parameter             :: n = 450, iend = 10**3

double precision, allocatable  :: as(:), bs(:,:)
double precision, allocatable  :: flux_dat(:,:), temp_dat(:,:), blum_dat(:,:), rs_dat(:,:)
double precision  :: ts(n+1), ts_s(n+1)


! Physical and mathematical constants (cgs units)
pi    = acos(-1.0d0)           ! Yummy
h     = 6.6261d0*10.0d0**(-27) ! Plancks constant
boltz = 1.3807d0*10.0d0**(-16) ! Boltzmann constant
sb    = 5.67d0*10.0d0**(-5)    ! Stefan Boltzmann constant
c     = 3d0*10d0**10           ! Speed of light
big_g = 6.67430d0*10d0**(-8)   ! Big Daddy G

w        = 10.0d0**9*1.989e33                    ! Supermassive BH mass
acc_rate = 3.154*10.0d0**(-7)*1.989d0*10.0d0**33 ! ~ 1 M_solar/yr
r_h      = 2.0d0*big_g*w/c**2                    ! Schwarz radius


r_min = 3.0d0*r_h   ! radius of innermost stable circular orbit (r_isco)
r_max = 40.0d0*r_h  ! upper bound for radius of accretion disk

freq_min = 10.0d0**(12) ! frequency range for black body spectrum
freq_max = 10.0d0**(16) ! 

phi_scrn = 10d0**(-3)   ! Screening value (approx no scalar tensor
                        ! contribution to gravity)


! To plot (j=1) or not to plot (j .ne. 1)
! flux(r), Temp(r), Luminosity(\nu) distributions in
! the Schwarzschild & scalar tensor cases as they
! are computed
j = 1

! 
! Create data files of Lunminsotiy, Flux, and Temperature distribution peaks
! as functions of parameters (and plot lum, flux, and temp if j=1)
! 

! Schwarzschild peaks for comparison
call flux_s_data(r_min,r_max,j,flx_smax)
call temp_s_data(r_min,r_max,j,temp_smax,rs_smax,ts_s)
call blum_s_data(r_min,r_max,freq_min,freq_max,ts_s,j,blum_smax,freqs_smax)

! Store the maximums of each distribution - Schwarz
aa = flx_smax
bb = rs_smax
dd = temp_smax
ee = blum_smax

! 
! Relative percent differences between peaks
! as functions of the parameters a,b (and r0?)
! used in the logistic scalaron model
! 

! Change these !

!!!!!!!!!!!!!!!!!
alpha = 100d0   ! Coupling constant for f1 model
gamma = 0.001d0 ! Coupling constant for f2 model
ll    = 1       ! ll = 1,2 for f1 or f2 models respectively
mm    = 1       ! mm = 1,2 for logistic or lorentzian model respectively
!!!!!!!!!!!!!!!!!

kend = 50
yend = kend

allocate(as(kend),bs(kend,kend))
allocate(flux_dat(kend,kend),temp_dat(kend,kend),blum_dat(kend,kend),rs_dat(kend,kend))

! Create a grid of 
do k  = 1,kend
! Progress bar
n_count = 100d0*dble(k)/dble(kend)
if (mod(n_count,10) .eq. 0) print*,n_count

do kk = 1,kend ! b loop
zk  = k
zkk = kk

! \phi model parameters
call remap(1.0d0,yend,0.1d0,0.5d0,zk,a)                ! a ~ height ! light field: a << sqrt(\alpha), a << 1/\gamma
if (mm.eq.1) call remap(1.0d0,yend,0.50d0,4.5d0,zkk,b) ! b ~ width
if (mm.eq.2) call remap(1.0d0,yend,1.50d0,3.0d0,zkk,b) ! around r0 for full width half max comparison

r0    = 2.5*r_h    ! Center
r_min = 3.0d0*r_h  ! reset initial integration starting point
r_max = 40.0d0*r_h ! and end point


! Find new r_min to integrate from
outer: do
do i = 1,iend
zi = i
call remap(1.0d0,dble(iend),r_min,r_max,zi,r)

! Check for numerical issues with logistic model
if (dphi(r,r0,a,b) .ne. dphi(r,r0,a,b) .or. phi(r,r0,a,b) .ne. phi(r,r0,a,b) .or. &
ddphi(r,r0,a,b) .ne. ddphi(r,r0,a,b) ) then
print*, 'phi numerical problem'
stop
end if

! Find new r_min
if (f_int(r) .lt. 0d0) then
r_min = r_min+0.00001*r_h
r_max = r_max+0.00001*r_h
cycle outer
endif
end do
exit outer
end do outer

! Find penetration radius defined to be where \phi(rp) ~ phi_scrn
call newt(phi0,dphi0,2.5d0,r_p,0)

! Status updates
print*,'a, b, r_p =',a,b,r_p
print*,'r_min =',r_min/r_h
print*,'r_max =',r_max/r_h
!print*,'Disk size =', (r_max/r_h-r_min/r_h)


! Store selected model parameters and compute flux, temp, and luminosity distributions
as(k)    = a
bs(k,kk) = r_p ! Use penetration radius to parameterize data
call flux_data(ll,a,b,r0,r_min,r_max,j,flx_max)
call temp_data(ll,a,b,r0,r_min,r_max,j,temp_max,rs_max,ts)
call blum_data(ll,a,b,r0,r_min,r_max,freq_min,freq_max,ts,j,blum_max,freqs_max)


! Check for problems
if ( (flx_max.eq.0) .or. (rs_max.eq.0) .or. (temp_max.eq.0) .or. (blum_max.eq.0) ) then
print*,'Error: Computed to zero'
print*,a,b
print*,real(flx_max),real(rs_max),real(temp_max),blum_max
stop
end if

! Store the maximums of each distribution
flux_dat(k,kk) = flx_max
rs_dat(k,kk)   = rs_max
temp_dat(k,kk) = temp_max
blum_dat(k,kk) = blum_max

! Check
print*,''
print*,'flux, temp, lum % changes:', 100d0*(flx_max-aa)/abs(aa), & 
       100d0*(temp_max-dd)/abs(dd),100d0*(blum_max-ee)/abs(ee)
print*,''


end do ! end a
end do ! and b loops

! Write to files - to access via python or whatever you like

! Flux
open(20, file='flux_peaks.dat')
call dump(as,bs,flux_dat,kend,aa,20)
close(20)


! Radius for peak flux
open(21, file='radius_peaks.dat')
call dump(as,bs,rs_dat,kend,bb,21)
close(21)

! Temperature
open(22, file='temperature_peaks.dat')
call dump(as,bs,temp_dat,kend,dd,22)
close(22)

! Luminosity
open(23, file='luminosity_peaks.dat')
call dump(as,bs,blum_dat,kend,ee,23)
close(23)

! Print model choices
print*,'ll',ll, 'mm',mm
if (ll.eq.1) print*,'alpha',alpha
if (ll.eq.2) print*,'gamma',gamma



contains



subroutine dump(xs,ys,zs,n,ref,n_unit)
implicit double precision (a-h,o-z)
double precision :: xs(n),ys(n,n),zs(n,n) ! z(x,y) but y is in inner loop

do k  = 1,n
do kk = 1,n
z = zs(k,kk)
x = xs(k)
y = ys(k,kk)

! Write the percent change between the Schwarzschild case and 
! scalar tensor case
write(n_unit,'(4g24.16)') x,y, 100d0*(z-ref)/abs(ref)

end do
end do


end subroutine



subroutine flux_s_data(r_min,r_max,j,flx_smax)
implicit double precision (a-h,o-z)
double precision, intent(out) :: flx_smax


flx_smax = 0
if (j.eq.1) open(7, file='flux_s.dat')


xend = iend

do i = 1,iend
zi = i
call remap(0.99d0,xend,r_min,r_max,zi,r)

!Grab the maxes
flx_s    = flux_s(r)
flx_smax = max(flx_smax,flx_s)
if (flx_smax .eq. flx_s) rs_smax = r

if (j.eq.1) write(7,'(4g24.16)') r/r_h, flx_s
end do
if (j.eq.1) close(7)

end subroutine



subroutine flux_data(ll,a,b,r0,r_min,r_max,j,flx_max)
implicit double precision (a-h,o-z)
double precision, intent(out) :: flx_max


flx_max = 0
if (j.eq.1) open(1, file='flux.dat')


xend = iend
open(71, file='test')
do i = 1,iend
zi = i
call remap(0.99d0,xend,r_min,r_max,zi,r)
!print*,r/r_h,-c**2*acc_rate/(4.0*pi*r) &
!     * domega_s(r)/(enrgy_s(r)-omega_s(r)*ang_s(r))**2, &
!     -c**2*acc_rate/(4.0*pi*sqrt(-g(r,ll))) &
!     * domega(r,ll)/(enrgy(r,ll)-omega(r,ll)*ang(r,ll))**2

!print*,r/r_h,ang(r,ll),(dgtt(r)*f2(r,gamma)+gtt(r)*df2(r,gamma))/ (2.0/r*gtt(r)-dgtt(r))

!print*,r/r_h,phi(r,r0,a,b),dphi(r,r0,a,b),ddphi(r,r0,a,b)
!print*,r/r_h,(1.0d0+exp(b/r_h*(r-r0))),(1.0d0+exp(b/r_h*(r-r0)))**2,(1.0d0+exp(b/r_h*(r-r0)))**3

!write(71,'(4g24.16)') r/r_h,1d0/(4.0*pi*sqrt(-g_s(r)))* domega_s(r)/(enrgy_s(r)-omega_s(r)*ang_s(r))**2, &
!1d0/(4.0*pi*sqrt(-g(r,ll)))* domega(r,ll)/(enrgy(r,ll)-omega(r,ll)*ang(r,ll))**2

!Grab the maxes
flx        = flux(r)
flx_max    = max(flx_max,flx)

if (flx_max .eq. flx) rs_max = r

if (j.eq.1) write(1,'(4g24.16)') r/r_h, flx
end do
if (j.eq.1) close(1)
close(71)

end subroutine



subroutine temp_s_data(r_min,r_max,j,temp_smax,rs_smax,ts_s)
implicit double precision (a-h,o-z)
double precision, intent(out)   :: temp_smax, rs_smax
double precision, intent(inout) :: ts_s(n+1)


temp_smax = 0
rs_smax   = 0
if (j.eq.1) open(8, file='temperature_s.dat')


xend = iend

do i = 1,iend
zi = i
call remap(0.99d0,xend,r_min,r_max,zi,r)

!Grab the maxes
temp_s    = T_s(r)
temp_smax = max(temp_smax,temp_s)
if (temp_smax .eq. temp_s) rs_smax = r 

if (j.eq.1) write(8,'(4g24.16)') r/r_h, temp_s
end do
if (j.eq.1) close(8)


! Get integration nodes 
hh = (r_max-r_min)/n
ts_s(1) = 1d0
do i = 2, n+1
ts_s(i) = T_s(r_min + hh*i)
end do

end subroutine


subroutine temp_data(ll,a,b,r0,r_min,r_max,j,temp_max,rs_max,ts)
implicit double precision (a-h,o-z)
double precision, intent(out)   :: temp_max, rs_max
double precision, intent(inout) :: ts(n+1)


temp_max = 0
rs_max   = 0
if (j.eq.1) open(3, file='temperature.dat')


xend = iend

do i = 1,iend
zi = i
call remap(0.99d0,xend,r_min,r_max,zi,r)

!Grab the maxes
temp     = T(r)
temp_max = max(temp_max,temp)
if (temp_max .eq. temp) rs_max = r ! If there is a new max temp, record the radius


if (j.eq.1) write(3,'(4g24.16)') r/r_h, temp
end do
if (j.eq.1) close(3)

! Get integration nodes
hh = (r_max-r_min)/n
ts(1) = 1d0
do i = 2, n+1
ts(i) = T(r_min + hh*i)
end do


end subroutine



subroutine blum_s_data(r_min,r_max,freq_min,freq_max,ts_s,j,blum_smax,freqs_smax)
implicit double precision (a-h,o-z)
double precision, intent(in)  :: ts_s(n+1)
double precision, intent(out) :: blum_smax,freqs_smax


blum_smax = 0
if (j.eq.1) open(9, file='luminosity_s.dat')


xend = iend

do i = 1,iend
zi = i
call remap(0.99d0,xend,freq_min,freq_max,zi,freq)


!Grab the maxes
blum_s    = brght_s(freq,ts_s)
blum_smax = max(blum_smax,blum_s)
if (blum_smax .eq. blum_s) freqs_smax = freq

if (j.eq.1) write (9,'(4g24.16)') freq, blum_s
end do
if (j.eq.1) close(9)

end subroutine


subroutine blum_data(ll,a,b,r0,r_min,r_max,freq_min,freq_max,ts,j,blum_max,freqs_max)
implicit double precision (a-h,o-z)
double precision, intent(in)   :: ts(n+1)
double precision, intent(out)  :: blum_max,freqs_max

blum_max = 0
if (j.eq.1) open(5, file='luminosity.dat')


xend = iend

do i = 1,iend

zi = i
call remap(0.99d0,xend,freq_min,freq_max,zi,freq)

!Grab the maxes
blum     = brght(freq,ts)
blum_max = max(blum_max,blum)
if (blum_max .eq. blum) freqs_max = freq


if (j.eq.1) write (5,'(4g24.16)') freq, blum
end do
if (j.eq.1) close(5)


end subroutine


! Metric components and derivatives
function gtt(r)
gtt  = -(1.0d0-2.0d0*big_g*w/(c**2*r))
end function

function dgtt(r)
dgtt = -2.0d0*big_g*w/(c**2*r**2)
end function

function ddgtt(r)
ddgtt = 4.0d0*big_g*w/(c**2*r**3)
end function

function dddgtt(r)
dddgtt = -12.0d0*big_g*w/(c**2*r**4)
end function


! Logistic scalaron model
! (static solutions in https://arxiv.org/pdf/1704.04114.pdf)
function phi(r,r0,a,b)
if (mm.eq.1) phi = a/(1.0d0+exp(b/r_h*(r-r0)))

if (mm.eq.2) then
if (r.lt.0) then
phi = a
else
phi = a/(1d0+(r/(b*r_h))**2)
endif
endif

end function

function dphi(r,r0,a,b)
if (mm.eq.1) dphi = -a*b/r_h*exp(b/r_h*(r-r0))/(1.0d0+exp(b/r_h*(r-r0)))**2

if (mm.eq.2) then
if (r.lt.0) then
dphi = 0
else
dphi = -2d0*a*(b*r_h)**2*r/(r**2+(b*r_h)**2)**2
endif
endif

end function

function ddphi(r,r0,a,b)
if (mm.eq.1) ddphi = a*(b/r_h)**2*(exp(b/r_h*(r-r0))-1.0d0)*exp(b/r_h*(r-r0)) &
        / (exp(b/r_h*(r-r0)) + 1.0d0)**3

if (mm.eq.2) then
if (r.lt.0) then
ddphi = 0
else
ddphi = 2d0*a*(b*r_h)**2*(3d0*r**2-(b*r_h)**2)/(r**2+(b*r_h)**2)**3
endif
endif

end function

function dddphi(r,r0,a,b)
dddphi = -a*(b/r_h)**3*(exp(2d0*b/r_h*(r-r0))-4d0*exp(b/r_h*(r-r0))+1d0) &
        / (exp(b/r_h*(r-r0)) + 1.0d0)**4
end function

! For finding penetration radius
function phi0(r)
if (mm.eq.1) phi0 = a/(1.0d0+exp(b*(r-2.5d0))) - phi_scrn

if (mm.eq.2) then
if (r.lt.0) then
phi0 = a - phi_scrn
else
phi0 = a/(1d0+(r/b)**2) - phi_scrn
endif
endif
end function

function dphi0(r)
if (mm.eq.1) dphi0 = -a*b*exp(b*(r-2.5d0))/(1.0d0+exp(b*(r-2.5d0)))**2

if (mm.eq.2) then
if (r.lt.0) then
dphi0 = 0
else
dphi0 = -2d0*a*(b)**2*r/(r**2+(b)**2)**2
endif
endif

end function


! Inflection point of phi
function ddphi0(r)
if (mm.eq.1) ddphi0 = a*(b)**2*(exp(b*(r-2.5d0))-1.0d0)*exp(b*(r-2.5d0)) &
        / (exp(b*(r-2.5d0)) + 1.0d0)**3

if (mm.eq.2) then
if (r.lt.0) then
ddphi0 = 0
else
ddphi0 = 2d0*a*(b)**2*(3d0*r**2-(b)**2)/(r**2+b**2)**3
endif
endif
end function

function dddphi0(r)
dddphi0 = -a*(b)**3*(exp(2d0*b*(r-2.5d0))-4d0*exp(b*(r-2.5d0))+1d0) &
        / (exp(b*(r-2.5d0)) + 1.0d0)**4
end function


! Conformal Transformation models !
! 
! If you change these, you'll have to find
! new model parameters that don't stop the routine
! due to conditions in line 95
! 


! f_1(\phi) = 1 + \phi^2 / \alpha
function f1(r,alpha)
f1 = 1.0d0 + phi(r,r0,a,b)**2/alpha
end function

function df1(r,alpha)
df1 = 2.0d0*phi(r,r0,a,b)/alpha * dphi(r,r0,a,b)
end function

function ddf1(r,alpha)
ddf1 = 2.0d0/alpha*(dphi(r,r0,a,b))**2 &
       + 2.0d0*phi(r,r0,a,b)/alpha*ddphi(r,r0,a,b)
end function

function dddf1(r,alpha)
dddf1 = 6.0d0/alpha*dphi(r,r0,a,b)*ddphi(r,r0,a,b) &
       + 2.0d0*phi(r,r0,a,b)/alpha*dddphi(r,r0,a,b)
end function


! f_2(\phi) = exp(\gamma*\phi)
function f2(r,gamma)
f2 = exp(gamma*phi(r,r0,a,b))
end function

function df2(r,gamma)
df2 = gamma*f2(r,gamma)*dphi(r,r0,a,b)
end function

function ddf2(r,gamma)
ddf2 = gamma**2*f2(r,gamma)*(dphi(r,r0,a,b))**2 &
       + gamma*f2(r,gamma)*ddphi(r,r0,a,b)
end function

function dddf2(r,gamma)
dddf2 = gamma**3*f2(r,gamma)*(dphi(r,r0,a,b))**3 &
       +3.0d0*gamma**2*f2(r,gamma)*dphi(r,r0,a,b)*ddphi(r,r0,a,b) &
       +gamma*f2(r,gamma)*dddphi(r,r0,a,b)
end function


! Determinant including conformal transform (\theta = \pi/2, \dot{\theta} = 0)
function g(r,ll)
if (ll.eq.1) g = -r**2*f1(r,alpha)**3
if (ll.eq.2) g = -r**2*f2(r,gamma)**3
end function

! Without transformation (\theta = \pi/2, \dot{\theta} = 0)
function g_s(r)
g_s = -r**2
end function



! Energy and angular momentum and derivatives with conformal factor
function enrgy(r,ll) ! ll = 1 -> f1, ll = 2 -> f2, units ergs/gram
if (ll.eq.1) enrgy = c**2*sqrt( gtt(r)**2*(2.0*f1(r,alpha)+df1(r,alpha)*r) &
                     /(r*dgtt(r)-2.0*gtt(r)) )

if (ll.eq.2) enrgy = c**2*sqrt( gtt(r)**2*(2.0*f2(r,gamma)+df2(r,gamma)*r) &
                     /(r*dgtt(r)-2.0*gtt(r)) )
end function

function denrgy(r,ll)
denrgy = omega(r,ll)*dang(r,ll)
end function


function ang(r,ll)
if (ll.eq.1) ang = c*r*sqrt( (dgtt(r)*f1(r,alpha)+gtt(r)*df1(r,alpha)) &
                   / (2.0/r*gtt(r)-dgtt(r)) )

if (ll.eq.2) ang = c*r*sqrt( (dgtt(r)*f2(r,gamma)+gtt(r)*df2(r,gamma)) &
                   / (2.0/r*gtt(r)-dgtt(r)) )
end function

function dang(r,ll)
if (ll.eq.1) dang = ang(r,ll)/(r) + (c*r)**2/(2.0*ang(r,ll)) &
                    *( (ddgtt(r)*f1(r,alpha)+2.0*dgtt(r)*df1(r,alpha) &
                    +gtt(r)*ddf1(r,alpha)) *(2.0/r*gtt(r)-dgtt(r)) &
                    -(dgtt(r)*f1(r,alpha)+gtt(r)*df1(r,alpha)) &
                    *(2.0/r*dgtt(r)-2.0/r**2*gtt(r)-ddgtt(r)) ) &
                    / (2.0/r*gtt(r)-dgtt(r))**2

if (ll.eq.2) dang = ang(r,ll)/(r)+(c*r)**2/(2.0*ang(r,ll)) &
                    *( (ddgtt(r)*f2(r,gamma)+2.0*dgtt(r)*df2(r,gamma) &
                    +gtt(r)*ddf2(r,gamma)) *(2.0/r*gtt(r)-dgtt(r)) &
                    -(dgtt(r)*f2(r,gamma)+gtt(r)*df2(r,gamma)) &
                    *(2.0/r*dgtt(r)-2.0/r**2*gtt(r)-ddgtt(r)) ) &
                    / (2.0/r*gtt(r)-dgtt(r))**2
end function

! angular velocity \omega = d\varphi/dt
function omega(r,ll)
omega = c**2*ang(r,ll)/enrgy(r,ll) * (-gtt(r)/r**2)
end function

function domega(r,ll)
domega = c**2*(2d0*gtt(r)/r-dgtt(r))*ang(r,ll)/(enrgy(r,ll)*r**2)
end function


! Energy and angular momentum schwarzschild
function enrgy_s(r)
enrgy_s = c**2*sqrt( 2.0*gtt(r)**2/(r*dgtt(r)-2.0*gtt(r)) )
end function

function denrgy_s(r)
denrgy_s = omega_s(r)*dang_s(r)
end function


function ang_s(r)
ang_s = c*r*sqrt( dgtt(r)/(2.0/r*gtt(r)-dgtt(r)) )
end function

function dang_s(r)
dang_s = ang_s(r)/r+(c*r)**2/(2.0*ang_s(r))*(ddgtt(r)*(2.0/r*gtt(r)-dgtt(r)) &
       - dgtt(r)*(2.0/r*dgtt(r)-2.0/r**2*gtt(r)-ddgtt(r)) ) &
       / (2.0/r*gtt(r)-dgtt(r))**2
end function



! angular velocity \omega for schwarzschild
function omega_s(r)
omega_s = c**2*ang_s(r)/enrgy_s(r) * (-gtt(r)/r**2)
end function

function domega_s(r)
domega_s = c**2*(2d0*gtt(r)/r-dgtt(r))*ang_s(r)/(enrgy_s(r)*r**2)
end function


! Flux, Temperature, and Luminosity for conformally related geometry
function f_int(r) ! flux integrand
f_int = (enrgy(r,ll)-omega(r,ll)*ang(r,ll))*dang(r,ll)
end function

function flux(r)
pp = -c**2*acc_rate/(4.0*pi*sqrt(-g(r,ll))) &  !sqrt(-g(r,ll))
     * domega(r,ll)/(enrgy(r,ll)-omega(r,ll)*ang(r,ll))**2

! Integrate
call simp(f_int,r_min,r,n,y)

flux = pp*y
end function

! Temperature - Steffan-Boltzmann Law
function T(r)
T = (flux(r)/sb)**(0.25)
end function

! Redshift factor for conformal metric
function gz(r,ll)
if (ll.eq.1) gz = sqrt(-f1(r,alpha)*(c**2*gtt(r)+r**2*omega(r,ll)**2)) &
                  / (1.0-ang(r,ll)*omega(r,ll)/enrgy(r,ll)) / c

if (ll.eq.2) gz = sqrt(-f2(r,gamma)*(c**2*gtt(r)+r**2*omega(r,ll)**2)) &
                  / (1.0-ang(r,ll)*omega(r,ll)/enrgy(r,ll)) / c
end function


function brght_int(r,freq,temp)
freq_e = freq/gz(r,ll)
brght_int = r/(exp(freq_e*h/(boltz*temp))-1.0d0)
end function

function brght(freq,t0) ! Luminsoity (inclination i=0)
double precision :: t0(n+1)
! freq = radiation freq in local rest frame (observed)
! freq_e = emission frequency in local rest frame
qq = 8d0*pi*h*2d0*pi/c**2

! Integrate
call simp_lum(brght_int,freq,t0,r_min,r_max,n,z)

brght =  freq*qq*(freq**3)*z
end function



! Flux, Temperature, and Luminosity for Schwarzschild metric
function f_int_s(r) ! flux integrand
f_int_s = (enrgy_s(r)-omega_s(r)*ang_s(r))*dang_s(r)
end function

function flux_s(r)
pp_s = -c**2*acc_rate/(4.0*pi*sqrt(-g_s(r))) &
     * domega_s(r)/(enrgy_s(r)-omega_s(r)*ang_s(r))**2

call simp(f_int_s,r_min,r,n,y_s)

flux_s = pp_s*y_s
end function

function T_s(r)
T_s = (flux_s(r)/sb)**(.25)
end function

! Redshift factor for Schwarzschild metric
function gz_s(r)
gz_s = sqrt(-(c**2*gtt(r)+r**2*omega_s(r)**2)) &
       / (1d0-ang_s(r)*omega_s(r)/enrgy_s(r)) / c

end function

function brght_int_s(r,freq,temp) ! integrand of luminosity
freq_e = freq/gz_s(r)
brght_int_s = r/(exp(freq_e*h/(boltz*temp))-1.0d0)
end function

function brght_s(freq,t0) ! Luminsoity (inclination i=0)
double precision :: t0(n+1)
! freq = radiation freq in local rest frame (observed)
! freq_e = emission frequency in local rest frame
qq_s  = 8.0d0*pi*h*2.0d0*pi/c**2



call simp_lum(brght_int_s,freq,t0,r_min,r_max,n,z_s)

brght_s =  freq*qq_s*(freq**3)*z_s
end function



end program
