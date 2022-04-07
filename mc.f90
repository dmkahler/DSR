module sorb_react
contains
function dsdt(con,sor,k_d,gamma)
integer p
parameter(p=kind(1.0D0))
real (p) :: con, sor, dsdt
! k_d is the "distribution coefficient" and is the equilibrium constant defining sorption
! gamma is the rate parameter, and also controls the reversability
real (p) :: k_d, gamma
! reversable transient sorption
dsdt = gamma * (k_d * con - sor)
end function

function drdt(con,a,b,dt)
integer p
parameter(p=kind(1.0D0))
real (p) :: con, dt, drdt
real (p) :: a, b
! 1625.6 kg m^-3 saturation of Ag- in water.
if (con < react) then
    ! ORIGINAL: 
    ! drdt = (react - con) / dt
    ! this results in a range from "react", which was 0.15 minus C, which ranged from 0 to 1e-5 
    ! all over dt = 0.001; therefore, drdt ranged from 150 to 149.95
    ! NEW:
    ! We seek a function that begins with a value of 150, but decreases as the concentration approaches 15, but
    ! continues to increase (not convergent).
    drdt = (1625.6 - con)**a)/b
    ! at first, con=0, drdt should equal 150
else
    drdt = 0
endif
end function
end module

program madicubes

use sorb_react

! This program solves the diffusion equation with a second-order finite
! difference model in an axi-symetric geometry for a bucket with the original
! silver-empregnated disk in the center in the bottom of the bucket.

! Full model in cylindrical coordinates
! Axi-symetric around z = 0
! Assuming radial symmetry d/dtheta=0

! REVISION HISTORY
! 07 Jan 2016   Adapted from madipieces.f90, hardcoded instantanious mixing (via volumetric averaging)
! 09 Feb 2016   Finished corrections on diffusion routines
! 13 Feb 2016   Finished diffusion routines and mixing outside the ceramic
! 14 Feb 2016   Problem with sorbtion only increasing total Ag
! 18 Feb 2016   Fix the dsdt term -- new rev. lin. kinetic sorb. function was dt-independent
! 07 Sep 2016   Trying averaging routine again.
! 10 Sep 2016   Rewrote most of program, found lots of bugs, fixed volume shells and averaging routine
! 18 Oct 2016   Fixed the use of c_ave, in the averaging - now using c1(rm,zm) for most.
! 02 Nov 2016   Converted to large MadiDrop geometry to test v0036 and v0037.
! 02 Nov 2016   Changing the mass transfer function

! RUN HISTORY
! set   date        dt      dr=dz   k_1     k_2     react   runtime     notes
! v0001 10Sep2016   0.0010  0.0001  0.2     0.18    0       17 min      checking run
! v0002 11Sep2016   0.0010  0.0001  0.2     0.18    0       13 min      Ci=100ppb
! set   date        dt      dr=dz   k_d     gamma   react   runtime     notes
! v0003 27Sep2016   0.0010  0.0001  2.5e-4  6.8e-8  0       14          Ci=1ppm, 1 hr, C=0.00006ppm TOO FAST
! v0004 27Sep2016   0.0010  0.0001  2.5e-4  6.8e-8  0       14          Ci=10ppm, 1 hr, C=0.0006ppm
! v0005 29Sep2016   0.0010  0.0001  2.5e-4  6.8e-10 0       29          Ci=10ppm, 2 hr, C=0.0006ppm
! v0006 29Sep2016   0.0010  0.0001  2.5e-1  6.8e-10 0       29          Ci=10ppm, 2 hr, C=0.0006ppm
! v0007 30Sep2016   0.0010  0.0001  2.5e-1  6.8e-10 0       34          Ci=10ppm, 2 hr, fixed reassignment and fluid computations, reached good equilibrium, but too fast
! v0008 30Sep2016   0.0010  0.0001  2.5e-1  1e-10   0       33          Ci=100ppm, 2 hr, uncertain
! v0009 30Sep2016   0.0010  0.0001  2.5e-1  1e-8    0       tt          Ci=100ppm, 2 hr, minute resolution
! v0010 30Sep2016   0.0010  0.0001  2.5e-1  1e-8    0       tt          Ci=100ppm only water
! v0011 30Sep2016   0.0010  0.0001  2.5e-4  1e-9    0       tt          Ci=100ppm only water
! v0012 30Sep2016   0.0010  0.0001  2.5e-4  1e-9    0       tt          Ci=1000ppm-adjusted for volume
! v0013 30Sep2016   0.0010  0.0001  1e-4    1e-9    0       32 min      Ci=1000ppm
! v0014 14Oct2016   0.0010  0.0001  1e-4    1e-9    0       tt min      Ci=1000ppm, just one second
! v0015 14Oct2016   0.0010  0.0001  0       0       0       stopped     Ci=1000ppm, just one hour
! v0016 14Oct2016   0.0010  0.0001  0       0       0       tt min      Ci=1000ppm, everywhere
! v0017 14Oct2016   0.0010  0.0001  0       0       0       tt min      Ci=1000ppm, everywhere, diagnostic
! v0018 18Oct2016   0.0010  0.0001  1e-4    1e-9    0       35 min      Ci=1000ppm, post averaging fix
! v0019 18Oct2016   0.0010  0.0001  1e+3    1e+2    0       tt min      Ci=1000ppm, 48 hours, NaN...
! v0020 20Oct2016   0.0010  0.0001  1       1       0       tt min      Ci=1000ppm, 8 hours, NaN...
! set   date        dt      dr=dz   k_d     gamma   react   initial C   notes
! v0021 24Oct2016   0.0010  0.0001  1e-1    1e-1    0       Ci=100ppm, 2 hours, normal run
! v0022 24Oct2016   0.0010  0.0001  5e-1    1e-1    0       Ci=100ppm, 2 hours, normal run
! v0023 24Oct2016   0.0010  0.0001  5e-1    5e-1    0       Ci=100ppm, 2 hours, normal run
! v0024 24Oct2016   0.0010  0.0001  5e-1    7e-1    0       Ci=100ppm, 24 hours, good result, about 10x too fast
! v0025 25Oct2016   0.0010  0.0001  5e-1    1e-1    0       Ci=100ppm, 24 hours, good result, still too fast
! v0026 25Oct2016   0.0010  0.0001  5e-1    1e-3    0       Ci=100ppm, 24 hours, stable but still too fast
! v0027 25Oct2016   0.0010  0.0001  5e-1    1e-5    0       Ci=100ppm, still too fast.  altering both kd&g
! v0028 25Oct2016   0.0010  0.0001  1e-1    1e-7    0       Ci=100ppm, way too slow
! v0029 25Oct2016   0.0010  0.0001  5e-1    5e-1    0       Ci=100ppm, upon further review, i was closer on v0024
! v0030 27Oct2016   0.0010  0.0001  5e-1    5e-4    0       Ci=100ppm, better, but still could be lower gamma
! v0031 27Oct2016   0.0010  0.0001  5e-1    5e-4    1       Ci=0, testing react=1, achieved over 10ppb in 3hrs.
! v0032 31Oct2016   0.0010  0.0001  5e-1    5e-4    1e-1    Ci=0,
! v0033 31Oct2016   0.0010  0.0001  5e-1    5e-4    1e-2    Ci=0, tau=0.4 (0.4 was used for all above)
! set   date        dt      dr=dz   k_d     gamma   react   tau     notes
! v0034 01Nov2016   0.0010  0.0001  5e-1    5e-4    1.5e-1  0.6     new tau from another diff exp.
! v0035 01Nov2016   0.0010  0.0001  7e-1    1e-4    1.5e-1  0.6     quick trial
! v0036 02Nov2016   0.0010  0.0001  7e-1    1e-4    1.5e-1  0.6     cubes, to compare       mc.v0036.f90
! v0037 02Nov2016   0.0010  0.0001  7e-1    1e-4    1.5e-1  0.6     MadiDrop, to compare    md.v0037.f90
! set   date        dt      dr=dz   k_d     gamma   a       b       tau     notes
! v0038 02Nov2016   0.0010  0.0001  7e-1    1e-4                    0.6     checking new mass trans func, no .pdt
! v0039 02Nov2016   0.0010  0.0001  7e-1    1e-4                    0.6

implicit none
integer p
parameter(p=kind(1.0D0))
real (p) :: r_bucket, r_md, h_bucket, h_md, v_bucket, v_md
real (p) :: dr, dz, dt
real (p) ::b_d, n, tau, k_d, gamma, a, b
real (p) :: diff, diff_md
real start, finish, elapsed
integer tm, hr, rm, zm, rmd, zmd
integer r, z, t, h, mi
real (p) :: ci, si
real (p) :: c1(57,101)
real (p) :: c2(57,101)
real (p) :: s1(57,101)
real (p) :: s2(57,101)
real (p) :: m, v, v_water, c_ave
real (p) :: vb(58), vu(58)
real (p) :: term1, term2, term3, term4, term5

call cpu_time(start)

!open(20, file='madi.v0036.pdt')
open(25, file='mc.v0036.adt')
open(40, file='mc.v0036.par')

!****************************************************************************
!*                                                                          *
!*                               GEOMETRY AND TIME                          *
!*                                                                          *
!****************************************************************************

dr = 0.0001                 ! m, dr = 0.1 mm yields 56 nodes in the MadiDrop plus two edge nodes
dz = 0.0001                 ! m, dz = 0.1 mm yields 100 nodes in the MadiDrop plus two edge nodes
rm = 57                     ! number of nodes in the radial direction in the domain
zm = 101                    ! number of nodes in the z-direction in the domain
rmd = 56                    ! number of nodes in the radial direction within MadiDrop
zmd = 100                   ! number of nodes in the z-direction within MadiDrop
r_md = rmd*dr               ! m, radius of the MadiDrop
h_md = zmd*dz               ! m, height of MadiDrop
! The "bucket" size is set to model 1/100th of a 10 L bucket as a tall cylinder.  This model assumes
! 100 adjacent and identical cylinders that will make up the 10 L.
r_bucket = 0.0150           ! m, modified radius of bucket to reflect many MadiDrop pieces (summer
                            ! 2015 form factor).  Radius = 1 cm as a perfectly mirrored column in 
                            ! bucket.
h_bucket = 0.1415           ! m, height of bucket or reactor
dt = 0.0010                 ! s
tm = 60000                  ! number of timesteps in one minute
hr = 24                     ! hours of simulation

! Calculate volume shells
! Lower elements
do r = 1,rm
    vb(r) = 3.14159 * ((r * dr)**2 - ((r-1) * dr)**2) * dz
enddo
vb(rm+1) = 3.14159 * (r_bucket**2 - ((r_md+dr)**2)) * dz
! Upper elements
do r = 1,rm
    vu(r) = 3.14159 * ((r * dr)**2 - ((r-1) * dr)**2) * (h_bucket - h_md + dz)
enddo
vu(rm+1) = 3.14159 * (r_bucket**2 - ((r_md+dr)**2)) * (h_bucket - h_md + dz)
! Other volume elements
v_bucket = (zm)*sum(vb) + sum(vu)                           ! changed from zmd+1 to zm 11 Sep 2016
v_md = zmd * (sum(vb) - vb(rm+1) - vb(rm))
v_water = zmd * (vb(rm+1) + vb(rm)) + sum(vb) + sum(vu)     ! changed from difference to summation

!****************************************************************************
!*                                                                          *
!*                               PARAMETERS                                 *
!*                                                                          *
!****************************************************************************

k_d = 5e-1                  ! m^3 kg^-1; 1 L kg^-1 = 10^-3 m^3 kg^-1 equilibrium sorption parameter
gamma = 5e-4                ! s^-1 rate parameter
a = 1
b = 1
!react = 1.5e-1                ! kg m^-3 s^-1
tau = 0.6 ! empirical value from tritium experiments.  Speculate that a
! higher tau is related to the pore structures in the ceramic having
! greater connectivity than a uniformly distributed homogeneous pore
! structure.

b_d = 0.0013 / (3.14159 * r_md**2 * h_md)   ! kg m^-3, bulk density (mass density of dry ceramic
n = 0.4                                     ! dimensionless, porosity
diff = 1.5d-9                               ! m^2 s^-1, value from literature:
!                                           Johnston and Spiro, 1967.
diff_md = tau * diff                        ! m^2 s^-1, adjusted for ceramic.
write(*,*) 'Geometry and physical parameters loaded'


write(40,*) 'r_bucket =                 ', r_bucket
write(40,*) 'r_md =                     ', r_md
write(40,*) 'h_bucket =                 ', h_bucket
write(40,*) 'h_md =                     ', h_md
write(40,*) 'dr =                       ', dr
write(40,*) 'dz =                       ', dz
write(40,*) 'dt =                       ', dt
write(40,*) 'rmd =                      ', rmd
write(40,*) 'zmd =                      ', zmd
write(40,*) 'tm =                       ', tm
write(40,*) 'hr =                       ', hr
write(40,*) 'tau =                      ', tau
write(40,*) 'diffusion (MadiDrop) =     ', diff_md
write(40,*) 'B_d =                      ', b_d
write(40,*) 'k_d =                      ', k_d
write(40,*) 'gamma =                    ', gamma
write(40,*) 'reaction rate =            ', react
write(40,*) 'Lower volume elements'
do r = 1,rm
    write(40,*) 'Element                ', r, vb(r)
enddo
write(40,*) 'Upper volume elements'
do r = 1,rm
    write(40,*) 'Element                ', r, vu(r)
enddo
write(40,*) 'v_bucket                   ', v_bucket
write(40,*) 'v_MD                       ', v_md
write(40,*) 'v_water                    ', v_water
write(40,*) 'Output files'
write(40,*) 'madi.vXXXX.pdt     concentration matrix'
write(40,*) 'headers:              hour, radius, height, concentration, sorption'
write(40,*) 'madi.vXXXX.adt     concentration average time series'
write(40,*) 'headers:              hour, average concentration'
write(40,*) 'madi.vXXXX.par     parameter list'
write(40,*) 'madi.vXXXX.chk     simulation output verification - NOT ALWAYS USED'

!****************************************************************************
!*                                                                          *
!*                          MODEL INITIALIZATION                            *
!*                                                                          *
!****************************************************************************

! The nodes within the MadiDrop will compute the diffusion directly.  The first node outside the
! MadiDrop will compute the diffusion out of the surface of the MadiDrop.  The second node outside the
! MadiDrop will track the concentration of the bulk water; while subsequent concentrations will not be
! computed, the value of the second node outside will be assigned to the remaining nodes not specified
! outside the MadiDrop.

ci = 0 !1d-4 ! kg m^-3, 10 mg/L (ppm) = 10^-2 kg m^-3, 10 ppb = 10^-5 kg m^-3
si = 0 !
! c preallocation of concentration array
! s preallocation of sorbed mass array - MadiDrop only
! r preallocation of reaction array - MadiDrop only
! NOTE: C(r,z) where c1 and c2 are progressive timesteps

! Assign only outer water to ci
!do r = 1,rm
!    c1(r,zm)=ci
!enddo
!do z = 1,zm
!    c1(rm,z)=ci
!enddo

! Assign all nodes to ci
do r = 1,rm
	do z = 1,zm
        c1(r,z) = ci
		s1(r,z) = si
	enddo
enddo
write(*,*) 'Model specifications loaded and arrays initialized'
write(*,*) 'Iniitial concentration ', ci

! Preparation: 12 FLOPS

!****************************************************************************
!*                                                                          *
!*                             MODEL EVALUATION                             *
!*                                                                          *
!****************************************************************************

write(*,*) 'Starting simulation'
do h = 1,hr
    do mi = 1,60
    do t = 1,tm
        ! INSIDE THE MadiDrop
        ! boundary point !	BOTTOM CENTER OF BUCKET
        ! no-flux boundary condition, phantom points created at r-1,z and r,z-1
        term1 = (c1(2,1)-c1(1,1))/(dr**2) ! note that the denomenator was simplified
        term2 = ((c1(2,1)-(2*c1(1,1))+c1(1,1))/(dr**2))
        term3 = ((c1(1,2)-(2*c1(1,1))+c1(1,1))/(dz**2))
        term4 = dsdt(c1(1,1),s1(1,1),k_d,gamma)                             ! SORPTION TERM
        term5 = drdt(c1(1,1),react,dt) * dr * dz                          ! REACTION TERM
        c2(1,1) = c1(1,1) + (dt * diff_md * (term1 + term2 + term3 )) - (dt * term4 * b_d / n) + (dt * term5)
        s2(1,1) = s1(1,1) + (dt * term4)
        ! boundary edge!	AXIS OF SYMMETRY
        do z = 2,zmd
            ! symmetry boundary condition, phantom points created such that C(r-1) = C(1)
            term1 = (c1(2,z)-c1(1,z))/(dr**2) ! note that the denomenator was simplified
            term2 = ((c1(2,z)-(2*c1(1,z))+c1(1,z))/(dr**2))
            term3 = ((c1(1,z+1)-(2*c1(1,z))+c1(1,z-1))/(dz**2))
            term4 = dsdt(c1(1,z),s1(1,z),k_d,gamma)                         ! SORPTION TERM
            term5 = drdt(c1(1,z),react,dt) * dr * dz                      ! REACTION TERM
            c2(1,z) = c1(1,z) + (dt * diff_md * ( term1 + term2 + term3 )) - (dt * term4 * b_d / n) + (dt * term5)
            s2(1,z) = s1(1,z) + (dt * term4)
        enddo
        ! boundary edge!	BOTTOM OF BUCKET
        do r = 2,rmd
            ! no-flux boundary condition, phantom point created such that C(z-1) = C(z+1)
            term1 = (c1(r+1,1)-c1(r-1,1))/(2*dr*(r-0.5)*dr)
            term2 = (c1(r+1,1)-(2*c1(r,1))+c1(r-1,1))/(dr**2)
            term3 = ((c1(r,2)-(2*c1(r,1))+c1(r,1))/(dz**2))
            term4 = dsdt(c1(r,1),s1(r,1),k_d,gamma)                         ! SORPTION TERM
            term5 = drdt(c1(r,1),react,dt) * dr * dz                      ! REACTION TERM
            c2(r,1) = c1(r,1) + (dt * diff_md * ( term1 + term2 + term3 )) - (dt * term4 * b_d / n) + (dt * term5)
            s2(r,1) = s1(r,1) + (dt * term4)
        enddo
        ! CERAMIC INTERIOR —- NON-BOUNDARY
        do r = 2,rmd
            do z = 2,zmd
                ! Product rule of the first term gives us term1 and term2
                term1 = (c1(r+1,z)-c1(r-1,z))/(2*dr*(r-0.5)*dr)         ! 1/r dc/dr
                term2 = (c1(r+1,z)-(2*c1(r,z))+c1(r-1,z))/(dr**2)       ! d2c/dr2
                ! there is no theta dependance, all terms with d/dtheta go to zero
                term3 = ((c1(r,z+1)-(2*c1(r,z))+c1(r,z-1))/(dz**2))     ! d2c/dz2
                term4 = dsdt(c1(r,z),s1(r,z),k_d,gamma)                     ! SORPTION TERM
                term5 = drdt(c1(r,z),react,dt) * dr * dz                  ! REACTION TERM
                ! Error corrected 09 Jan 2016, term1 was missing from the c2= line.
                c2(r,z) = c1(r,z) + (dt * diff_md * ( term1 + term2 + term3 )) - (dt * term4 * b_d / n) + (dt * term5)
                s2(r,z) = s1(r,z) + (dt * term4)
            enddo
        enddo

        ! OUTSIDE THE MadiDrop
        ! boundary edge!	BOTTOM OF BUCKET
        ! no-flux boundary condition, phantom point created such that C(z-1) = C(z)
        term1 = (c1(rm,1)-c1(rmd,1))/(2*dr*(rmd+0.5)*dr)
        term2 = ((c1(rm,1)-(2*c1(rm,1))+c1(rmd,1))/(dr**2))
        term3 = ((c1(rm,2)-(2*c1(rm,1))+c1(rm,1))/(dz**2))     ! d2c/dz2
        c2(rm,1) = c1(rm,1) + (dt * diff * ( term1 + term2 +term3 ))
        ! FLUID —- LOWER SECTION
        do z = 2,zmd
            term1 = (c1(rm,z)-c1(rmd,z))/(2*dr*(rmd+0.5)*dr)
            term2 = ((c1(rm,z)-(2*c1(rm,z))+c1(rmd,z))/(dr**2))
            term3 = ((c1(rm,z+1)-(2*c1(rm,z))+c1(rm,z-1))/(dz**2))     ! d2c/dz2
            c2(rm,z) = c1(rm,z) + (dt * diff * ( term1 + term2 +term3 ))
        enddo
        ! FLUID -- rm,zm
        term1 = (c1(rm,zm)-c1(rmd,zm))/(2*dr*(rmd+0.5)*dr)
        term2 = ((c1(rm,zm)-(2*c1(rm,zm))+c1(rmd,zm))/(dr**2))
        term3 = ((c1(rm,zm)-(2*c1(rm,zm))+c1(rm,zmd))/(dz**2))     ! d2c/dz2
        c2(rm,zm) = c1(rm,zm) + (dt * diff * ( term1 + term2 +term3 ))
        ! FLUID —- disk above MadiCube
        ! z = (zmd+1) = (zm-1)
        do r = 2,rmd
            term1 = (c1(r+1,zm)-c1(r-1,zm))/(2*dr*(rmd+0.5)*dr)
            term2 = ((c1(r+1,zm)-(2*c1(r,zm))+c1(r-1,zm))/(dr**2))
            term3 = ((c1(rm,zm)-(2*c1(rm,zm))+c1(rm,zmd))/(dz**2))     ! d2c/dz2
            c2(r,zm) = c1(r,zm) + (dt * diff * ( term1 + term2 +term3 ))
        enddo
        ! FLUID, boundary edge -- 1,zm
        ! z = (zmd+1),(zm-1)
        term1 = (c1(2,zm)-c1(1,zm))/(2*dr*(rmd+0.5)*dr)
        term2 = ((c1(2,zm)-(2*c1(1,zm))+c1(1,zm))/(dr**2))
        term3 = ((c1(1,zm)-(2*c1(1,zm))+c1(1,zmd))/(dz**2))     ! d2c/dz2
        c2(1,zm) = c1(1,zm) + (dt * diff * ( term1 + term2 +term3 ))

        ! MIXING ROUTINE OUTSIDE MADIDROP
        ! This routine substitutes the average concentration (volume shell weighted)
        ! outside the MadiDrop and reassigns it to the entire volume outside the
        ! MadiDrop.  This simulates instantaneous mixing outside the ceramic.
        m = 0
        v = 0
        ! lower ring
        do z = 1,zmd
            m = m + (c2(rm,z) * vb(rm))
            v = v + vb(rm)
            m = m + (c1(rm,z) * vb(rm+1))
            v = v + vb(rm+1)
        enddo
        ! disk just above ceramic
        do r = 1,rmd
            m = m + (c2(r,zm) * vb(r))
            v = v + vb(r)
        enddo
        m = m + (c1(rm,zm) * vb(rm))
        v = v + vb(rm)
        m = m + (c1(rm,zm) * vb(rm+1))
        v = v + vb(rm+1)
        ! disk at top
        do r = 1,rmd
            m = m + (c1(rm,zm) * vu(r))
            v = v + vu(r)
        enddo
        m = m + (c1(rm,zm) * vu(rm))
        v = v + vu(rm)
        m = m + (c1(rm,zm) * vu(rm+1))
        v = v + vu(rm+1)

        ! Averaging
        c_ave = m / v

        ! Reassignment
        do z = 1,zmd
            c1(rm,z) = c_ave
        enddo
        do r = 1,rm
            c1(r,zm) = c_ave
        enddo
        do r = 1,rmd
            do z = 1,zmd
                c1(r,z) = c2(r,z)
                s1(r,z) = s2(r,z)
            enddo
        enddo
!write(25,*) h, mi, c_ave ! DELETE THIS LINE!
    enddo
!    write(*,*) (h-1), 'hour(s)', mi, 'minute(s), average concentration:', c_ave
!    write(25,*) (h-1), mi, c_ave
! Point data, madi.vXXXX.pdt
!    do r = 1,rm
!        do z = 1,zm
!            write(20,*) (h-1), ',', mi, ',', r, ',', z, ',', c2(r,z), ',', s2(r,z)
!        enddo
!    enddo
    enddo
write(*,*) h, 'hour(s) ended, average concentration:', c_ave
write(25,*) h, c_ave
enddo

call cpu_time(finish)
write(40,*) 'Time: ', (finish-start), 'seconds'
elapsed = finish - start
if (elapsed < 60) then
    write(*,*) 'Simulation completed in', elapsed, 'seconds'
else
    elapsed = elapsed / 60
    if (elapsed < 60) then
        write(*,*) 'Simulation completed in', elapsed, 'minutes'
    else
        elapsed = elapsed / 60
        write(*,*) 'Simulation completed in', elapsed, 'hours'
    end if
end if

end program madicubes