module sorb_react
contains

! ********** SORPTION MODEL **********
function dsdt(con,sor,k_d,k_s,gamma)
integer p
parameter(p=kind(1.0D0))
real (p) :: con, sor, dsdt
! k_d is the "distribution coefficient" and is the equilibrium constant defining sorption
! gamma is the rate parameter, and also controls the reversability
real (p) :: k_d, k_s, gamma
! reversable nonlinear kinetic sorption from Fetter, C. W. (2008), Contaminant Hydrogeology, Waveland Press, Inc., Long Grove, Il.
dsdt=gamma*(k_d*con-k_s*sor)
end function

! ********** NAIVE REACTION MODEL **********
function dndt(con,a,b,dt)
integer p
parameter(p=kind(1.0D0))
real (p) :: dndt
real (p) :: con, c1, c0, a, b, dt, t
c0=con
t=(c0/a)**(1/b)
t=t+dt
c1=a*(t**b)
dndt=(c1-c0)/dt
end function

end module

program madi

use sorb_react

! This program solves the diffusion equation, as modified, by the Euler scheme with a second-order finite
! difference model in an axi-symetric geometry for a bucket for a cylindrical representation of the ceramic
! cubes or the original silver-empregnated disk in the center of the bottom of the bucket.

! Full model in cylindrical coordinates
! Axi-symetric around z = 0
! Assuming radial symmetry d/dtheta=0

! REVISION HISTORY
! 04 Nov 2017   Program uploaded for reference to manuscript

implicit none
integer p
parameter(p=kind(1.0D0))
real (p) :: r_bucket, r_md, h_bucket, h_md, v_bucket, v_md
real (p) :: dr, dz, dt
real (p) ::b_d, n, tau, k_d, k_s, gamma, a, b
real (p) :: diff, diff_md
real start, finish, elapsed
integer tm, hr, rm, zm, rmd, zmd
integer r, z, t, h, mi
real (p) :: ci, si, oi
real (p) :: m, v, v_water, c_ave
real (p) :: term1, term2, term3, term4, term5
!********** FOR CUBE GEOMETERY **********
real (p) :: c1(57,101)
real (p) :: c2(57,101)
real (p) :: s1(57,101)
real (p) :: s2(57,101)
real (p) :: ox1(57,101)
real (p) :: ox2(57,101)
real (p) :: vb(58), vu(58) ! this should have one greater than the other dimension to account for the skin layer
!********** FOR DROP GEOMETERY **********
!real (p) :: c1(282,401)
!real (p) :: c2(282,401)
!real (p) :: s1(282,401)
!real (p) :: s2(282,401)
!real (p) :: ox1(282,401)
!real (p) :: ox2(282,401)
!real (p) :: vb(283), vu(283) ! this should have one greater than the other dimension to account for the skin layer
!********** END GEOMETRY **********

call cpu_time(start)

open(20, file='madi.pdt')
open(25, file='madi.adt')
open(40, file='madi.par')

!****************************************************************************
!*                                                                          *
!*                               GEOMETRY AND TIME                          *
!*                                                                          *
!****************************************************************************

dr = 0.0001                 ! m, dr = 0.1 mm yields 56 nodes in the MadiDrop plus two edge nodes
dz = 0.0001                 ! m, dz = 0.1 mm yields 100 nodes in the MadiDrop plus two edge nodes
dt = 0.0010                 ! s
tm = 60000                  ! number of timesteps in one minute; validated diffusion in a narrow rod,
                            ! see: Kahler, D. M. (2011), The Acceleration of the Diffusion-Limited Pump-and-Treat Aquifer Remediation with Pulsed Pumping that Generates Deep Sweeps and Vortex Ejections in Dead-End Pores, Duke University. http://dukespace.lib.duke.edu/dspace/handle/10161/3915
hr = 24                     ! hours of simulation

!********** FOR CUBE GEOMETERY **********
rm = 57                     ! number of nodes in the radial direction in the domain
zm = 101                    ! number of nodes in the z-direction in the domain
rmd = 56                    ! number of nodes in the radial direction within MadiDrop
zmd = 100                   ! number of nodes in the z-direction within MadiDrop
! The "bucket" size is set to model 1/100th of a 10 L bucket as a tall cylinder.  This model assumes
! 100 adjacent and identical cylinders that will make up the 10 L.
r_bucket = 0.0150           ! m, modified radius of bucket to reflect many MadiDrop pieces (summer
! 2015 form factor).  Radius = 1 cm as a perfectly mirrored column in
! bucket.
h_bucket = 0.1415           ! m, height of bucket or reactor
!********** FOR DROP GEOMETERY **********
!rm = 282                    ! number of nodes in the radial direction in the domain (was 57)
!zm = 401                    ! number of nodes in the z-direction in the domain (was 101)
!rmd = 281                   ! for full MadiDrop: r_md should be about 3.2cm (was 56)
!zmd = 400                   ! for full MadiDrop: z_md should be about 3cm (was 100)
!r_bucket = 0.1500           ! m, radius of model bucket
!h_bucket = 0.1415           ! m, height of bucket or reactor
!********** END GEOMETRY **********

r_md = rmd*dr               ! m, radius of the MadiDrop
h_md = zmd*dz               ! m, height of MadiDrop


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

k_d = 2.28                          ! m^3 kg^-1; 1 L kg^-1 = 10^-3 m^3 kg^-1 distribution coefficent
k_s = 0.01                          ! k_s gamma == reverse sorption parameter
gamma = 0.000018                    ! s^-1 rate parameter

a=3.00d-3                           ! reaction parameter
b=0.631                             ! reaction exponent

tau = 0.76 ! empirical value from tritium experiments.  Speculate that a
! higher tau is related to the pore structures in the ceramic having
! greater connectivity than a uniformly distributed homogeneous pore
! structure (i.e., soil).

b_d = 1319.5256                             ! kg m^-3, bulk density (mass density of dry ceramic)
n = 0.3                                     ! dimensionless, porosity
diff = 1.5d-9                               ! m^2 s^-1, value from literature:
! Johnston, R. R. M., and M. Spiro (1967), Diffusion coefficients of the silver ion and the disulfitosilver (I) ion by the rotating disk method, J. Phys. Chem., 71(12), 3784–3790, doi:10.1021/j100871a011.
diff_md = tau * diff                        ! m^2 s^-1, adjusted for ceramic.

write(*,*) 'Geometry and physical parameters loaded'

write(40,*) 'Radius of "bucket,r_bucket =   ', r_bucket
write(40,*) 'r_md =                         ', r_md
write(40,*) 'h_bucket =                     ', h_bucket
write(40,*) 'h_md =                         ', h_md
write(40,*) 'dr =                           ', dr
write(40,*) 'dz =                           ', dz
write(40,*) 'dt =                           ', dt
write(40,*) 'rmd =                          ', rmd
write(40,*) 'zmd =                          ', zmd
write(40,*) 'Time steps in a second tm =    ', tm
write(40,*) 'Total simulation time hr =     ', hr
write(40,*) 'DIFFUSION'
write(40,*) 'tau =                          ', tau
write(40,*) 'diffusion (MadiDrop) =         ', diff_md
write(40,*) 'SORPTION'
write(40,*) 'B_d =                          ', b_d
write(40,*) 'k_d =                          ', k_d
write(40,*) 'gamma =                        ', gamma
write(40,*) 'REACTION'
write(40,*) 'reaction rate =                ', a ! this is the reaction multiplier
write(40,*) 'reaction exponent =            ', b ! this it the reaction exponent
write(40,*) 'Lower volume elements'
do r = 1,rm
    write(40,*) 'Element                    ', r, vb(r)
enddo
write(40,*) 'Upper volume elements'
do r = 1,rm
    write(40,*) 'Element                    ', r, vu(r)
enddo
write(40,*) 'v_bucket                       ', v_bucket
write(40,*) 'v_MD                           ', v_md
write(40,*) 'v_water                        ', v_water
write(40,*) 'Output files'
write(40,*) 'madi.pdt     concentration matrix'
write(40,*) 'headers:              hour, radius, height, concentration, sorption'
write(40,*) 'madi.adt     concentration average time series'
write(40,*) 'headers:              hour, average concentration'
write(40,*) 'madi.par     parameter list'

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

! NOTE: C(r,z) where c1 and c2 are progressive timesteps

! PORE WATER (inside of ceramic -- this is especially useful for multiple days)
ci=1d-12 ! c_i = 1.14d-1 or 6.93d-3 kg m^-3
!                10 mg/L (ppm) = 10^-2 kg m^-3, 10 ppb = 10^-5 kg m^-3
si=1e-12 ! used 3e-5 to model day 2, should be rechecked.
do r=1,rmd
    do z=1,zmd
        c1(r,z)=ci
        s1(r,z)=si
    enddo
enddo

! BULK WATER (outside the ceramic)
! for the kinetic sorption trials, we start with 114 ppm = 1.14e-1 kg m^-3 for ci everywhere.
ci=1d-12
si=0
do r=1,rm
    c1(r,zm)=ci
    s1(r,zm)=si
enddo
do z=1,zm
    c1(rm,z)=ci
    s1(rm,z)=si
enddo


write(*,*) 'Model specifications loaded and arrays initialized'

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
        term4 = dsdt(c1(1,1),s1(1,1),k_d,k_s,gamma)
        term5=dndt(c1(1,1),a,b,dt)
        c2(1,1)=c1(1,1)+(dt*(diff_md*(term1+term2+term3)-(term4*b_d/n)+term5))
        s2(1,1) = s1(1,1) + (dt * term4)

        ! boundary edge!	AXIS OF SYMMETRY
        do z = 2,zmd
            ! symmetry boundary condition, phantom points created such that C(r-1) = C(1)
            term1 = (c1(2,z)-c1(1,z))/(dr**2) ! note that the denomenator was simplified
            term2 = ((c1(2,z)-(2*c1(1,z))+c1(1,z))/(dr**2))
            term3 = ((c1(1,z+1)-(2*c1(1,z))+c1(1,z-1))/(dz**2))
            term4 = dsdt(c1(1,z),s1(1,z),k_d,k_s,gamma)
            term5=dndt(c1(1,z),a,b,dt)
            c2(1,z) = c1(1,z)+(dt*(diff_md*(term1+term2+term3)-(term4*b_d/n)+term5))
            s2(1,z) = s1(1,z) + (dt * term4)
        enddo

        ! boundary edge!	BOTTOM OF BUCKET
        do r = 2,rmd
            ! no-flux boundary condition, phantom point created such that C(z-1) = C(z+1)
            term1 = (c1(r+1,1)-c1(r-1,1))/(2*dr*(r-0.5)*dr)
            term2 = (c1(r+1,1)-(2*c1(r,1))+c1(r-1,1))/(dr**2)
            term3 = ((c1(r,2)-(2*c1(r,1))+c1(r,1))/(dz**2))
            term4 = dsdt(c1(r,1),s1(r,1),k_d,k_s,gamma)
            term5=dndt(c1(r,1),a,b,dt)
            c2(r,1) = c1(r,1)+(dt*(diff_md*(term1+term2+term3)-(term4*b_d/n)+term5))
            s2(r,1) = s1(r,1) + (dt * term4)
        enddo

        ! CERAMIC INTERIOR —- NON-BOUNDARY
        do r = 2,rmd
            do z = 2,zmd
                ! Ag DIFFUSION
                ! Product rule of the first term gives us term1 and term2
                term1 = (c1(r+1,z)-c1(r-1,z))/(2*dr*(r-0.5)*dr)             ! 1/r dc/dr
                term2 = (c1(r+1,z)-(2*c1(r,z))+c1(r-1,z))/(dr**2)           ! d2c/dr2
                ! there is no theta dependance, all terms with d/dtheta go to zero
                term3 = ((c1(r,z+1)-(2*c1(r,z))+c1(r,z-1))/(dz**2))         ! d2c/dz2
                term4 = dsdt(c1(r,z),s1(r,z),k_d,k_s,gamma)                  ! sorption term
                term5=dndt(c1(r,z),a,b,dt)                                 ! equilibrium limit term
                c2(r,z) = c1(r,z)+(dt*(diff_md*(term1+term2+term3)-(term4*b_d/n)+term5))
                s2(r,z) = s1(r,z) + (dt * term4)
            enddo
        enddo

        ! OUTSIDE THE MadiDrop
        ! boundary edge!	BOTTOM OF BUCKET
        ! no-flux boundary condition, phantom point created such that C(z-1) = C(z)
        term1 = (c1(rm,1)-c1(rmd,1))/(2*dr*(rmd+0.5)*dr)
        term2 = ((c1(rm,1)-(2*c1(rm,1))+c1(rmd,1))/(dr**2))
        term3 = ((c1(rm,2)-(2*c1(rm,1))+c1(rm,1))/(dz**2))
        c2(rm,1) = c1(rm,1) + (dt * diff * ( term1 + term2 +term3 ))

        ! FLUID —- LOWER SECTION
        do z = 2,zmd
            term1 = (c1(rm,z)-c1(rmd,z))/(2*dr*(rmd+0.5)*dr)
            term2 = ((c1(rm,z)-(2*c1(rm,z))+c1(rmd,z))/(dr**2))
            term3 = ((c1(rm,z+1)-(2*c1(rm,z))+c1(rm,z-1))/(dz**2))
            c2(rm,z) = c1(rm,z) + (dt * diff * ( term1 + term2 +term3 ))
        enddo

        ! FLUID -- rm,zm
        term1 = (c1(rm,zm)-c1(rmd,zm))/(2*dr*(rmd+0.5)*dr)
        term2 = ((c1(rm,zm)-(2*c1(rm,zm))+c1(rmd,zm))/(dr**2))
        term3 = ((c1(rm,zm)-(2*c1(rm,zm))+c1(rm,zmd))/(dz**2))
        c2(rm,zm) = c1(rm,zm) + (dt * diff * ( term1 + term2 +term3 ))

        ! FLUID —- disk above MadiCube
        ! z = (zmd+1) = (zm-1)
        do r = 2,rmd
            term1 = (c1(r+1,zm)-c1(r-1,zm))/(2*dr*(rmd+0.5)*dr)
            term2 = ((c1(r+1,zm)-(2*c1(r,zm))+c1(r-1,zm))/(dr**2))
            term3 = ((c1(rm,zm)-(2*c1(rm,zm))+c1(rm,zmd))/(dz**2))
            c2(r,zm) = c1(r,zm) + (dt * diff * ( term1 + term2 +term3 ))
        enddo

        ! FLUID, boundary edge -- 1,zm
        ! z = (zmd+1),(zm-1)
        term1 = (c1(2,zm)-c1(1,zm))/(2*dr*(rmd+0.5)*dr)
        term2 = ((c1(2,zm)-(2*c1(1,zm))+c1(1,zm))/(dr**2))
        term3 = ((c1(1,zm)-(2*c1(1,zm))+c1(1,zmd))/(dz**2))
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
        c_ave=m/v

        ! Reassignment
        ! Outer water
        do z=1,zmd
            c1(rm,z)=c_ave
            ox1(rm,z)=oi ! conserve saturation outer water
        enddo
        do r=1,rm
            c1(r,zm) =c_ave
            ox1(r,zm)=oi
        enddo
        ! Pore water
        do r=1,rmd
            do z=1,zmd
                c1(r,z)=c2(r,z)
                s1(r,z)=s2(r,z)
                ox1(r,z)=ox2(r,z)
            enddo
        enddo
    enddo
enddo

write(*,*) h, 'hour(s) ended'
write(*,*) 'average concentration:', c_ave
write(25,*) h, c_ave

enddo
! Point data, vXXXX.pdt
!do r = 1,rm
!    do z = 1,zm
!        ! This is the node data for the last time point
!        write(20,*) r, ',', z, ',', c2(r,z), ',', s2(r,z)!, ',', ox2(r,z)
!    enddo
!enddo

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

end program madi
