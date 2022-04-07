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
! 01 Nov 2016   Changing the mass transfer function
! 23 Nov 2016   Adding second diffusion tracking for "oxidant"
! 27 Nov 2016   Codified reaction term (for comparison) using inverse root function
! 28 Nov 2016   Fixed dynamic bulk density calculation -- all previous MD runs had an incorrect B_d
! 30 Nov 2016   For naive reaction term (.re.) runs, c_i was in wrong unites (0.001 formerly) changed to 1e-8
! 30 Nov 2016   Gave up on oxidant (for now) because the consumption of the oxidant would have to be unrealistic
! 30 Nov 2016   Fixed the naive reaction term.  Was developed with ppb levels, requires brut force for units.
! 01 Dec 2016	Fixed naive reaction term.  Divided when I should have multiplied: ignore v58.
! 13 Dec 2016   Changed dndt function back to pre-ppb adjustment, much closer now.
! 21 Dec 2016   Fixed porosity to n=0.3.
! 22 Dec 2016   Fixed other parameters, developing equilibrium limit function starting at v80.
! 23 Dec 2016   Made kinetic sorption model nonlinear
! 23 Dec 2016   Removed *dr*dz from reaction term; v80 wasn't making sense.  Rerunning same parameters: v81
! 23 Dec 2014   Added exponent to equilibrium limit function.
! 28 Mar 2017   removed me exponent from sorption term and included k_s, recast gamma
! 28 Mar 2017   Changed the open file command location to update in realtime.
! 29 Mar 2017   Significantly changed the naive reaction function to match manuscript.
! 15 May 2017   Attempting parallelization

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
! v0036 02Nov2016   0.0010  0.0001  7e-1    1e-4    1.5e-1  0.6     cubes, to compare
! v0037 02Nov2016   0.0010  0.0001  7e-1    1e-4    1.5e-1  0.6     MadiDrop, to compare
! set   date        dt      dr=dz   k_d     gamma   react   tau     notes
! v0036 01Nov2016   0.0010  0.0001  7e-1    1e-4    1.5e-1  0.6     checking new mass trans func, no .pdt
! v0037 01Nov2016   0.0010  0.0001  7e-1    1e-4    1.5e-1  0.6
! v0038 21Nov2016   0.0010  0.0001  7e-1    1e-4    1.5e-1  0.6     Just checking the sorption
! v0039 22Nov2016   0.0010  0.0001  7e-1    1e-4    1.5e-1  0.6     Second day run, sorption sites from day 1, v0038
! set   date        dt      dr=dz   k_d     gamma   oxida   tau     notes
! v0040 24Nov2016   0.0010  0.0001  7e-1    1e-4    1e-3    0.6     Pilot test for new oxidation term
! v0041 25Nov2016   0.0010  0.0001  7e-1    1e-4    1e-2    0.6     Testing DO levels and Ox rates - too low
! v0042 25Nov2016   0.0010  0.0001  7e-1    1e-4    5e-1    0.6     Testing DO levels and Ox rates - too low
! v0043 25Nov2016   0.0010  0.0001  7e-1    1e-4    1e-1    0.6     low, but approaching the correct curve
! v0044 25Nov2016   0.0010  0.0001  7e-1    1e-4    5e-1    0.6     Testing DO levels and Ox rates
! v0045 25Nov2016   0.0010  0.0001  7e-1    1e-4    1       0.6     Testing DO levels and Ox rates
! v0046 25Nov2016   0.0010  0.0001  7e-1    1e-4    5       0.6     Testing DO levels and Ox rates
! v0047 25Nov2016   0.0010  0.0001  7e-1    1e-4    1e+1    0.6     Testing DO levels and Ox rates
! v0048 26Nov2016   0.0010  0.0001  7e-1    1e-4    1e+2    0.6     Changed DO to 8e-3, 24 hour
! set   date        geo dt      dr=dz   tau     k_d     gamma   ox      a       b       notes
! v49   28Nov2016   mc  0.0010  0.0001  0.6     7e-1    1e-4            1.5dt   -0.2    local, low, but good
! v50   28Nov2016   md  0.0010  0.0001  0.6     7e-1    1e-4            1.5dt   -0.2    remote, still running
! v51   28Nov2016   mc  0.0010  0.0001  0.6     7e-1    1e-4    1e+2                    NOT RUN
! v52   28Nov2016   md  0.0010  0.0001  0.6     7e-1    1e-4    1e+2                    NOT RUN
! v53   27Nov2016   mc  0.0010  0.0001  0.6     7e-1    1e-4    1e+2                    changed ratio of DO/Ag, but still linear-looking Ag increase
! v54   29Nov2016   mc  0.0010  0.0001  0.6     7e-1    1e-4    1e+5                    original ox ratio, increased rate
! v55   28Nov2016   mc  0.0010  0.0001  0.6     7e-4    1e-4            52.96dt -0.2    attempting to match a (local, a=5196dt, remote, a=51.96dt), still way too low, appears to be dropping from seed c_i
! set   date        geo dt      dr=dz   tau     k_d     gamma   a       b       ci      si      notes
! v56   30Nov2016   mc  0.0010  0.0001  0.6     7e-4    1e-4    52.96dt -0.2    1e-8    0       fixed c_i
! v57   30Nov2016   mc  0.0010  0.0001  0.6     7e-4    1e-4    52.96dt -0.2    1e-8    5e-9    w/ c_i & s_i
! v57 was an error, the beginning values were smaller than the initial concentration.  This was because the 
! reaction function was designed for ppb values and not kg m^-3 values.  The result was catastrophic.  I am 
! surprised it didn't create an unstable mess; actually, it might have.  The naive reaction model was changed 
! and should produce a stable result in v58.  v59 will be used to develop the comparison run for MadiCubes and,
! hopefully, v60 will be the last MadiDrop run.  Still to do, sensitivity analysis; however, we know the rough 
! results of that.
! v58   30Nov2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    52.96dt -0.2    1e-9    5e-10   w/ corrected dr
! v59   30Nov2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    52.96dt -0.2    1e-9    5e-10   bug fixed
! Still having problems... examining in Matlab, v60
! v60   12Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    52.96dt -0.2    1e-9    5e-10   check each dt
! v61   13Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    15dt    -0.2    1e-9    5e-10   dndt, 2 hrs
! v62   13Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    150dt   -0.2    1e-9    5e-10
! v63   13Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    15      -0.2    1e-9    5e-10
! v64   14Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    14      -0.2    1e-9    5e-10   better, 24hr
! v65   14Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    15      -0.2    1e-9    5e-10
! v66   14Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    16      -0.2    1e-9    5e-10
! v67   14Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    17      -0.2    1e-9    5e-10   Wrong shape
! v68   15Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    20      -0.2    1e-9    5e-10   Too low
! v69   15Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    15      -0.9    1e-9    5e-10   Way too high
! v70   19Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    15      -0.7    1e-9    5e-10   8 hour
! v71   19Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    15      -0.5    1e-9    5e-10   8 hour
! v72   19Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    15      -0.3    1e-9    5e-10
! v73   20Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    10      -0.3    1e-9    5e-10
! v74   20Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    1       -0.4    1e-9    5e-10
! v75   20Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    0.5     -0.4    1e-9    5e-10
! v76   20Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    0.1     -0.4    1e-9    5e-10
! v77   21Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    5       -0.35   1e-9    5e-10
! v78   22Dec2016   mc  0.0010  0.0001  0.6     5e-1    5e-4    15      -0.1    1e-9    5e-10
! Fix parameters, tau=0.66, k_d=0.28
! To analyze:
! new set of gamma and M (exponent) to fit sorption
! maybe a new reaction function (1000-C)
! set   date        geo dt      dr=dz   tau     k_d     gamma   M   a       b       eqr ci      si      time
! v79   22Dec2016   mc  0.0010  0.0001  0.66    0.28    5e-4    1   0       0       0   1.14e-1 5e-10   8
! v80   22Dec2016   mc  0.0010  0.0001  0.66    0.28    5e-4    1   0       0       0.1  0      5e-10   8
! v81   22Dec2016   mc  0.0010  0.0001  0.66    0.28    5e-4    1   0       0       0.1  0      5e-10   8
! Altering the gamma sorption rate parameter based on results from v79 and the eqr oxidation/equilibrium rate
! based on v81:
! v82   23Dec2016   mc  0.0010  0.0001  0.66    0.28    2e-4    1   0       0       0.01 0      5e-10   8
! set   date        dt      dr=dz   tau     k_d     gamma   M   a   b   eqr     eqe ci      si      time
! v83   23Dec2016   0.0010  0.0001  0.66    0.28    2e-4    1   0   0   1d-3    1.2 0       5e-10   8
! v84/5r23Dec2016   0.0010  0.0001  0.66    0.28    8e-5    1   0   0   0       1.2 6.95d-3 5e-10   8
! eqr =/= 0 => NaN error.  Instead, hard-code term5=0.
! v86   23Dec2016   0.0010  0.0001  0.66    0.28    8e-5    1.2 0   0   0       1.2 6.95d-3 5e-10   8
! v87   23Dec2016   0.0010  0.0001  0.66    0.28    8e-5    0.8 0   0   0       1.2 6.95d-3 5e-10   8   unstab
! v87 failed, v86 still running
! v88   11Jan2017   0.0010  0.0001  0.66    0.28    3.5e-4  1.1 0   0   0       1.2 6.95d-3 5e-10   12  unstab
! v89   11Jan2017   0.0010  0.0001  0.66    0.28    3.5e-4  1.1 0   0   1d-3    .25 0       5e-10   12  unstab
! v90   11Jan2017   0.0010  0.0001  0.66    0.28    3.5e-4  1.1 0   0   1d-3    .50 0       5e-10   12  unstab
! v91   13Jan2017   0.0010  0.0001  0.66    0.28    6.5e-4  1   0   0   0       -   6.95d-3 5e-10   12
! v92   14Jan2017   0.0010  0.0001  0.66    0.35    6.5e-4  1   0   0   0       -   6.95d-3 5e-10   12
! v93   14Jan2017   0.0010  0.0001  0.66    0.35    1.5e-3  1   0   0   0       -   6.95d-3 5e-10   12
! v94   14Jan2017   0.0010  0.0001  0.66    0.55    1.5e-3  1   0   0   0       -   6.95d-3 5e-10   12
! v95   14Jan2017   0.0010  0.0001  0.66    0.65    1.5e-3  1   0   0   0       -   6.95d-3 5e-10   16
! v96   23Mar2017   0.0010  0.0001  0.65    0.30    1.5e-2  1   0   0   0       -   6.95d-3 5e-10   8
! v97   24Mar2017   0.0010  0.0001  0.65    0.30    1.5e-1  1   0   0   0       -   6.95d-3 5e-10   8
! v98   24Mar2017   0.0010  0.0001  0.65    0.30    5       1   0   0   0       -   6.95d-3 5e-10   8   unstab
! v99   28Mar2017   0.0010  0.0001  0.65    2.28    1.5e-2  1   0   0   0       -   6.95d-3 5e-10   8
! set   date        dt      dr=dz   tau     k_d     k_s     gamma   a       b       ci      si      time
! w01   28Mar2017   0.0010  0.0001  0.65    2.28    0.1     0.008   ?       ?       6.95d-3 5e-10   8
! w02   28Mar2017   0.0010  0.0001  0.65    2.28    0.1     0.008   ?       ?       1.13d-2 5e-10   8   invalid
! w03   28Mar2017   0.0010  0.0001  0.65    2.28    0.1     0.001   ?       ?       6.95d-3 5e-10   8
! w04   28Mar2017   0.0010  0.0001  0.65    2.28    0.1     0.001   4.039   0.631   1d-12   1e-12   24
! w05   30Mar2017   0.0010  0.0001  0.65    2.28    0.1     0.0001  0       0       6.95d-3 1e-12   8
! w06   30Mar2017   0.0010  0.0001  0.65    2.28    0       0.001   0       0       6.95d-3 1e-12   8
! w07   30Mar2017   0.0010  0.0001  0.65    2.28    0.01    0.00001 0       0       6.95d-3 1e-12   8
! w08   30Mar2017   0.0010  0.0001  0.65    2.28    0.01    0.00003 0       0       6.95d-3 1e-12   8
! w09   30Mar2017   0.0010  0.0001  0.65    2.28    0.01    0.00005 0       0       6.95d-3 1e-12   8
! w10   30Mar2017   0.0010  0.0001  0.65    2.28    0.01    0.00007 0       0       6.95d-3 1e-12   8
! w11   31Mar2017   0.0010  0.0001  0.65    2.28    0.01    0.00002 0       0       6.95d-3 1e-12   24
! w12   02Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  0       0       6.95d-3 1e-12   24
! w13   02Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  0       0       1.14d-1 1e-12   24
! w14   02Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  4.039   0.631   1d-12   1e-12   8
! w15   09Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  5d-3    0.631   1d-12   1e-12   8
! w16   11Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  4.1d-3  0.631   1d-12   1e-12   24
! w17   11Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  4.0d-3  0.631   1d-12   1e-12   24
! w18   11Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.9d-3  0.631   1d-12   1e-12   24
! w19   13Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.5d-3  0.720   1d-12   1e-12   24
! w20   18Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.5d-3  0.550   1d-12   1e-12   24
! w21   21Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.5d-3  0.600   1d-12   1e-12   24
! w22   21Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.75d-3 0.600   1d-12   1e-12   24
! w23   21Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.0d-3  0.631   1d-12   1e-12   24
! w24   24Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.5d-3  0.631   1d-12   1e-12   24
! w25   27Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.4d-3  0.650   1d-12   1e-12   24
! Variations: tortuosity
! t04   27Apr2017   0.0010  0.0001  0.40    2.28    0.01    1.8d-5  3.4d-3  0.650   1d-12   1e-12   8
! t10   27Apr2017   0.0010  0.0001  1.00    2.28    0.01    1.8d-5  3.4d-3  0.650   1d-12   1e-12   8
! Variations: sorption (k_d)
! s02   27Apr2017   0.0010  0.0001  0.65    0.228   0.01    1.8d-5  3.4d-3  0.650   1d-12   1e-12   8
! s22   27Apr2017   0.0010  0.0001  0.65    228     0.01    1.8d-5  3.4d-3  0.650   1d-12   1e-12   8
! Variations: reaction
! r-4   27Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.4d-4  0.650   1d-12   1e-12   8
! r-2   27Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.4d-2  0.650   1d-12   1e-12   8
! Variations: "drop" geometry
! d01   27Apr2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.4d-4  0.650   1d-12   1e-12   24
! dp    15May2017   0.0010  0.0001  0.65    2.28    0.01    1.8d-5  3.4d-4  0.650   1d-12   1d-12   8 parallel

implicit none
integer p
parameter(p=kind(1.0D0))
real (p) :: r_bucket, r_md, h_bucket, h_md, v_bucket, v_md
real (p) :: dr, dz, dt
real (p) ::b_d, n, tau, k_d, k_s, gamma, me, a, b, eqr, eqe
real (p) :: diff, diff_md
real start, finish, elapsed
integer tm, hr, rm, zm, rmd, zmd
integer r, z, t, h, mi
real (p) :: ci, si
real (p) :: m, v, v_water, c_ave
real (p) :: term1, term2, term3, term4, term5
!********** FOR CUBE GEOMETERY **********
!real (p) :: c1(57,101)
!real (p) :: c2(57,101)
!real (p) :: s1(57,101)
!real (p) :: s2(57,101)
!real (p) :: vb(58), vu(58) ! this should have one greater than the other dimension to account for the skin layer
!********** FOR DROP GEOMETERY **********
real (p) :: c1(282,401)
real (p) :: c2(282,401)
real (p) :: s1(282,401)
real (p) :: s2(282,401)
real (p) :: vb(283), vu(283) ! this should have one greater than the other dimension to account for the skin layer
!********** END GEOMETRY **********

call cpu_time(start)

open(20, file='md.pdt')
open(25, file='md.adt')
!open(40, file='current.par')

!****************************************************************************
!*                                                                          *
!*                               GEOMETRY AND TIME                          *
!*                                                                          *
!****************************************************************************

dr = 0.0001                 ! m, dr = 0.1 mm yields 56 nodes in the MadiDrop plus two edge nodes
dz = 0.0001                 ! m, dz = 0.1 mm yields 100 nodes in the MadiDrop plus two edge nodes
dt = 0.0010                 ! s
tm = 60000                  ! number of timesteps in one minute
hr = 7                      ! hours of simulation

!********** FOR CUBE GEOMETERY **********
!rm = 57                     ! number of nodes in the radial direction in the domain
!zm = 101                    ! number of nodes in the z-direction in the domain
!rmd = 56                    ! number of nodes in the radial direction within MadiDrop
!zmd = 100                   ! number of nodes in the z-direction within MadiDrop
! The "bucket" size is set to model 1/100th of a 10 L bucket as a tall cylinder.  This model assumes
! 100 adjacent and identical cylinders that will make up the 10 L.
!r_bucket = 0.0150           ! m, modified radius of bucket to reflect many MadiDrop pieces (summer
! 2015 form factor).  Radius = 1 cm as a perfectly mirrored column in
! bucket.
!h_bucket = 0.1415           ! m, height of bucket or reactor
!********** FOR DROP GEOMETERY **********
rm = 282                    ! number of nodes in the radial direction in the domain (was 57)
zm = 401                    ! number of nodes in the z-direction in the domain (was 101)
rmd = 281                   ! for full MadiDrop: r_md should be about 3.2cm (was 56)
zmd = 400                   ! for full MadiDrop: z_md should be about 3cm (was 100)
r_bucket = 0.1500           ! m, radius of model bucket
h_bucket = 0.1415           ! m, height of bucket or reactor
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
me=1                                ! exponent

!react = 1.5e-1                ! kg m^-3 s^-1  OLD PARAMETER
a=3.40d-3                           ! reaction parameter
b=0.650                             ! reaction exponent

eqr=1d-3                            ! equilibrium oxidataion rate parameter (constant multiplier)
eqe=0.50                            ! equilibrium oxidataion rate exponent

! Oxidation rate corrections of: 15999 and 215.736 was tested in v53...
! Equilibrium oxidation ratio of consume O2(aq) at 15.999 for every Ag(aq) produced at 215.736
! based on Liu, J., and R. H. Hurt (2010), Ion release kinetics and particle persistence in aqueous nano-silver colloids., Environ. Sci. Technol., 44(6), 2169–75, doi:10.1021/es9035557.

tau = 0.65 ! empirical value from tritium experiments.  Speculate that a
! higher tau is related to the pore structures in the ceramic having
! greater connectivity than a uniformly distributed homogeneous pore
! structure (i.e., soil).

b_d = 1319.5256                             ! kg m^-3, bulk density (mass density of dry ceramic)
n = 0.3                                     ! dimensionless, porosity
diff = 1.5d-9                               ! m^2 s^-1, value from literature:
! Johnston, R. R. M., and M. Spiro (1967), Diffusion coefficients of the silver ion and the disulfitosilver (I) ion by the rotating disk method, J. Phys. Chem., 71(12), 3784–3790, doi:10.1021/j100871a011.
diff_md = tau * diff                        ! m^2 s^-1, adjusted for ceramic.

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
ci=7.83840d-03 ! c_i = 1.14d-1 or 6.93d-3 kg m^-3
!                10 mg/L (ppm) = 10^-2 kg m^-3, 10 ppb = 10^-5 kg m^-3
si=6.92540d-03 ! used 3e-5 to model day 2, should be rechecked.
do r=1,rmd
    do z=1,zmd
        c1(r,z)=ci
        s1(r,z)=si
    enddo
enddo

! BULK WATER (outside the ceramic)
! for the kinetic sorption trials, we start with 114 ppm = 1.14e-1 kg m^-3 for ci everywhere.
ci=1.11818d-03
! for all simulations, ci=0
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
do h = 7,hr
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
        enddo
        do r=1,rm
            c1(r,zm) =c_ave
        enddo
        ! Pore water
        do r=1,rmd
            do z=1,zmd
                c1(r,z)=c2(r,z)
                s1(r,z)=s2(r,z)
            enddo
        enddo
    enddo
    enddo
write(*,*) h, 'hour(s) ended'
write(*,*) 'average concentration:', c_ave
write(25,*) h, c_ave    !, ox2(25,50)

enddo
! Point data, vXXXX.pdt
do r = 1,rm
    do z = 1,zm
        ! This is the node data for the last time point
        write(20,*) r, ',', z, ',', c2(r,z), ',', s2(r,z)!, ',', ox2(r,z)
    enddo
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
