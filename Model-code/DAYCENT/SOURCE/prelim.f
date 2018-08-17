
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine prelim

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc' 
      include 't0par.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ... Initialize variables and parameters

c ... Function declarations
      double precision catanf, line
      external  catanf, line

c ... Local variables
      integer   ii, lyr
      double precision dely, delx, fcbd(6), fccl(6), fcin(6), fcom(6),
     &          fcsa(6), fcsi(6), fcwp(6), ompc, xslope,
     &          textur, wpbd(6), wpcl(6), wpin(6),
     &          wpom(6), wpsa(6), wpsi(6), wpwp(6), yint
      double precision soildepth, adepsum

c ... swflag lets the model user choose between using actual data 
c ... for awilt and afiel or equations from Gupta and Larson (1979) 
c ... or Rawls et al (1982).
c ...
c ... swflag=0 Use actual data
c ... swflag=1 Use G&L for both awilt (-15 bar) and afiel (-0.33 bar)
c ... swflag=2 Use G&L for both awilt (-15 bar) and afiel (-0.10 bar)
c ... swflag=3 Use Rawls for both awilt (-15 bar) and afiel (-0.33 bar)
c ... swflag=4 Use Rawls for both awilt (-15 bar) and afiel (-0.10 bar)
c ... swflag=5 Use Rawls for afiel (-0.33 bar) and actual data for awilt
c ... swflag=6 Use Rawls for afiel (-0.10 bar) and actual data for awilt
c ...
c ...     swflag   1          2          3        4       5       6
      data fcsa / 0.3075,    0.5018,   -0.20,   -0.30,  -0.19,   0.31/
      data fcsi / 0.5886,    0.8548,    0.0,     0.0,    0.0,    0.0/
      data fccl / 0.8039,    0.8833,    0.36,    0.23,   0.0,    0.0/
      data fcom / 2.208E-03, 4.966E-03, 0.0299,  0.0317, 0.0210, 0.026/
      data fcbd /-0.1434,   -0.2423,    0.0,     0.0,    0.0,    0.0/
      data fcwp / 0.0,       0.0,       0.0,     0.0,    0.72,   0.41/
      data fcin / 0.0,       0.0,       0.2576, 0.4118, 0.2391, 0.4103/
      data wpsa /-0.0059,   -0.0059,    0.0,     0.0,    0.0,    0.0/
      data wpsi / 0.1142,    0.1142,    0.0,     0.0,    0.0,    0.0/
      data wpcl / 0.5766,    0.5766,    0.50,    0.50,   0.0,    0.0/
      data wpom / 2.228E-03, 2.228E-03, 0.0158,  0.0158, 0.0,    0.0/
      data wpbd / 0.02671,   0.02671,   0.0,     0.0,    0.0,    0.0/
      data wpwp / 0.0,       0.0,       0.0,     0.0,    1.0,    1.0/
      data wpin / 0.0,       0.0,       0.0260,  0.0260, 0.0,    0.0/

c ... Time initializations -  time step is one month
      dt = 1.0/12.0
      time = strtyr
      month = 0
      cntstep = 0
      time2 = strtyr

c ... Allow for time step < 1 month for running decomp
c ... ntspm is the number of time steps per month for decomp
c ... (read from the fix.100 file)

c ... Initializations
      crpgrw = 0
      seedl = 0
      forgrw = 0
      falprc = 0

c ... Initialize volitalization accumulators
      volgma = 0.0
      volexa = 0.0
      volpla = 0.0
 
c ... Initialize erosion variables
      scloss = 0.0
      sclosa = 0.0
 
c ... Initialize accumulators
      call annacc

c ... Open the c14 data file 
      if (labtyp .eq. 1) then
        open(unit=10,file='c14data',status='OLD')
      endif

c ... Field capacity and wilting point.  Computations based on
c ... Gupta and Larson 1979, 'Estimating soil and water retention
c ... characteristics from particle size distribution, organic 
c ... matter percent and bulk density'. Water Resources Research 15:1633
c ... or Rawls et al (1982) 'Estimation of soil water properties'
c ... Trans. ASAE ???:1316
c ... Field capacity options of -0.1 or -0.33 bar.
c ... Wilting point assumed to be water content at -15 bars.
c ... Calculate organic matter from initial conditions, ivauto or 
c ... value at the beginning of an extend
c ... Note that Gupta and Larson and Rawls use % for texture
c ... but values here are fractions.
      soildepth = 0.0
      if (swflag .ne. 0) then
        do 11 lyr = 1, nlayer
          soildepth = soildepth + adep(lyr)
11      continue
c ..... For really deep soils, don't use rock fraction -mdh 6/29/99
        if (soildepth .gt. 150) then
          rock = 0.0
        endif
c ..... Set somsc using initial values read from the <site>.100 file
        somsc = som1ci(2,1) + som1ci(2,2) + som2ci(2,1) + som2ci(2,2) +
     &          som3ci(1) + som3ci(2)
        ompc = somsc*1.724/(10000*bulkd*edepth)
        do 10 lyr = 1, nlayer
          afiel(lyr) = fcsa(swflag)*sand  + fcsi(swflag)*silt +
     &                 fccl(swflag)*clay  + fcom(swflag)*ompc +
     &                 fcbd(swflag)*bulkd + fcwp(swflag)*awilt(lyr) +
     &                 fcin(swflag)
          awilt(lyr) = wpsa(swflag)*sand  + wpsi(swflag)*silt +
     &                 wpcl(swflag)*clay  + wpom(swflag)*ompc +
     &                 wpbd(swflag)*bulkd + wpwp(swflag)*awilt(lyr) +
     &                 wpin(swflag)
          ompc = ompc * 0.85
c ....... modifiy afiel and awilt according to fractional volume of rock -mdh 5/27/99
          afiel(lyr) = afiel(lyr) * (1.0 - rock)
          awilt(lyr) = awilt(lyr) * (1.0 - rock)
10      continue
      endif
        
      open(unit=80,file='daily.csv')
      write(80,85) 'time', 'dayofyr', 'PET(cm)', 'tfunc',
     &             'agwfunc', 'bgwfunc', 'agdefac', 'bgdefac',
     &             'stemp(C)', 'snow', 'snlq', 'srad'
85    format(a10,',',a8,',',7(a12,','),a7,',',a7,',',a12)

      open(unit=90,file='summary.csv')
      write(90,95) 'time','dayofyr','tmax', 'tmin', 'ppt', 'N2Oflux', 
     &             'NOflux', 'CH4', 'NIT', 'CO2resp'
95    format(a10,',',1x,a8,',',1x,3(a8,',',1x),5(a12,',',1x))

      open(unit=123,file='cbalance.csv')
      write(123,96) 'time','dayofyr','balance', 'Cstate1', 'Cstate2', 
     &             'Cgains', 'Closses', 'omadC', 'hrespCO2', 
     &             'co2losUV', 'strm5u', 'strm5l', 'CH4'
96    format(a10,',',1x,a8,',',1x,11(a10,',',1x))

      open(unit=124,file='nbalance.csv')
      write(124,97) 'time','dayofyr','balance', 'Nstate1', 'Nstate2', 
     &              'Ngains', 'Nlosses', 'atmosNdep', 'nonSymSoilNfix',
     &              'omadN', 'fertN', 'orgNlch','inorgNlch',   
     &              'NOflux', 'Nn2oflux', 'Dn2oflux', 'Dn2flux', 
     &              'NOabsorp_g', 'NOabsorp_t' 
97    format(a10,',',1x,a8,',',1x,17(a10,',',1x))

      open(unit=125,file='decomp.csv')
      write(125,98) 'cntstep', 'time', 'dayofyr', 'month', 'timestep', 
     &              'strucc(1)', 'metabc(1)', 'strucc(2)', 'metabc(2)',
     &              'som1c(1)', 'som2c(1)', 'som1c(2)',' som2c(2)',   
     &              'som3c', 'struce(1_1)', 'metabe(1_1)','struce(2_1)', 
     &              'metabe(2_1)', 'som1e(1_1)', 'som2e(1_1)',
     &              'som1e(2_1)', 'som2e(2_1)', 'som3e(1)', 'tnetmn(1)', 
     &              'agdefac', 'bgdefac', 'somsc', 'somtc', 'Rh'
98    format(a10, ',', 1x, 28(a11,',',1x))

      open(unit=128,file='diagnostics.txt')

c ... Calculate available water holding capacity in top 30 cm of soil
c ... as this is the zone where the majority of plant roots occur
      awhc  = 0.0
      adepsum = 0.0
      do 12 lyr = 1, nlayer
        adepsum = adepsum + adep(lyr)
        if (adepsum .le. 30.0) then
          awhc = awhc + (afiel(lyr) - awilt(lyr)) * adep(lyr)
        endif
12    continue

c ... Calculate total water holding capacity in top 30 cm of soil
c ... as this is the zone where the majority of plant roots occur
      twhc = 0.0
      adepsum = 0.0
      do 15 lyr = 1, nlayer
        adepsum = adepsum + adep(lyr)
        if (adepsum .le. 30.0) then
          twhc = twhc + (afiel(lyr) * adep(lyr))
        endif
15    continue

c ... Compute ORGLCH for use in SOMDEC.
      orglch = omlech(1) + omlech(2) * sand

c ... Intercept for the texture equation of secondary P depends upon
c ... pH input.
      if (ph .le. phesp(1)) then
        texesp(2) = phesp(2)
      else if (ph .ge. phesp(3)) then
        texesp(2) = phesp(4)
      else
        dely = phesp(4) - phesp(2)
        delx = phesp(3) - phesp(1)
        xslope = dely / delx
        yint = phesp(2) - (xslope*phesp(1))
        texesp(2) = (xslope*ph) + yint
      endif

      if (micosm .eq. 0) then
c ..... Preset array which will contain monthly values of defac
        do 20 ii = 1, MONTHS
          agdefacm(ii) = -1.
          bgdefacm(ii) = -1.
20      continue
        aagdefac = 0.
        abgdefac = 0.
        agdefac = 0.
        bgdefac = 0.
      else
        call message(' ')
        call message('Microcosms are not implemented in this version')
        call message('of Daily Century')
        STOP
      endif

c ... Effect of soil texture on the microbe decomposition rate
      eftext = peftxa+peftxb*sand

c ... Compute parameters which control decomposition of som1
c ... p1co2 must be computed for surface and soil.   vek  08/91
c ... Note that p1co2b(1) must equal 0 because there is no
c ... soil texture effect on the surface.
      p1co2(SRFC) = p1co2a(SRFC)
      p1co2(SOIL) = p1co2a(SOIL)+p1co2b(SOIL)*sand

c ... Decomposition of som1 to som3 is a function of clay content
c ... vek june90
      fps1s3 = ps1s3(1) + ps1s3(2) * clay
      fps2s3 = ps2s3(1) + ps2s3(2) * clay

      if (texepp(1) .eq. 1.0) then
c ..... Calculate pparmn(2)
c ..... Include effect of texture; weathering factor should be per year
c ..... Note that this code changes the value of a 'fixed' parameter
c ..... (pparmn(2))
        textur = clay + silt
        pparmn(2) = 12.0 * catanf(textur, texepp(2), texepp(3),
     &                            texepp(4), texepp(5))
      endif

      if (texesp(1) .eq. 1.0) then
c ..... Calculate psecmn(2)
c ..... Include effect of texture
c ..... Note that this code changes the value of a 'fixed' parameter
c ..... (psecmn(2))
        psecmn(2) = 12.0 * (texesp(2) + texesp(3) * sand)
      endif

c ... Compute VLOSSG as a function of soil texture based on clay content
c ... vlossg_m is the VLOSSG parameter value as read from the fix.100 file,
      if (clay .lt. 0.10) then
c        vlossg = 0.015
        vlossg = 0.03
c      else if (clay .gt. 0.40) then
      else if (clay .gt. 0.30) then
c        vlossg = 0.003
        vlossg = 0.01
      else
c        vlossg = line(clay, 0.10, 0.015, 0.40, 0.003)
        vlossg = line(clay, 0.10, 0.03, 0.30, 0.01)
      endif
      vlossg = vlossg * vlossg_m

c ... Save initial values for printing or plotting
      call savarp

c ... Clear the flow stack.
      call floclr

      return
      end
