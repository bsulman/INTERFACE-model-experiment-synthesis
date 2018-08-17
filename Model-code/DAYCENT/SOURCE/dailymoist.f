
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine dailymoist(curday, sradKJ, soilsrad, 
     &                      grass_lai, tree_lai, isdecid, isagri,
     &                      tfunc, agwfunc, bgwfunc, petdly, rpeff)

      implicit none
      include 'cflows.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'jday.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'timvar.inc'
      include 't0par.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

c ... FORMAL PARAMETERS
      integer curday
      integer isdecid, isagri
      double precision sradKJ, soilsrad
      double precision grass_lai, tree_lai
      double precision tfunc
      double precision agwfunc, bgwfunc
      double precision petdly, rpeff

c ... Subroutine dailymoist
c ... 
c ... This routine is called once a day: it calls 
c ... the water budget routine (watrflow), 
c ... the photo decomposition routine (photodecomp), 
c ... the soil organic mattter decomposition routine (decomp),   
c ... and the trace_gas_model.
c ...
c ... FORMAL PARAMETERS
c ... Input (used in this subroutine)
c ...   curday 		- day of the year (1..366)
c ...   isdecid 	- 1 if the system is a deciduous forest, 0 otherwise
c ...   isagri 		- 1 if the system has ever been cultivated for agriculture, 0 otherwise
c ...   sradKJ    	- amount of solar radiation (KJ/m^2/day)
c ...   soilsrad	- amount of solar radiation below the plant canopy (KJ/m^2/day)
c ...   grass_lai	- LAI of crop/grass (m^2/m^2)
c ...   tree_lai	- LAI of tree canopy (m^2/m^2)
c ... Output (for output files)
c ...   tfunc		- temperature effect on decomposition (0.0 - 1.0)
c ...   agwfunc		- above ground moisture effect on decomposition (0.0 - 1.0)
c ...   bgwfunc		- below ground moisture effect on decomposition (0.0 - 1.0)
c ...   petdly		- potential evapotranspiration that has occurred in the current day (cm H2O)
c ... 
c ... LOCAL VARIABLES
c ... Temperature and Moisture
c ...   avgwfps 	- average WFPS in top 15 cm of soil (0-1)
c ...   fertN	        - amount of N fertilizer added (gN/m^2)
c ...   omadC	        - amount of C from organic matter addition (gC/m^2)
c ...   omadN	        - amount of N from organic matter addition (gN/m^2)
c ...   rprpet		- ratio of rainfall (including melt) to PET 
c ...   amovdly(10)	- total net flux of water thru the bottom of each Century soil layer each day (cm H2O)
c ...   wfluxout[] 	- total net flux of water thru the bottom of each DayCent soil layer each day (cm H2O)
c ...              	  (positive is downward, negative is upward)
c ... Surface/Soil CO2 losses
c ...   CO2losUV	- CO2 loss from photodegradation of standing dead biomass and 
c ...                     surface litter in the current day (gC/m^2)
c ...   CO2resp		- amount of SOIL microbial respiration in the current day (gC/m^2)
c ...   newCO2 		- modified value of CO2resp passed to the trace gas submodel (gC/m^2)
c ...   hrespCO2	- amount of microbial respiration (surface, soil, dead wood) in the current day (gC/m^2)
c ...   hresp1		- sum of ALL microbial respiration accumulators before call to decomp (gC/^2)
c ...   hresp2		- sum of ALL microbial respiration accumulators after call to decomp (gC/^2)
c ...   soilresp1	- sum of SOIL microbial respiration accumulators before call to decomp (gC/^2)
c ...   soilresp2	- sum of SOIL microbial respiration accumulators after call to decomp (gC/^2)
c ...   inorgNlch	- inorganic N leaching that has occurred in the current day (gN/m2)
c ...   orgNlch		- organic N leaching that has occurred in the current day (gN/m2)
c ...   newminrl	- mineralization that has occurred in the current day (gN/m2)
c ... Soil Leaching
c ...   fsol 		- fraction of soluble mineral in soil mineral layers (0.0 - 1.0)
c ...   frlech(MAXIEL)	- inorganic N, P, and S leaching fractions (0.0 - 1.0), used in leachdly.
c ...   texeff		- effect of soil texture (sand) in the leaching of inorganic N, P, and S
c ... Trace Gas Model Outputs
c ...   CH4		- methane consumption in the current day (gC/m^2)
c ...   Dn2oflux	- denitrification N2O flux in the current day (gN/m^2)
c ...   Dn2flux		- denitrification N2 flux in the current day (gN/m^2)
c ...   Nn2oflux	- nitrification N2O flux in the current day (gN/m^2)
c ...   NOflux		- total soil NOx flux in the current day (gN/m^2)
c ...   Note: the following 2 fluxes should increase grass and tree N storage, but for CLM linkage
c ...      esrsnk(N) will be increased instead.
c ...   NOabsorp_grass	- total NOx flux absorbed by grass in the current day (gN/m^2)
c ...   NOabsorp_tree	- total NOx flux absorbed by tree canopy in the current day (gN/m^2)
c ...   nit_amt		- amount of nitrification in the current day (gN/m^2)
c ...   dN2lyr() 	- denitrification N2 flux in the current day by soil layer (gN/m^2)
c ...   dN2Olyr()	- denitrification N2O flux in the current day by soil layer (gN/m^2)
c ... N FIXATION
c ...   atmosNdep 	- atmospheric N depostion in the current day (gN/m^2)
c ...   nonSymSoilNfix	- non-symbiotic soil N fixation in the current day (gN/m^2)
c ... C and N balance
c ...   totCstate	- total of C state variables (gC/m^2)
c ...   totNstate	- total of N state variables (gC/m^2)
c ...
c ... OTHER VARIABLES defined in include files
c ...
c ... plot1.inc
c ...   aminrl(iel)	- available mineral in the top Century soil layer (gE/m^2)
c ...   asmos(lyr)	- total water in each Century soil layer (cm_
c ...   avh2o(3)	- available water in the soil(cm)
c ...   minerl(lyr,iel) - mineral soil for Century layers (assumes nitrate + ammonium) (gE/m^2) 
c ...   rwcf(lyr)	- relative water content fraction of Century layers (0.0 - 1.0)
c ...   snow		- snow water equivalent of the snowpack (cm H2O)
c ...   stream(*)       - water, inorganic stream flow, organic stream flow 
c ...
c ... plot2.inc
c ...   esrsnk(iel)	- elemental source/sink (gE/m^2)
c ...
c ... npool.inc
c ...   ammonium	- ammonium-nitrogen pool, no layer structure (gN/m^2)
c ...   nitrate[]	- nitrate-nitrogen pool for DayCent soil layers (gN/m^2)
c ...   texture		- COARSE, MEDIUM, or FINE (see n2o_model.h)
c ...
c ...  dovars.inc
c ...    dofert		- true if it is time for a fertilization event
c ...    fertday	- scheduled day of year for fertilization event
c ...    doomad		- true if it is time for a organic matter addition event
c ...    omadday	- scheduled day of year for organic matter addition event
c ...
c ...  fertil.inc
c ...    nreduce
c ...
c ...  param.inc
c ...    afiel(1)	- field capacity in the top Century soil layer (volumetric fration 0.0 - 1.0)
c ...    basef
c ...    bulkd 		- avgerage soil bulk density for the soil profile (g/cm^3)
c ...    stormf
c ...    drain
c ...    nelem
c ...    pslsrb
c ...    sorpmx
c ...
c ...  parfx.inc
c ...    aneref
c ...    ntspm		- number of time steps per day (not month) for the decomp model
c ...    idef
c ...    fleach(*)	- leaching parameters
c ...    minlch
c ...    teff(*)	- 
c ...
c ...  seq.inc
c ...    decsys		- 1 if no wood is present to decompose, 2 otherwise
c ...
c ...  site.inc
c ...    sand, silt, clay
c ...
c ...  t0par.inc
c ...    dt 		- time step of a month (fraction of a year)
c ...  timvar.inc
c ...    decodt		- time step for decomposition model (fraction of a year)
c ...    month
c ...
c ... wth.inc
c ...   maxt 		- long term average maximum monthly temperature (C) used in the trace gas model
c ...
c ... wthdaily.inc
c ...   avgtemp    
c ...   tempmin
c ...   tempmax
c ...
c ... zztime.inc
c ...   time
c ...   time2


c ... LOCAL VARIABLES
      integer iel, kts, ntspd, ii
      integer iyr, imo, idy, ikts
      double precision cdi, amtN
      double precision fsol 
      double precision texeff
      double precision frlech(MAXIEL)
      double precision fertN, omadN, omadC
      double precision inorgNlch
      double precision orgNlch
      double precision newminrl
      double precision amovdly(10)
      double precision rprpet
      double precision wfluxout(SWMAXLYR)
      double precision hresp1, hresp2
      double precision soilresp1, soilresp2
      double precision Nn2oflux, Dn2oflux, Dn2flux, NOflux
      double precision NOabsorp_grass, NOabsorp_tree
      double precision avgwfps
      double precision CH4, nit_amt
      double precision CO2losUV
      double precision CO2resp
      double precision newCO2
      double precision hrespCO2
      double precision dN2lyr(SWMAXLYR), dN2Olyr(SWMAXLYR)
      double precision atmosNdep, nonSymSoilNfix
      double precision NgasFlux
      double precision totCstate, totNstate
      double precision t1, t2, timeincrmt
      character subname*10
      integer experiment

c ... FUNCTION DECLARATIONS
      double precision anerob, fsfunc
      external anerob, fsfunc

      subname = 'dailymoist'

c ... experiment = 1, 2, 3 for litter experiments.  Set to 0 for normal execution.
c     0 = normal DayCent execution
c         compute atmospheric N dep and non-symbiotic N fixation
c ... For all non-zero experiments, set pH = 7, soilsrad = 0
c     No atmospheric N dep or non-sybiotic N fixation
c     1 = anerob = 1.0
c     2 = read cdi file instead of calculating defac
c         Reset minerl(*,*) and aminrl(*) to values specified in the site.100 file
c         Don't run tracegas model
c         anerob = 1.0
c     3 = 

      experiment = 2
      
      totCstate = 0
      omadC = 0
      CO2resp = 0
      CO2losUV = 0
      CH4 = 0
      totNstate = 0
      fertN = 0
      omadN = 0
      orgNlch = 0
      inorgNlch = 0
      Nn2oflux = 0
      Dn2oflux = 0 
      Dn2flux = 0
      NOflux = 0
      NOabsorp_grass = 0
      NOabsorp_tree = 0
      atmosNdep = 0
      nonSymSoilNfix = 0

      call cbalance(time, curday, omadC, hrespCO2, CO2losUV, 
     &              CH4, totCstate)
      call nbalance(time, curday, omadN, fertN, 
     &              orgNlch, inorgNlch, 
     &              NOflux, Nn2oflux, Dn2oflux, Dn2flux, 
     &              NOabsorp_grass, NOabsorp_tree, 
     &              atmosNdep, nonSymSoilNfix, totNstate)

      decodt = dt/dble((dysimo(month)*ntspm))
c ... daily timestep (t1) has 24 * 60 minutes, the new timestep (t2) has 30 minutes
      t1 = 24 * 60
      t2 = 30
      ntspd = ntspm * 48
c ... decomposition time increment as a fraction of a year (48 decomposition timesteps per day)
      timeincrmt = 1.0 / (365.0 * 48.0)

      if (experiment .ne. 0) then
        soilsrad = 0
        ph = 7.0
      endif

      newminrl = 0.0

c ... Initialize output variables that are tracking the daily carbon
c ... flows due to decomposition (cflows.inc).
      metc1tosom11  = 0.0
      metc2tosom12  = 0.0
      struc1tosom11 = 0.0
      struc1tosom21 = 0.0
      struc2tosom12 = 0.0
      struc2tosom22 = 0.0
      som11tosom21  = 0.0
      som12tosom22  = 0.0
      som12tosom3   = 0.0
      som21tosom11  = 0.0
      som21tosom22  = 0.0
      som22tosom12  = 0.0
      som22tosom3   = 0.0
      som3tosom12   = 0.0
      wood1tosom11  = 0.0
      wood1tosom21  = 0.0
      wood2tosom11  = 0.0
      wood2tosom21  = 0.0
      wood3tosom12  = 0.0
      wood3tosom22  = 0.0

c ... Atmospheric N deposition and non-symbiotic soil N fixation 
      if (experiment .eq. 0) then
        call nfixday(ppt(curday), atmosNdep, nonSymSoilNfix)
      endif

c ... Organic matter addition
      if (doomad .and. (omadday .eq. curday)) then
        call addomad(omadC, omadN)
        write(*,*) 'DOOMAD: ', time, curday, omadC, omadN
        doomad = .FALSE.
      endif

c.... Add plant litter
c     call addlit()

c ... Compute frlech, the leaching fraction.  
      texeff = fleach(1) + fleach(2) * sand
      do 50 iel = 1, nelem
        if (iel .eq. P) then
          fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
        else
          fsol = 1.0
        endif
        frlech(iel) = texeff * fleach(iel+2) * fsol
50    continue

c ... Fertilization option
      if (dofert .and. curday .eq. fertday) then
        call addfert(nelem, month, fertN)
      endif

      call watrflow(nlayer, 
     &              avgtemp(curday), tempmin(curday), tempmax(curday),  
     &              ppt(curday), petdly, 
     &              wfluxout, amovdly, rwcf, asmos, avh2o, rprpet)

c ... Calculate the effect impact of anerobic conditions on decomposition
      anerb = anerob(aneref,drain,rprpet,petdly,experiment)

c ... Combined effects of temperature and moisture on decomposition
      if (experiment .ne. 2) then
        call calcdefac(texture, tfunc, agwfunc, bgwfunc, 
     &                 agdefac, bgdefac, avgwfps, teff, rprpet,  
     &                 idef, ppt(curday), snow, experiment)
      endif

      if (experiment .eq. 2) then
c...    Reset minerl(*,*) and aminrl(*) to values specified in the site.100 file
        do 19 iel = 1, 3
          do 17 ii = 1, nlayer
            amtN = iminerl(ii,N)-minerl(ii,N)
            minerl(ii, iel) = iminerl(ii, iel)
            if (iel .eq. N) then
              call update_npool(ii, amtN, frac_nh4_fert,
     &            frac_no3_fert, ammonium, nitrate, subname)
            endif
17        continue
          if (iel .eq. P) then
            fsol = fsfunc(minerl(1,P), pslsrb, sorpmx)
          else
            fsol = 1.0
          endif
          aminrl(iel) = minerl(1,iel) * fsol
19      continue
      endif

c ... *********************************************************************
c ... Decomposition Submodel
c ... Call decomp routines ntspd times per day

c ... Change CO2 respiration value passed to the trace gas submodel
c ... so that we are passing only soil respiration, cak - 10/22/2007
      soilresp1 = mt2c2(UNLABL) + mt2c2(LABELD) +
     &            st2c2(UNLABL) + st2c2(LABELD) +
     &            s12c2(UNLABL) + s12c2(LABELD) +
     &            s22c2(UNLABL) + s22c2(LABELD) +
     &            s3c2(UNLABL)  + s3c2(LABELD)
      hresp1 = soilresp1 + 
     &            mt1c2(UNLABL) + mt1c2(LABELD) + 
     &            st1c2(UNLABL) + st1c2(LABELD) +
     &            s11c2(UNLABL) + s11c2(LABELD) +
     &            s21c2(UNLABL) + s21c2(LABELD) +
     &            wd1c2(UNLABL) + wd1c2(LABELD) +
     &            wd2c2(UNLABL) + wd2c2(LABELD) +
     &            wd3c2(UNLABL) + wd3c2(LABELD) 

c ... Photodecomposition of standing dead and litter
c     call photodecomp(sradKJ, soilsrad, CO2losUV)

      orgNlch = stream(6)

      do 40 kts = 1, ntspd

        if (experiment .eq. 2) then
c ........ Read defac from the cdi file, rewind the file at EOF
181        read(18, *, end=182) iyr, imo, idy, ikts, cdi
           goto 183

182        rewind(18)
           goto 181

183        continue
           agdefac = cdi
           bgdefac = cdi
           if (month .ne. imo) then
             write(*,*) 'Warning: month = ', month, ' imo = ', imo
           endif
           if (kts .ne. ikts) then
             write(*,*) 'Warning: kts = ', kts, ' ikts = ', ikts
           endif
        endif

c ... Organic matter decomposition
        call decomp(decodt,t1,t2,decsys,amovdly,newminrl,soilsrad,rpeff)

        if (nelem .ge. P) then
          call pschem(decodt, t1, t2)
        endif

c ..... Update decomposition and nitrogen fixation flows.
        call flowup(time)
        call sumcar

c ..... aminrl contains the average amount of N, P, and S available
c ..... in the top layer for the time period covered by dt/ntspd.  
c ..... minerl contains the current value of mineral N, P, and S by layer.

        do 30 iel = 1, nelem
          if (iel .eq. P) then
            fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
          else
            fsol = 1.0
          endif
          aminrl(iel) = (aminrl(iel) + minerl(SRFC,iel)*fsol)/2.0
30      continue
c       write(*,*) cntstep, ': minerl(1,N) = ', minerl(1,N), 
c    &      'aminrl(N) = ', aminrl(N)

      cntstep = cntstep + 1
      time2 = time2 + timeincrmt

40    continue

c ... newCO2 - soil CO2 respiration value passed to the trace gas submodel
      soilresp2 = mt2c2(UNLABL) + mt2c2(LABELD) +
     &            st2c2(UNLABL) + st2c2(LABELD) +
     &            s12c2(UNLABL) + s12c2(LABELD) +
     &            s22c2(UNLABL) + s22c2(LABELD) +
     &            s3c2(UNLABL)  + s3c2(LABELD)
      hresp2 = soilresp2 + 
     &            mt1c2(UNLABL) + mt1c2(LABELD) + 
     &            st1c2(UNLABL) + st1c2(LABELD) +
     &            s11c2(UNLABL) + s11c2(LABELD) +
     &            s21c2(UNLABL) + s21c2(LABELD) +
     &            wd1c2(UNLABL) + wd1c2(LABELD) +
     &            wd2c2(UNLABL) + wd2c2(LABELD) +
     &            wd3c2(UNLABL) + wd3c2(LABELD) 
      newCO2 = soilresp2 - soilresp1      

      write(125,98) cntstep, time2, curday, month, kts, 
     &              strucc(1), metabc(1), strucc(2), metabc(2),
     &              som1c(1), som2c(1), som1c(2), som2c(2),   
     &              som3c, struce(1,1), metabe(1,1), struce(2,1), 
     &              metabe(2,1), som1e(1,1), som2e(1,1),
     &              som1e(2,1), som2e(2,1), som3e(1), tnetmn(1), 
     &              agdefac, bgdefac, somsc, somtc, newCO2
98    format(i10,',',f12.8,',',1x,3(i4,',',1x),24(f11.2,',',1x))

c ... dshresp is daily soil heterotrophic respiration, needed to calculate 
c ... priming effect for the next day (see rpeff calculation in simsom.f).
      dshresp = newCO2
      hrespCO2 = hresp2 - hresp1
      if (newCO2 .le. 0.0) then
        newCO2 = 0.000001
      endif
      CO2resp = newCO2
      if (avgwfps .gt. 0.60 .and. bgwfunc .gt. 0.0) then
        newCO2 = newCO2 / bgwfunc
      endif

c ... *********************************************************************
c ... Trace Gas Model
   
      if (experiment .ne. 2) then
      call trace_gas_model(curday, time, newminrl, ammonium, nitrate, 
     &                     texture, sand, silt, clay, afiel(1), bulkd, 
     &                     maxt, ppt(curday), snow, avgwfps, 
     &                     stormf, basef, frlech, stream, inorgNlch,
     &                     minlch, wfluxout, newCO2, NOflux,
     &                     Nn2oflux, Dn2oflux, Dn2flux, CH4, isdecid,
     &                     isagri, grass_lai, tree_lai, NOabsorp_grass,
     &                     NOabsorp_tree, nit_amt, nreduce, 
     &                     dN2lyr, dN2Olyr)
      endif

      NgasFlux = NOflux + Nn2oflux + Dn2oflux + Dn2flux + 
     &           NOabsorp_grass + NOabsorp_tree
      esrsnk(N) = esrsnk(N) + NgasFlux

c ... Update state variables and accumulators and sum carbon isotopes
      call flowup(time)
      call sumcar

c ... Now check for N balance and rebalance nh4 and no3 pools with minerl N
      minerl(1,N) = minerl(1,N) - Nn2oflux - Dn2oflux - Dn2flux -
     &                NOflux 

      call bal_npool(nlayer, minerl, ammonium, nitrate, inorgNlch)

c ... *********************************************************************

c ... Organic N leaching for the day
      orgNlch = stream(6) - orgNlch

      call cbalance(time, curday, omadC, hrespCO2, CO2losUV, 
     &              CH4, totCstate)
      call nbalance(time, curday, omadN, fertN, 
     &              orgNlch, inorgNlch, 
     &              NOflux, Nn2oflux, Dn2oflux, Dn2flux, 
     &              NOabsorp_grass, NOabsorp_tree, 
     &              atmosNdep, nonSymSoilNfix, totNstate)


      return
      end
