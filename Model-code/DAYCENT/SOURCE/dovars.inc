
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Scheduling option variables
c ... docult      time for a cultivation event
c ... doerod      time for an erosion event
c ... dofert      time for a fertilization event
c ... dofire      time for a fire event
c ...             (3) - grass/crop,forest,savannah
c ... doflst      time for a forest last month growth event
c ... dofone      time for a forest first month growth event
c ... dofrst      time for a crop first month growth event
c ... dograz      time for a grazing event
c ... dohrvt      time for a harvest event
c ... doirri      time for an irrigation event
c ... dolast      time for a crop last month growth event
c ... doomad      time for a organic matter addition event
c ... doplnt      time for a planting event
c ... dosene      time for a senescence event
c ... dotrem      time for a forest removal event
c ... cultday     scheduled day of year for cultivation event
c ... erodday     scheduled day of year for erosion event
c ... fertday     scheduled day of year for fertilization event
c ... fireday     scheduled day of year for fire event
c ... flstday     scheduled day of year for forest last month growth event
c ... foneday     scheduled day of year for forest first month growth event
c ... frstday     scheduled day of year for crop first month growth event
c ... grazday     scheduled day of year for grazing event
c ... hrvtday     scheduled day of year for harvest event
c ... irriday     scheduled day of year for irrigation event
c ... lastday     scheduled day of year for crop last growth month event
c ... omadday     scheduled day of year for organic matter addition event
c ... plntday     scheduled day of year for planting event
c ... seneday     scheduled day of year for senescence event
c ... tremday     scheduled day of year for forest removal event
c ... cultcnt     number of days that cultivation effect on decomposition
c ...             has occurred
c ... erodcnt     number of days that the erosion event has occurred
c ... fertcnt     number of days since the fertilization event has occurred
c ... grazcnt     number of days that the grazing event has occurred
c ... irricnt     number of days that the irrigation event has occurred
c ... plntcnt     number of days since the planting event has occurred
c ... senecnt     number of days since the senescence event has occurred
c ... frstschd    set to .true. when a FRST event is scheduled, .false. when
c ...             the FRST event occurs
c ... harvschd    set to .true. when a HARV event is scheduled, .false. when
c ...             the harvest event occurs
c ... plntschd    set to .true. when a PLNT event is scheduled, .false. when
c ...             the planting event occurs
c ... senmschd    set to .true. when a SENM event is scheduled, .false. when
c ...             the senensence event occurs
c ... savefrstday scheduled day of year for crop first month growth event as
c ...             read from the schedule file
c ... saveplntday scheduled day of year for planting event as read from
c ...             schedule file

      common/dovars/docult, doerod, dofert, dofire(3), doflst,
     &              dofone, dofrst, dograz, dohrvt, doirri,
     &              dolast, doomad, doplnt, dosene, dotrem, cultday,
     &              erodday, fertday, fireday, flstday, foneday,
     &              frstday, grazday, hrvtday, irriday, lastday,
     &              omadday, plntday, seneday, tremday, cultcnt,
     &              erodcnt, fertcnt, grazcnt, irricnt, plntcnt,
     &              senecnt, frstschd, harvschd, plntschd, senmschd,
     &              savefrstday, saveplntday
      logical docult, doerod, dofert, dofire, doflst, dofone,
     &        dofrst, dograz, dohrvt, doirri, dolast, doomad,
     &        doplnt, dosene, dotrem, frstschd, harvschd,
     &        plntschd, senmschd
      integer  cultday, erodday, fertday, fireday, flstday, foneday,
     &         frstday, grazday, hrvtday, irriday, lastday, omadday,
     &         plntday, seneday, tremday
      integer  cultcnt, fertcnt, erodcnt, grazcnt, irricnt, plntcnt,
     &         senecnt
      integer  savefrstday, saveplntday

      save /dovars/
