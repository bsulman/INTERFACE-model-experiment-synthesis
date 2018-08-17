
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

c***********************************************************************
c**
c**  FILE:     getwth.f
c**
c**  PURPOSE:  Retrieve a day's worth of weather for the weather file.
c**            Compute monthly average temperature using a
c**            circular arary over a 30 day period.
c**
c**  This routine was developed for the RAMS / Daily Century linkage
c**
c**  Melannie D. Hartman
c**  12/5/96
c**
c**  Add more robust checking for valid weather data values.
c**  CAK - 04/09/01
c**
c**  INPUTS:
c**     curday       - the day of the year to read weather from the file
c**     month        - current month of the year (1..12)
c**     tmn2m        - average minimum air temperature for the month
c**                    as read from <site>.100 file (deg C - 2m)
c**     tmx2m        - average maximum air temperature for the month
c**                    as read from <site>.100 file (deg C - 2m)
c**  OUTPUTS:
c**     (From weather file):
c**     tempmax - maximum air temperature for the day (deg C)
c**     tempmin - minimum air temperature for the day (deg C)
c**     avgtemp - average air temperature for the day (deg C)
c**     ppt     - precipitation for the day (cm)
c**     solrad  - total incoming shortwave radiation (langleys/day)
c**     srad    - total incoming shortwave radiation (W/m^2/day)
c**     rhumid  - average relative humidity for the day (% 1..100)
c**     vpd     - vapor pressure deficit (kPa/day)
c**     windsp  - average daily windspeed at 2 meters (mph)
c** 
c**  Called by:  simsom.f
c**
c**  Calls:  none
c**
c***********************************************************************

      subroutine getwth(curday, month, tmn2m, tmx2m, 
     &                  tempmax, tempmin, avgtemp, ppt,
     &                  solrad, rhumid, windsp, srad, vpd)

      implicit none
      include 'dconst.inc'
      include 'jday.inc'

c ... Formal parameters

      integer curday, month
      double precision tempmax(NDAY+1),tempmin(NDAY+1),avgtemp(NDAY+1),
     &  ppt(NDAY+1),solrad(NDAY+1),rhumid(NDAY+1),windsp(NDAY+1)
      double precision tmn2m(NMONTH), tmx2m(NMONTH)
      double precision srad(NDAY+1),vpd(NDAY+1)

c ... Local Variables

      integer ndy, nyr, njday, nmth, imo
      double precision dailyPrecip, maxTemp, minTemp
      double precision tdew
      logical debug

      debug = .false.

c ... Return here if we have reached the end of the weather data file
10    continue

c ... Use extra weather drivers for PET calculations
      if (usexdrvrs .eq. 1) then
        read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(curday),
     &                   tempmin(curday),ppt(curday),solrad(curday),
     &                   rhumid(curday), windsp(curday)
        if ((solrad(curday) .le. -99.0) .or. 
     &      (rhumid(curday) .le. -99.0) .or. 
     &      (windsp(curday) .le. -99.0)) then
          write(*,*) 'Invalid weather data, day ', curday
          STOP
        endif
        srad(curday) = -999.0
        vpd(curday) = -999.0
c ... Use extra weather drivers for photosynthesis calculations
      elseif (usexdrvrs .eq. 2) then
        read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(curday),
     &                   tempmin(curday),ppt(curday), srad(curday),
     &                   vpd(curday)
        solrad(curday) = -999.0
        rhumid(curday) = -999.0
        windsp(curday) = -999.0
c ... Use extra weather drivers for both PET and photosynthesis
c ... calculations
      elseif (usexdrvrs .eq. 3) then
        read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(curday),
     &                   tempmin(curday),ppt(curday),solrad(curday),
     &                   rhumid(curday), windsp(curday), srad(curday),
     &                   vpd(curday)
        if ((solrad(curday) .le. -99.0) .or. 
     &      (rhumid(curday) .le. -99.0) .or. 
     &      (windsp(curday) .le. -99.0)) then
          write(*,*) 'Invalid weather data, day ', curday
          STOP
        endif
c ... No extra weather drivers used
      else
        read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(curday),
     &                   tempmin(curday),ppt(curday)
        solrad(curday) = -999.0
        rhumid(curday) = -999.0
        windsp(curday) = -999.0
        srad(curday) = -999.0
        vpd(curday) = -999.0
      endif

c ... Checks for valid weather data
      if (tempmax(curday) .le. -99.0) then
        if (debug) then
          write(*,*) 'Warning: missing maximum temperature data, ',
     &               'day ', curday, ' year ', nyr
        endif
        tempmax(curday) = tmx2m(month)
      endif
      if (tempmin(curday) .le. -99.0) then
        if (debug) then
          write(*,*) 'Warning: missing minimum temperature data, ',
     &               'day ', curday, ' year ', nyr
        endif
        tempmin(curday) = tmn2m(month)
      endif
      if (ppt(curday) .le. -99.0) then
        if (debug) then
          write(*,*) 'Warning:  missing precipitation data, day ',
     &               curday, ' year ', nyr
        endif
        ppt(curday) = 0
      endif
      if (tempmax(curday) .lt. tempmin(curday)) then
        write(*,*) 'Warning:  invalid weather data, day ', curday,
     &             ' year ', nyr, ', tmax < tmin'
        tempmax(curday) = tmx2m(month)
        tempmin(curday) = tmn2m(month)
      endif

      goto 30

c ... If necessary, start reading the weather file again from the beginning      
20    rewind(9)
      goto 10

30    continue

      if (njday .ne. curday) then
        write(*,*) 'Expect day ', curday, ' got day ', njday
        write(*,*) 'in weather file.'
        STOP
      endif

      if (nmth .ne. month) then
        write(*,*) 'Expect month ', month, ' got month ', nmth
        write(*,*) 'in weather file.'
        STOP
      endif

c ... Check for leap year
      if (curday .eq. 1) then 
        if ((mod(nyr,400) .eq. 0)) then
          leapyr = .TRUE.
        else if ((mod(nyr,100) .eq. 0)) then
          leapyr = .FALSE.
        else if ((mod(nyr,4) .eq. 0)) then
          leapyr = .TRUE.
        else
          leapyr = .FALSE.
        endif
        if (leapyr) then
          dysimo(2) = idysimo(2)+1
          do 40 imo = 2, 12
            lstdy(imo) = ilstdy(imo)+1
            if (imo .gt. 2) frstdy(imo) = ifrstdy(imo)+1
40        continue
        endif
      endif

      avgtemp(curday) = (tempmax(curday) + tempmin(curday)) / 2.0

c ... If the solar radiation value has not been read from the weather data
c ... file calculate the solar radiation for the day using Peter Thornton's
c ... code, cak - 06/18/2009
      if (srad(curday) .le. -99.0) then
        minTemp = tempmin(curday)
        maxTemp = tempmax(curday)
        dailyPrecip = ppt(curday)
        tdew = 0.0
        call calc_srad(maxTemp, minTemp, dailyPrecip, curday, tdew,
     &                 srad(curday))
c ..... Adjust the total daily radiation value returned from Peter
c ..... Thornton's code to reflect the the transmission coefficient
c ..... and cloud cover at the site.
        srad(curday) = srad(curday) * sradadj(month)
      endif

      return
      end
