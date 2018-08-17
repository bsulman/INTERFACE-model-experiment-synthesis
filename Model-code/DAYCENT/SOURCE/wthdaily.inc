      include 'dconst.inc'

      common /wthdaily/ srad(NDAY+1),
     &                  vpd(NDAY+1),
     &                  avgtemp(NDAY+1),
     &                  tempmax(NDAY+1),
     &                  tempmin(NDAY+1),
     &                  ppt(NDAY+1),
     &                  solrad(NDAY+1),
     &                  rhumid(NDAY+1),
     &                  windsp(NDAY+1)

      double precision srad, vpd

      double precision avgtemp, tempmin, tempmax, ppt
      double precision solrad, rhumid, windsp

      save /wthdaily/
