
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      common/param/afiel(10),amov(10),autoresp1(2),autoresp2(2),
     &    awilt(10),basef,bioabsorp,bulkd,
     &    co2ipr(2),co2ice(2,2,3),co2irs(2),co2itr(2),co2sys,co2tm(2),
     &    crpindx,crpcmn,crpcmx,crpmnmul,crpmxmul,
     &    drain,dshresp,epnfa(2),epnfs(2),falprc,fracro,
     &    hpttr(12),htran(12),ivauto,labtyp,labyr,
     &    ckmrspmx(3),cmrspnpp(6),fkmrspmx(6),fmrsplai(6),
     &    cgresp(3),fgresp(6),
     &    maxphoto,maxtmp(12),mctemp,micosm,mintmp(12),
     &    nelem,newautoresp(2),Ninput,nlayer,nlaypg,no3pref(2),
     &    Nstart,OMADinput,OMADstart,
     &    ph,phsys,phtm,
     &    ppdf(4,2),prcskw(12),prcstd(12),prdx(2),
     &    ps2mrsp(2),precip(12),precro,psloss,pslsrb,
     &    rcelit(2,3),rces1(2,3),rces2(2,3),rces3(3),remwsd,rock,
     &    satmos(2),satmt,sirri,snfxmx(2),sorpmx,stamt,stormf,
     &    strm5l,strm5u,ststart,stsys,swflag,trbasl,
     &    claypg,claypg_const,tlaypg,cmix,tmix

      integer falprc,ivauto,labtyp,labyr,micosm,nelem,
     &        Ninput,nlayer,nlaypg,Nstart,OMADinput,OMADstart,phsys,
     &        phtm,swflag,claypg,claypg_const,tlaypg,crpindx

      double precision    afiel,amov,autoresp1,autoresp2,
     &        awilt,basef,bioabsorp,bulkd,
     &        co2ipr,co2ice,co2irs,co2itr,co2sys,co2tm,
     &        crpcmn,crpcmx,crpmnmul,crpmxmul,
     &        drain,dshresp,epnfa,epnfs,fracro,
     &        hpttr,htran,ckmrspmx,cmrspnpp,fkmrspmx,fmrsplai,
     &        cgresp,fgresp,
     &        maxphoto,maxtmp,mctemp,mintmp,
     &        newautoresp,no3pref,
     &        ph,
     &        ppdf,prcskw,prcstd,prdx,
     &        ps2mrsp,precip,precro,psloss,pslsrb,
     &        rcelit,rces1,rces2,rces3,remwsd,rock,
     &        satmos,satmt,sirri,snfxmx,sorpmx,stamt,stormf,
     &        strm5l,strm5u,ststart,stsys,trbasl,
     &        cmix,tmix

      save /param/
