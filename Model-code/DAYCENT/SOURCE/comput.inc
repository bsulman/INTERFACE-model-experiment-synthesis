
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Variables used internally which must be shared among routines.
c ...
c ... Change p1co2 to p1co2(2) to represent surface and soil layers.
c ...
c ... baseNdep - base amount of atmospheric N deposition for the year
c ...            based on average annual precipitation only
c ... fps1s3 - the fraction of som1 decomposing to som3
c ... fps2s3 - the fraction of som2 decomposing to som3
c ... orglch - the fraction of organics that leach from soil som1 when
c ...          there is sufficient water flow.  Computed in predec.
c ... ratnew(iel,1) - the C/E ratio for new material created when a
c ...                 lignin component decomposes to som1.
c ... ratnew(iel,2) - the C/E ratio for new material created when a
c ...                 lignin component decomposes to som2.


      common/comput/agdefacm(12),bgdefacm(12),
     &  eftext,fps1s3,fps2s3,
     &  orglch,p1co2(2),
     &  ratnew(3,2)

      double precision agdefacm,bgdefacm,
     &  eftext,fps1s3,fps2s3,orglch,p1co2,
     &  ratnew

      save /comput/
