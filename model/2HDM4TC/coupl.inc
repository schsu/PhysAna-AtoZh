

c Define common block containing all coupling constants and masses
c  which are used in the HELAS routines.
c
c  These include masses, widths and real/complex couplings.
c
c  This file can be built automatically from particles.dat and
c  interactions.dat
c
c====================================================================
c
c
c
c     QCD 
c
      double complex       gg(2)
      double precision     g
      common /COUPL_QCD/   g,gg   
c
c    kinematical masses
c
      double precision     wmass, zmass, amass,
     &                     tmass, bmass, lmass,
     &                     cmass, smass, mmass
      common /COUPL_MASS/  wmass, zmass, amass,
     &                     tmass, bmass, lmass,
     &                     cmass, smass, mmass
c
c   widths
c
      double precision     wwidth, zwidth, 
     &                     twidth, lwidth, awidth
      common /COUPL_WIDTH/ wwidth, zwidth, 
     &                     twidth, lwidth, awidth
c
c  couplings in the feynman rules
c
      double complex       gal(2), gad(2), gau(2), gwf(2),
     &                     gzn(2), gzl(2), gzd(2), gzu(2)
      double precision     gw, gwwa, gwwz
      common /COUPL_GAUGE/ gal   , gad   , gau   , gwf   ,
     &                     gzn   , gzl   , gzd   , gzu   ,
     &                     gw, gwwa, gwwz
c
      double complex       gwfc(2),  gwfs(2), gwfm(2)
      common /coupl_ckm/   gwfc,     gwfs   , gwfm	

      double complex   GH1LL(2),GH2LL(2),GH3LL(2),GH1MM(2),GH2MM(2),GH3MM(2),GH1TT(2),GH2TT(2),GH3TT(2),GH1CC(
     +2),GH2CC(2),GH3CC(2),GH1BB(2),GH2BB(2),GH3BB(2),GH1SS(2),GH2SS(2),GH3SS(2),GHMDC(2),GHMSU(2),GHMSC(2),GHMBT(2)
     +,GHPUS(2),GHPCD(2),GHPCS(2),GHPTB(2),GHMMUVM(2),GHMTAVT(2),GHPVMMU(2),GHPVTTA(2),GWWH1,GWWH2,GWWH3,GZZH1,GZZH2
     +,GZZH3,GZH1H3,GZH2H3,GZH1H2,GAHCHC,GZHCHC,GWHCH1,GWHCH2,GWHCH3,GH111,GH112,GH122,GH222,GH133,GH233,GH1HMHP,GH2
     +HMHP,GWWH1H1,GWWH2H2,GWWH3H3,GWWHCHC,GZZH1H1,GZZH2H2,GZZH3H3,GZZHCHC,Gggh1(2),Gggh2(2),Gggh3(2)
      common /DC_COUPL/ GH1LL,GH2LL,GH3LL,GH1MM,GH2MM,GH3MM,GH1TT,GH2TT,GH3TT,GH1CC,GH2CC,GH3CC,GH1BB,GH2BB,GH
     +3BB,GH1SS,GH2SS,GH3SS,GHMDC,GHMSU,GHMSC,GHMBT,GHPUS,GHPCD,GHPCS,GHPTB,GHMMUVM,GHMTAVT,GHPVMMU,GHPVTTA,GWWH1,GW
     +WH2,GWWH3,GZZH1,GZZH2,GZZH3,GZH1H3,GZH2H3,GZH1H2,GAHCHC,GZHCHC,GWHCH1,GWHCH2,GWHCH3,GH111,GH112,GH122,GH222,GH
     +133,GH233,GH1HMHP,GH2HMHP,GWWH1H1,GWWH2H2,GWWH3H3,GWWHCHC,GZZH1H1,GZZH2H2,GZZH3H3,GZZHCHC,Gggh1,Gggh2,Gggh3
      double precision MH1,MH2,MH3,MHC
      double precision WH1,WH2,WH3,WHC
      common /USR_MASS/  MH1,MH2,MH3,MHC
      common /USR_DECAY/ WH1,WH2,WH3,WHC
