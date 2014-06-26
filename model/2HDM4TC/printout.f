      subroutine printout
      implicit none
c
c     local
c
      integer i,iformat
      character*2 ab(2)
      real*8 ene
      double precision  Zero, One, Two, Three, Four, Half, Rt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
c
c     include
c
      include 'coupl.inc'
      include 'input.inc'
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS,msMS,mmuMS!MSbar masses
      double precision  Vud,Vus             !CKM matrix elements
      common/values/    alpha,gfermi,alfas,
     &                  mtMS,mbMS,mcMS,mtaMS,msMS,mmuMS,
     &                  Vud
c
c output all info
c
 10   format( 1x,a16,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV' )
 11   format( 1x,a13,2x,f11.5,2x,f11.5,2x,a13,2x,f11.5,2x,f11.5 )
 12   format( 1x,a13,2x,f6.2,a )
 13   format( 1x,a13,2x,f6.4,a )
 14   format( 1x,2(a13,2x,f10.7,2x,f10.7) )
 15   format( 1x,a13,2x,f9.5,a )
 16   format( 1x,a13,2x,f7.5 )
 17   format( 1x,a13,2x,f8.4 )
 18   format( 1x,a13,2x,f8.4,' GeV' )
 19   format( 1x,a13,2x,f6.4,a13,2x,f6.4 )
 20   format( 1x,a13,2x,f11.5,1x,f11.5 )
 21   format( 1x,a13,2x,f8.4,' GeV',1x,a13,2x,f8.4,' GeV' )
 22   format( 1x,a13,2x,f10.8,a13,2x,f6.4 )
 23   format( 1x,a13,2x,f8.4)
 24   format( 1x,a16,2x,f9.3,' GeV        ',a16,2x,f7.4,' GeV  (calc @ LO)')
 25   format( 1x,a16,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV')


      write(*,*) '*****************************************************'
      write(*,*) '*                    MadEvent                       *'
      write(*,*) '*        --------------------------------           *'
      write(*,*) '*          http://madgraph.hep.uiuc.edu             *'
      write(*,*) '*          http://madgraph.phys.ucl.ac.be           *'
      write(*,*) '*          http://madgraph.roma2.infn.it            *'
      write(*,*) '*        --------------------------------           *'
      write(*,*) '*                                                   *'
      write(*,*) '*         INTEGRATION CHANNEL LOG FILE              *'
      write(*,*) '*                                                   *'
      write(*,*) '*****************************************************'
      write(6,*)
      write(*,*) '*****************************************************'
      write(*,*) '*          SUMMARY OF THE SM PARAMETERS             *'
      write(*,*) '*****************************************************'
      write(6,*)
      write(6,*)  ' EW Parameters:'
      write(6,*)  '---------------'
      write(6,*)
      write(6,23) ' GF (10^-5*GeV^-2) = ',gfermi*1d5
      write(6,23) ' 1/alpha           = ',1d0/alpha
      write(6,23) ' M_Z   (GeV)       = ',zmass
      write(6,*)
      write(6,*)
      write(6,*)  'Boson masses and widths:'
      write(6,*)  '------------------------'
      write(6,*)
      write(6,24) 'Z mass  =  ',zmass, 'Z width  = ',zwidth
      write(6,24) 'W mass  =  ',wmass, 'W width  = ',wwidth
      write(6,*)
      write(6,*)  'Fermion masses and widths:'
      write(6,*)  '--------------------------'
      write(6,*)
      write(6,24) 'top     mass  =  ', tmass, 'top     width  = ', twidth
      write(6,10) 'bottom  mass  =  ', bmass, 'bottom  width  = ', Zero
      write(6,10) 'charm   mass  =  ', cmass, 'charm   width  = ', Zero
      write(6,10) 'strange mass  =  ', smass, 'strange width  = ', Zero
      write(6,10) 'mu      mass  =  ', mmass, 'mu      width  = ', Zero
      write(6,10) 'tau     mass  =  ', lmass, 'tau     width  = ', Zero
      write(6,*)  'all other quark and lepton masses set to zero'
      write(6,*)
      write(6,*)  'User Model couplings:'
      write(6,*)  '---------------------'
      write(6,*)
      write(6,11) 'GH1LL(L)  =  ', GH1LL(1),'GH1LL(R)  =  ', GH1LL(2)
      write(6,11) 'GH2LL(L)  =  ', GH2LL(1),'GH2LL(R)  =  ', GH2LL(2)
      write(6,11) 'GH3LL(L)  =  ', GH3LL(1),'GH3LL(R)  =  ', GH3LL(2)
      write(6,11) 'GH1MM(L)  =  ', GH1MM(1),'GH1MM(R)  =  ', GH1MM(2)
      write(6,11) 'GH2MM(L)  =  ', GH2MM(1),'GH2MM(R)  =  ', GH2MM(2)
      write(6,11) 'GH3MM(L)  =  ', GH3MM(1),'GH3MM(R)  =  ', GH3MM(2)
      write(6,11) 'GH1TT(L)  =  ', GH1TT(1),'GH1TT(R)  =  ', GH1TT(2)
      write(6,11) 'GH2TT(L)  =  ', GH2TT(1),'GH2TT(R)  =  ', GH2TT(2)
      write(6,11) 'GH3TT(L)  =  ', GH3TT(1),'GH3TT(R)  =  ', GH3TT(2)
      write(6,11) 'GH1CC(L)  =  ', GH1CC(1),'GH1CC(R)  =  ', GH1CC(2)
      write(6,11) 'GH2CC(L)  =  ', GH2CC(1),'GH2CC(R)  =  ', GH2CC(2)
      write(6,11) 'GH3CC(L)  =  ', GH3CC(1),'GH3CC(R)  =  ', GH3CC(2)
      write(6,11) 'GH1BB(L)  =  ', GH1BB(1),'GH1BB(R)  =  ', GH1BB(2)
      write(6,11) 'GH2BB(L)  =  ', GH2BB(1),'GH2BB(R)  =  ', GH2BB(2)
      write(6,11) 'GH3BB(L)  =  ', GH3BB(1),'GH3BB(R)  =  ', GH3BB(2)
      write(6,11) 'GH1SS(L)  =  ', GH1SS(1),'GH1SS(R)  =  ', GH1SS(2)
      write(6,11) 'GH2SS(L)  =  ', GH2SS(1),'GH2SS(R)  =  ', GH2SS(2)
      write(6,11) 'GH3SS(L)  =  ', GH3SS(1),'GH3SS(R)  =  ', GH3SS(2)
      write(6,11) 'GHMDC(L)  =  ', GHMDC(1),'GHMDC(R)  =  ', GHMDC(2)
      write(6,11) 'GHMSU(L)  =  ', GHMSU(1),'GHMSU(R)  =  ', GHMSU(2)
      write(6,11) 'GHMSC(L)  =  ', GHMSC(1),'GHMSC(R)  =  ', GHMSC(2)
      write(6,11) 'GHMBT(L)  =  ', GHMBT(1),'GHMBT(R)  =  ', GHMBT(2)
      write(6,11) 'GHPUS(L)  =  ', GHPUS(1),'GHPUS(R)  =  ', GHPUS(2)
      write(6,11) 'GHPCD(L)  =  ', GHPCD(1),'GHPCD(R)  =  ', GHPCD(2)
      write(6,11) 'GHPCS(L)  =  ', GHPCS(1),'GHPCS(R)  =  ', GHPCS(2)
      write(6,11) 'GHPTB(L)  =  ', GHPTB(1),'GHPTB(R)  =  ', GHPTB(2)
      write(6,11) 'GHMMUVM(L)  =  ', GHMMUVM(1),'GHMMUVM(R)  =  ', GHMMUVM(2)
      write(6,11) 'GHMTAVT(L)  =  ', GHMTAVT(1),'GHMTAVT(R)  =  ', GHMTAVT(2)
      write(6,11) 'GHPVMMU(L)  =  ', GHPVMMU(1),'GHPVMMU(R)  =  ', GHPVMMU(2)
      write(6,11) 'GHPVTTA(L)  =  ', GHPVTTA(1),'GHPVTTA(R)  =  ', GHPVTTA(2)
      write(6,20) 'GWWH1  =  ', GWWH1
      write(6,20) 'GWWH2  =  ', GWWH2
      write(6,20) 'GWWH3  =  ', GWWH3
      write(6,20) 'GZZH1  =  ', GZZH1
      write(6,20) 'GZZH2  =  ', GZZH2
      write(6,20) 'GZZH3  =  ', GZZH3
      write(6,20) 'GZH1H3  =  ', GZH1H3
      write(6,20) 'GZH2H3  =  ', GZH2H3
      write(6,20) 'GZH1H2  =  ', GZH1H2
      write(6,20) 'GAHCHC  =  ', GAHCHC
      write(6,20) 'GZHCHC  =  ', GZHCHC
      write(6,20) 'GWHCH1  =  ', GWHCH1
      write(6,20) 'GWHCH2  =  ', GWHCH2
      write(6,20) 'GWHCH3  =  ', GWHCH3
      write(6,20) 'GH111  =  ', GH111
      write(6,20) 'GH112  =  ', GH112
      write(6,20) 'GH122  =  ', GH122
      write(6,20) 'GH222  =  ', GH222
      write(6,20) 'GH133  =  ', GH133
      write(6,20) 'GH233  =  ', GH233
      write(6,20) 'GH1HMHP  =  ', GH1HMHP
      write(6,20) 'GH2HMHP  =  ', GH2HMHP
      write(6,20) 'GWWH1H1  =  ', GWWH1H1
      write(6,20) 'GWWH2H2  =  ', GWWH2H2
      write(6,20) 'GWWH3H3  =  ', GWWH3H3
      write(6,20) 'GWWHCHC  =  ', GWWHCHC
      write(6,20) 'GZZH1H1  =  ', GZZH1H1
      write(6,20) 'GZZH2H2  =  ', GZZH2H2
      write(6,20) 'GZZH3H3  =  ', GZZH3H3
      write(6,20) 'GZZHCHC  =  ', GZZHCHC
      write(6,11) 'Gggh1(L)  =  ', Gggh1(1),'Gggh1(R)  =  ', Gggh1(2)
      write(6,11) 'Gggh2(L)  =  ', Gggh2(1),'Gggh2(R)  =  ', Gggh2(2)
      write(6,11) 'Gggh3(L)  =  ', Gggh3(1),'Gggh3(R)  =  ', Gggh3(2)
      write(6,*)
      write(6,*) 'Boson couplings:'
      write(6,*) '----------------'
      write(6,*)
      write(6,20) 'gwwa  = ', gwwa
      write(6,20) 'gwwz  = ', gwwz
      write(6,*)
      write(6,*) 'FFV couplings:'
      write(6,*) '--------------'
      write(6,*)
      write(6,11) 'gal(L)   =  ',gal(1), 'gal(R)   =  ',gal(2)
      write(6,11) 'gau(L)   =  ',gau(1), 'gau(R)   =  ',gau(2)
      write(6,11) 'gad(L)   =  ',gad(1), 'gad(R)   =  ',gad(2)
      write(6,*)
      write(6,11) 'gwf(L)   =  ',gwf(1), 'gwf(R)   =  ',gwf(2)
      write(6,*)
      write(6,11) 'gzn(L)   =  ',gzn(1), 'gzn(R)   =  ',gzn(2)
      write(6,11) 'gzl(L)   =  ',gzl(1), 'gzl(R)   =  ',gzl(2)
      write(6,11) 'gzu(L)   =  ',gzu(1), 'gzu(R)   =  ',gzu(2)
      write(6,11) 'gzd(L)   =  ',gzd(1), 'gzd(R)   =  ',gzd(2)
      write(6,*)
      write(6,*) 'Strong couplings:'
      write(6,*) '-----------------'
      write(6,*)
      write(6,14) 'gg(1)    =  ',gg(1)   , 'gg(2)    =  ',gg(2)
      write(6,*)
      write(6,*) 'User Model masses and widths:'
      write(6,*) '-----------------------------'
      write(6,*)
      write(6,24) 'MH1  =  ', MH1, 'WH1    = ', WH1
      write(6,24) 'MH2  =  ', MH2, 'WH2    = ', WH2
      write(6,24) 'MH3  =  ', MH3, 'WH3    = ', WH3
      write(6,24) 'MHC  =  ', MHC, 'WHC    = ', WHC
      write(6,*)
      write(6,*) 'User Model Parameters:'
      write(6,*) '----------------------'
      write(6,24) 'TanBeta   =  ', TanBeta 
      write(6,24) 'SinAlpha   =  ', SinAlpha 
      write(6,24) 'SinGamma   =  ', SinGamma 
      write(6,*)

      return
      end

