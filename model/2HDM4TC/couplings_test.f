      subroutine testcoupl
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
      double precision  alpha, sin2w, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS,msMS,mmuMS!MSbar masses
      double precision  Vud,Vus             !CKM matrix elements
      common/values/    alpha,sin2w,gfermi,alfas,
     &                  mtMS,mbMS,mcMS,mtaMS,msMS,mmuMS,
     &                  Vud
      open(unit=1,file="couplings_check.txt")
c
c output all info
c
 10   format( 1x,a10,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV' )
 11   format( 1x,a10,2x,f11.5,2x,f11.5,a3,f11.5,2x,f11.5 )
 12   format( 1x,a10,2x,f6.2,a )
 13   format( 1x,a10,2x,f6.4,a )
 14   format( 1x,2(a10,2x,f10.7,2x,f10.7) )
 15   format( 1x,a10,2x,f9.5,a )
 16   format( 1x,a10,2x,f7.5 )
 17   format( 1x,a10,2x,f8.4 )
 18   format( 1x,a10,2x,f8.4,' GeV' )
 19   format( 1x,a10,2x,f6.4,a13,2x,f6.4 )
 20   format( 1x,a10,2x,f11.5,1x,f11.5 )
 21   format( 1x,a10,2x,f8.4,' GeV',1x,a13,2x,f8.4,' GeV' )
 22   format( 1x,a10,2x,f10.8,a13,2x,f6.4 )
 23   format( 1x,a10,2x,f8.4)
 24   format( 1x,a10,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV  (calc @ LO)')
 25   format( 1x,a10,2x,f7.3,' GeV        ',a16,2x,f7.4,' GeV')


      write(1,11) 'GH1LL ', GH1LL(1),' ', GH1LL(2)
      write(1,11) 'GH2LL ', GH2LL(1),' ', GH2LL(2)
      write(1,11) 'GH3LL ', GH3LL(1),' ', GH3LL(2)
      write(1,11) 'GH1MM ', GH1MM(1),' ', GH1MM(2)
      write(1,11) 'GH2MM ', GH2MM(1),' ', GH2MM(2)
      write(1,11) 'GH3MM ', GH3MM(1),' ', GH3MM(2)
      write(1,11) 'GH1TT ', GH1TT(1),' ', GH1TT(2)
      write(1,11) 'GH2TT ', GH2TT(1),' ', GH2TT(2)
      write(1,11) 'GH3TT ', GH3TT(1),' ', GH3TT(2)
      write(1,11) 'GH1CC ', GH1CC(1),' ', GH1CC(2)
      write(1,11) 'GH2CC ', GH2CC(1),' ', GH2CC(2)
      write(1,11) 'GH3CC ', GH3CC(1),' ', GH3CC(2)
      write(1,11) 'GH1BB ', GH1BB(1),' ', GH1BB(2)
      write(1,11) 'GH2BB ', GH2BB(1),' ', GH2BB(2)
      write(1,11) 'GH3BB ', GH3BB(1),' ', GH3BB(2)
      write(1,11) 'GH1SS ', GH1SS(1),' ', GH1SS(2)
      write(1,11) 'GH2SS ', GH2SS(1),' ', GH2SS(2)
      write(1,11) 'GH3SS ', GH3SS(1),' ', GH3SS(2)
      write(1,11) 'GHMDC ', GHMDC(1),' ', GHMDC(2)
      write(1,11) 'GHMSU ', GHMSU(1),' ', GHMSU(2)
      write(1,11) 'GHMSC ', GHMSC(1),' ', GHMSC(2)
      write(1,11) 'GHMBT ', GHMBT(1),' ', GHMBT(2)
      write(1,11) 'GHPUS ', GHPUS(1),' ', GHPUS(2)
      write(1,11) 'GHPCD ', GHPCD(1),' ', GHPCD(2)
      write(1,11) 'GHPCS ', GHPCS(1),' ', GHPCS(2)
      write(1,11) 'GHPTB ', GHPTB(1),' ', GHPTB(2)
      write(1,11) 'GHMMUVM ', GHMMUVM(1),' ', GHMMUVM(2)
      write(1,11) 'GHMTAVT ', GHMTAVT(1),' ', GHMTAVT(2)
      write(1,11) 'GHPVMMU ', GHPVMMU(1),' ', GHPVMMU(2)
      write(1,11) 'GHPVTTA ', GHPVTTA(1),' ', GHPVTTA(2)
      write(1,20) 'GWWH1 ', GWWH1
      write(1,20) 'GWWH2 ', GWWH2
      write(1,20) 'GWWH3 ', GWWH3
      write(1,20) 'GZZH1 ', GZZH1
      write(1,20) 'GZZH2 ', GZZH2
      write(1,20) 'GZZH3 ', GZZH3
      write(1,20) 'GZH1H3 ', GZH1H3
      write(1,20) 'GZH2H3 ', GZH2H3
      write(1,20) 'GZH1H2 ', GZH1H2
      write(1,20) 'GAHCHC ', GAHCHC
      write(1,20) 'GZHCHC ', GZHCHC
      write(1,20) 'GWHCH1 ', GWHCH1
      write(1,20) 'GWHCH2 ', GWHCH2
      write(1,20) 'GWHCH3 ', GWHCH3
      write(1,20) 'GH111 ', GH111
      write(1,20) 'GH112 ', GH112
      write(1,20) 'GH122 ', GH122
      write(1,20) 'GH222 ', GH222
      write(1,20) 'GH133 ', GH133
      write(1,20) 'GH233 ', GH233
      write(1,20) 'GH1HMHP ', GH1HMHP
      write(1,20) 'GH2HMHP ', GH2HMHP
      write(1,20) 'GWWH1H1 ', GWWH1H1
      write(1,20) 'GWWH2H2 ', GWWH2H2
      write(1,20) 'GWWH3H3 ', GWWH3H3
      write(1,20) 'GWWHCHC ', GWWHCHC
      write(1,20) 'GZZH1H1 ', GZZH1H1
      write(1,20) 'GZZH2H2 ', GZZH2H2
      write(1,20) 'GZZH3H3 ', GZZH3H3
      write(1,20) 'GZZHCHC ', GZZHCHC
      write(1,11) 'Gggh1 ', Gggh1(1),' ', Gggh1(2)
      write(1,11) 'Gggh2 ', Gggh2(1),' ', Gggh2(2)
      write(1,11) 'Gggh3 ', Gggh3(1),' ', Gggh3(2)
      write(1,20) 'gwwa ', gwwa
      write(1,20) 'gwwz ', gwwz
      write(1,11) 'gal ',gal(1),' ',gal(2)
      write(1,11) 'gau ',gau(1),' ',gau(2)
      write(1,11) 'gad ',gad(1),' ',gad(2)
      write(1,11) 'gwf ',gwf(1),' ',gwf(2)
      write(1,11) 'gzn ',gzn(1),' ',gzn(2)
      write(1,11) 'gzl ',gzl(1),' ',gzl(2)
      write(1,11) 'gzu ',gzu(1),' ',gzu(2)
      write(1,11) 'gzd ',gzd(1),' ',gzd(2)
      write(1,14) 'gg ',gg(1),' ',gg(2)

      return
      end
