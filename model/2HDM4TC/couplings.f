c----------------------------------------------------------------------
C  couplings.f 
c----------------------------------------------------------------------
c  This files takes the inputs of the standard model from a Les Houches 
c  file (param_card.dat) and calculates all the couplings that end up
c  in the Feynman rules, i.e., in the HELAS routines.
c   
c  With the convention of the New Web Generation MadEvent in this
c  part no non-trivial computation is made. The SM model secondary
c  parameters have been all calculated by the SM calculator, SMCalc
c  which produced param_card.dat.
c
c  The only exception to the above rule is for alpha_S. In the case
c  where a pdf is used in the initial state, then the values as(MZ)
c  set in the les houches card is superseeded.
c  Scales and pdf's are set in the run_card.dat.
c
c This file contains the following routines:
c 
c- madgraph original call to set the parameters
c- lh_readin is called from here.
c  It talks to the lh_readin through the common block values.
c      subroutine setpara
c
c-routine to read the LH parameters in
c      subroutine lh_readin
c
c-to set
c      subroutine set_it(npara,ivalue,value,name,id,
c      subroutine case_trap(string,length)
c      subroutine no_spaces(buff,nchars)
c---------------------------------------------------------------------- 


      subroutine setpara(param_name,readlha)
c***********************************************************************
c This subroutine sets up the HELAS couplings of the STANDARD MODEL.
c***********************************************************************
      implicit none
c
c local
c
      character*(*) param_name
      logical readlha
      integer i
      real*8 dum
c
c     common file with the couplings
c
      include 'coupl.inc'
      include 'input.inc'
c
c     local
c
      double precision  v
      double precision  ee, ee2, ez, ey, sw, cw, sc2, sin2w, wm
      double precision  gwne, gwud, lambda, lam4, xt, rew, rqcd
      double precision  alphas, alfa, alfaw, mfrun
      external          alphas, alfa, alfaw, mfrun
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS,mmuMS,msMS!MSbar masses
      double precision  Vud,Vus             !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,mmuMS,msMS,
     &                  Vud
c
c constants
c
      double complex  ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
      double precision  Zero, One, Two, Three, Four, Half, Rt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
      parameter( Three = 3.0d0, Four = 4.0d0, Half = 0.5d0 )
      parameter( Rt2   = 1.414213562d0 )
      double precision  Pi, Fourpi
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Fourpi = Four * Pi )	  
	  
c	  mmuMS=Zero
c	  msMS=Zero
c
c 2HDM parameters
c	  
	double precision alpha2HDM,beta2HDM, gamma2HDM
	double precision ca,cb,sa,sb,cab,sab,s2b,s2a
	double precision cag, sag, cbg, sbg
	double precision H1FF,H2FF,H3FF,HCFF
c
c    GGF values
c
	double precision  themass, shatval, thecoeff
	double precision  thelogterm, realterm, imageterm
	double precision  lambdaparam
c
c     alfas and run
c************
c Uncomment the following lines in order to use alphas from the PDF
c      include '../alfas.inc'
c      include '../run.inc'
c***********
c
c------------------------------------------
c Start calculating the couplings for HELAS
c------------------------------------------
c
      if(readlha) then 
         call lh_readin(param_name)
         G = DSQRT(4d0*PI*ALFAS) ! use setting of the param_card.dat @ NLO
      endif
c     

      GG(1) = -G
      GG(2) = -G     
c
c auxiliary local values
c
      wm = sqrt(zmass**2/Two+
     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2))
      sin2w  = One-(wm/zmass)**2
      cw  = sqrt( One - sin2w )
      ee2 = alpha * Fourpi
      sw  = sqrt( sin2w )
      ee  = sqrt( ee2 )
      ez  = ee/(sw*cw)
      ey  = ee*(sw/cw)
      sc2 = sin2w*( One - sin2w )
      v   = Two*wm*sw/ee   ! the wmass is used to calculate v
c
c vector boson couplings
c
      gw   = ee/sw
      gwwa = ee
      gwwz = ee*cw/sw
c
c fermion-fermion-vector couplings
c
      gal(1) = dcmplx(  ee          , Zero )
      gal(2) = dcmplx(  ee          , Zero )
      gau(1) = dcmplx( -ee*Two/Three, Zero )
      gau(2) = dcmplx( -ee*Two/Three, Zero )
      gad(1) = dcmplx(  ee/Three    , Zero )
      gad(2) = dcmplx(  ee/Three    , Zero )

      gwf(1) = dcmplx( -ee/sqrt(Two*sin2w), Zero )
      gwf(2) = dcmplx(  Zero              , Zero )

      gzn(1) = dcmplx( -ez*Half                     , Zero )
      gzn(2) = dcmplx(  Zero                        , Zero )
      gzl(1) = dcmplx( -ez*(-Half + sin2w)          , Zero )
      gzl(2) = dcmplx( -ey                          , Zero )
      gzu(1) = dcmplx( -ez*( Half - sin2w*Two/Three), Zero )
      gzu(2) = dcmplx(  ey*Two/Three                , Zero )
      gzd(1) = dcmplx( -ez*(-Half + sin2w/Three)    , Zero )
      gzd(2) = dcmplx( -ey/Three                    , Zero )


c     Full CKM code currently unwritten     
c
c     CKM matrix: 
c     symmetric 3x3 matrix, Vud=Vcs, Vus=Vcd Vcb=Vub=0
c
c     >>>>>>>>>>>>>>>***** NOTE****<<<<<<<<<<<<<<<<<<<<<<<<<
c     these couplings matter only when interaction_CKM.dat
c     is used to generate all the diagrams with off-diagonal
c     couplings. The default of MadEvent is a diagonal
c     CKM matrix.

	  Vus=DSQRT(1d0-Vud**2)
      do i=1,2
         gwfc(i) = gwf(i)*Vud
         gwfs(i) = gwf(i)*Vus
         gwfm(i) =-gwf(i)*Vus
      enddo

c---------------------------------------------------------
c Set Photon Width to Zero, used by symmetry optimization
c---------------------------------------------------------

      awidth = 0d0
c************************************            
c 2HDM Couplings
c************************************
c Parameter Definitions
      beta2HDM = atan(tanbeta)
	  alpha2HDM = asin(sinalpha)	
	  gamma2HDM = asin(singamma)	  
	  
	  sb = sin(beta2HDM)
	  cb = cos(beta2HDM)
	  sa = sin(alpha2HDM)
	  ca = cos(alpha2HDM)
	  sab = sin(alpha2HDM-beta2HDM)
	  cab = cos(alpha2HDM-beta2HDM)
	  s2a = sin(2*alpha2HDM)
	  s2b = sin(2*beta2HDM)
	  
	  sbg = sin(gamma2HDM+beta2HDM)
	  cbg = cos(gamma2HDM+beta2HDM)
	  sag = sin(alpha2HDM+gamma2HDM)
	  cag = cos(alpha2HDM+gamma2HDM)
	  
	  H1FF = cag/sbg
	  H2FF = sag/sbg
	  H3FF = cbg/sbg
	  HCFF = cbg/sbg

c Fermion-fermion Higgs Coupling	
c Neutral  
      if(mtMS.gt.1d0) then
         gh1tt(1) = dcmplx( -mtMS/v*H1FF, Zero )
		 gh2tt(1) = dcmplx( -mtMS/v*H2FF, Zero )
		 gh3tt(1) = dcmplx( -mtMS/v*H3FF, Zero )
      else
         gh1tt(1) = dcmplx( Zero,Zero)
		 gh2tt(1) = dcmplx( Zero,Zero)
		 gh3tt(1) = dcmplx( Zero,Zero)
      endif
      gh1tt(2) = gh1tt(1)
	  gh2tt(2) = gh2tt(1)
	  gh3tt(2) = -gh3tt(1)

	  if(mbMS.gt.1d0) then
         gh1bb(1) = dcmplx( -mbMS/v*H1FF, Zero )
		 gh2bb(1) = dcmplx( -mbMS/v*H2FF, Zero )
		 gh3bb(1) = dcmplx( -mbMS/v*H3FF, Zero )
      else
         gh1bb(1) = dcmplx( Zero,Zero)
		 gh2bb(1) = dcmplx( Zero,Zero)
		 gh3bb(1) = dcmplx( Zero,Zero)
      endif
      gh1bb(2) = gh1bb(1)
	  gh2bb(2) = gh2bb(1)
	  gh3bb(2) = -gh3bb(1)
	  
	  if(mcMS.gt.1d-2) then
         gh1cc(1) = dcmplx( -mcMS/v*H1FF, Zero )
		 gh2cc(1) = dcmplx( -mcMS/v*H2FF, Zero )
		 gh3cc(1) = dcmplx( -mcMS/v*H3FF, Zero )
      else
         gh1cc(1) = dcmplx( Zero,Zero)
		 gh2cc(1) = dcmplx( Zero,Zero)
		 gh3cc(1) = dcmplx( Zero,Zero)
      endif
      gh1cc(2) = gh1cc(1)
	  gh2cc(2) = gh2cc(1)
	  gh3cc(2) = -gh3cc(1)

	  if(msMS.gt.1d-2) then
         gh1ss(1) = dcmplx( -msMS/v*H1FF, Zero )
		 gh2ss(1) = dcmplx( -msMS/v*H2FF, Zero )
		 gh3ss(1) = dcmplx( -msMS/v*H3FF, Zero )
      else
         gh1ss(1) = dcmplx( Zero,Zero)
		 gh2ss(1) = dcmplx( Zero,Zero)
		 gh3ss(1) = dcmplx( Zero,Zero)
      endif
      gh1ss(2) = gh1ss(1)
	  gh2ss(2) = gh2ss(1)
	  gh3ss(2) = -gh3ss(1)
	  
	  if(mmuMS.gt.1d-2) then
	    gh1mm(1) = dcmplx( -mmuMS/v*H1FF, Zero )
	    gh2mm(1) = dcmplx( -mmuMS/v*H2FF, Zero )
	    gh3mm(1) = dcmplx( -mmuMS/v*H3FF, Zero )
	  else
	    gh1mm(1) = dcmplx( Zero, Zero )
	    gh2mm(1) = dcmplx( Zero, Zero )
	    gh3mm(1) = dcmplx( Zero, Zero )
	  endif 
	  gh1mm(2) = gh1mm(1)
	  gh2mm(2) = gh2mm(1)
	  gh3mm(2) = -gh3mm(1) 
	  
	  if(mtaMS.gt.1d-2) then
	    gh1ll(1) = dcmplx( -mtaMS/v*H1FF, Zero )
	    gh2ll(1) = dcmplx( -mtaMS/v*H2FF, Zero )
	    gh3ll(1) = dcmplx( -mtaMS/v*H3FF, Zero )
	  else
	    gh1ll(1) = dcmplx( Zero, Zero )
	    gh2ll(1) = dcmplx( Zero, Zero )
	    gh3ll(1) = dcmplx( Zero, Zero )
	  endif
	  gh1ll(2) = gh1ll(1)
	  gh2ll(2) = gh2ll(1)
	  gh3ll(2) = -gh3ll(1) 
	  
c  Charged	  
	  ghmdc(1) = dcmplx(Zero,Zero)
	  ghmsu(1) = dcmplx(Zero,Zero)
	  ghmsc(1) = dcmplx(-(msMS+mcMS)/v*HCFF,Zero)
	  ghmbt(1) = dcmplx(-(mtMS+mbMS)/v*HCFF,Zero)
	  ghpus(1) = dcmplx(Zero,Zero)
	  ghpcd(1) = dcmplx(Zero,Zero)
	  ghpcs(1) = dcmplx(-(msMS+mcMS)/v*HCFF,Zero)
	  ghptb(1) = dcmplx(-(mtMS+mbMS)/v*HCFF,Zero)
	  ghmdc(2) = ghmdc(1)
	  ghmsu(2) = ghmsu(1)
	  ghmsc(2) = ghmsc(1)
	  ghmbt(2) = ghmbt(1)
	  ghpus(2) = ghpus(1)
	  ghpcd(2) = ghpcd(1)
	  ghpcs(2) = ghpcs(1)
	  ghptb(2) = ghptb(1)
	  ghmmuvm(1) = dcmplx(-mmuMS/v*HCFF,Zero)
	  ghmtavt(1) = dcmplx(-mtaMS/v*HCFF,Zero)
	  ghpvmmu(1) = dcmplx(-mmuMS/v*HCFF,Zero)
	  ghpvtta(1) = dcmplx(-mtaMS/v*HCFF,Zero)
	  ghmmuvm(2) = ghmmuvm(1)
	  ghmtavt(2) = ghmtavt(1)
	  ghpvmmu(2) = ghpvmmu(1)
	  ghpvtta(2) = ghpvtta(1)
	  
c
c  VVS couplings
c
      gwwh1  = dcmplx( ee2/sin2w*Half*v*sab, Zero )
      gwwh2  = dcmplx( ee2/sin2w*Half*v*cab, Zero )
      gwwh3  = dcmplx( Zero, Zero )
      gzzh1  = dcmplx( ee2/sc2*Half*v*sab, Zero )
      gzzh2  = dcmplx( ee2/sc2*Half*v*cab, Zero )
      gzzh3  = dcmplx( Zero, Zero )

c
c  VSS couplings
c	  
	  gzh1h3 = dcmplx(-Half*ez*cab,Zero)
	  gzh2h3 = dcmplx(-Half*ez*sab,Zero)
	  gzh1h2 = dcmplx(Zero,Zero) 
	  
c   charged higgs	  
	  gahchc = dcmplx(Zero,Zero)
	  gzhchc = dcmplx(Zero,Half*ez)
	  gwhch1 = dcmplx(Zero,Half*gw*cab)
	  gwhch2 = dcmplx(Zero,Half*gw*sab)
	  gwhch3 = dcmplx(Zero,Half*gw)
	  
	
c
c  SSS couplings
c
	  gh111 = Three*MH1**2 * (cos(Three**alpha2HDM-beta2HDM)+Three*cos(alpha2HDM+beta2HDM))/(Two*v*s2b)
	  gh112 = (Two*MH1**2 + MH2**2)*cab*s2a/(v*s2b)
	  gh122 = (MH1**2 + Two*MH2**2)*sab*s2a/(v*s2b)
	  gh222 = Three*MH2**2 * (sin(Three**alpha2HDM-beta2HDM)-Three*sin(alpha2HDM+beta2HDM))/(2*v*s2b)
	  gh133 = ( (MH1**2 - Two*MH3**2)*cos(alpha2HDM-Three*beta2HDM) +(Three*MH1**2 + Two*MH3**2)*cos(alpha2HDM+beta2HDM) )/(Two*v*s2b)
	  gh233 = ( (MH2**2 - Two*MH3**2)*sin(alpha2HDM-Three*beta2HDM) +(Three*MH2**2 + Two*MH3**2)*sin(alpha2HDM+beta2HDM) )/(Two*v*s2b)
	  gh1hmhp = ( (MH1**2 - Two*MH3**2)*cos(alpha2HDM-Three*beta2HDM) +(Three*MH1**2 + Two*MH3**2)*cos(alpha2HDM+beta2HDM) )/(Two*v*s2b)
	  gh2hmhp = ( (MH2**2 - Two*MH3**2)*sin(alpha2HDM-Three*beta2HDM) +(Three*MH2**2 + Two*MH3**2)*sin(alpha2HDM+beta2HDM) )/(Two*v*s2b)
	  
c VVSS 
	  gwwh1h1 = gw*gw*Half
	  gwwh2h2 = gw*gw*Half
	  gwwh3h3 = gw*gw*Half
	  gwwhchc = gw*gw*Half
	  gzzh1h1 = ez*ez*Half
	  gzzh2h2 = ez*ez*Half
	  gzzh3h3 = ez*ez*Half
	  gzzhchc = ez*ez*Half
	  
c **********************************************************************
c *                          GGH Couplings                             *
c **********************************************************************
c Following code permits access to momenta of particles in event
c Only needed for GGF couplings, but even there still not essential
c See below for comments
c
c Comment out when running 'make testprog', then uncomment after made 
c broad (three lines)
c      include '../genps.inc'
c      double precision pp(0:3,max_particles)	  
c      common/momenta_pp/pp 
c
c Scalar and Pseudoscalar glue-glue couplings 
c Taken from Gunion & Haber - Higgs Bosons in Supersymmetic Models (II)
c Dimension-5 couplings, only exact in collider GGF process without ISR 
c Presented here with shat rather than resonance mass 
c Using resonance mass will be the same except for broad widths
c If the broad width implementation is desired uncomment the lines 
c preceded by with "c broad", then comment the lines preceded by 
c "c narrow", for narrow width, can be used as is.
c
c   Comment out when running 'make testprog', then uncomment after made
c broad (one line)
c	  shatval = Two*(pp(0,1)*pp(0,2) - pp(3,1)*pp(3,2))	  

c     H1 
c narrow (one line)
	  shatval = MH1**2
	  lambdaparam =  mtMS**2 / shatval
	  thecoeff = G**2/Fourpi**2 * lambdaparam/v * H1FF * Four
	  if (lambdaparam < 0.25) then
		thelogterm = log(( One + sqrt(One - Four * lambdaparam))/ ( One - sqrt(One - Four * lambdaparam)))
		realterm = thecoeff*( Half*(One - Four * lambdaparam)*(thelogterm**2 - Pi**2) - Two)
		imageterm = thecoeff*thelogterm*Pi*(One - Four * lambdaparam)
	  else
		realterm = Two*thecoeff*((Four * lambdaparam - One)*asin(One/(Two*sqrt(lambdaparam)))**2 - One)
		imageterm = Zero
	  endif
	  gggh1(1)=dcmplx( realterm,imageterm)	
	  gggh1(2)=dcmplx(Zero,Zero)

c     H2 
c narrow (one line) 
	  shatval = MH2**2
	  lambdaparam =  mtMS**2 / shatval
	  thecoeff = G**2/Fourpi**2 * lambdaparam/v * H2FF * Four
	  if (lambdaparam < 0.25) then
		thelogterm = log(( One + sqrt(One - Four * lambdaparam))/ ( One - sqrt(One - Four * lambdaparam)))
		realterm = thecoeff*( Half*(One - Four * lambdaparam)*(thelogterm**2 - Pi**2) - Two)
		imageterm = thecoeff*thelogterm*Pi*(One - Four * lambdaparam)
	  else
		realterm = Two*thecoeff*((Four * lambdaparam - One)*asin(One/(Two*sqrt(lambdaparam)))**2 - One)
		imageterm = Zero
	  endif
	  gggh2(1)=dcmplx( realterm,imageterm)	
	  gggh2(2)=dcmplx(Zero,Zero)

c     H3 
c narrow (one line)	  
      shatval = MH3**2
	  lambdaparam = mtMS**2 / shatval
      thecoeff = G**2/Fourpi**2 * lambdaparam/v  * H3FF * Four
	  if (lambdaparam < 0.25) then
	    thelogterm = log(( One + sqrt(One - Four * lambdaparam)) / ( One - sqrt(One - Four * lambdaparam)))
        gggh3(1)=dcmplx(Zero,Zero)
	    gggh3(2)=dcmplx(thecoeff*Half*(thelogterm**2 - Pi**2) ,thecoeff*thelogterm*Pi)
	  else
	    gggh3(1)=dcmplx(Zero,Zero)
		gggh3(2)=dcmplx(-Two*thecoeff*asin(One/(Two*sqrt(lambdaparam)))**2,Zero)
	  endif
		
c----------------------------
c end subroutine coupsm
c----------------------------


      return
      end

      
