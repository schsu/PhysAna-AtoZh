//#############################################################
//#                                                          ##
//#                  2HDM for TC Calculator                  ##
//#                                                          ##
//# FILE : 2HDM4TCCalc.c                                     ##
//# VERSION : 1.0                                            ##
//# DATE : 25 August 2010                                    ##
//# AUTHOR : Jared A. Evans (UC Davis)                       ##
//#                                                          ##
//# Description:  Convert a file into a param_card for       ##
//#               use with the 2HDM4TC model                 ##
//#                                                          ##
//# Note:  Code heavily borrowed and adapted from            ##
//#        M. Herquet's TwoHiggsCalc.c                       ##
//#                                                          ##
//#############################################################


// GENERAL INCLUDE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <complex.h>
#include <time.h>

// INPUT FILE PARAM
// Maximum length of a input file line
#define MAXLG 1000
// Maximum length of model name
#define MAXBLNM 20

// MATH
#define PI M_PI
#define Sqrt2 sqrt(2)

// NUMERICAL TOLERANCE FOR TESTS
#define TOL pow(10,-3)


// FUNCTION GETRAND
// ================
double getrand(int acc)
{
	//time based random number generator
	//not great, but very fast and good enough for our purposes
	double x=0; int i=0;
	for(i = 0; i<acc;i++)
	{
		x += (rand() % 10)* pow(10,i);
	}
	x = x/pow(10,acc);
	return x;
}

// FUNCION KAPPA (KINEMATICS)
// ==========================
double kappa(double x, double y, double z) 
{
	return sqrt(pow(x,4)+pow(y,4)+pow(z,4)-2.*(pow(x,2)*pow(y,2)+pow(x,2)*pow(z,2)+pow(y,2)*pow(z,2)));
}

// FUNCION KINSVV (KINEMATICS)
// ===========================
// Kinematic parameter for massive vectors
double kinsvv(double MS, double MV1, double MV2) 
{
	return 1.+pow(pow(MS,2)-pow(MV1,2)-pow(MV2,2),2)/(8.*pow(MV1*MV2,2));
}

// FUNCTION R (SVV below threshold)
// ================================

double R(double x) {
	return 3.*(1.-8.*x+20.*pow(x,2.))/sqrt(4.*x-1.)*acos((3.*x-1.)/(2.*pow(x,1.5)))-(1.-x)/(2.*x)*(2.-13.*x+47.*pow(x,2.))-1.5*(1.-6.*x+4.*pow(x,2.))*log(x);
}

// FUNCTION WSFF
// =============
// PARAM: number of color (3 or 1), scalar mass, fermions masses, coupling, vector or axial(0 or 1)
// RETURN: decay width of a scalar into fermions

double WSFF(double N, double MS, double MF1, double MF2, double COUPL, double AV) {
	
	if (MS==0.) return 0.;
	if (MF1+MF2>MS) return 0.;
	return (N*pow(COUPL,2)/(8.*PI))*kappa(MS,MF1,MF2)*(pow(MS,2)-pow(MF1+2.*(0.5-AV)*MF2,2))/pow(MS,3);
}


// FUNCTION WSVV
// =============

// PARAM: dv=1 for Z and 2 for W, scalar mass, boson mass, coupling, sin theta_W
// RETURN: decay width of a scalar into vector bosons

double WSVV(double dv, double MS, double MV, double COUPL, double sw) {
	
	if (MS==0.) return 0.;
	double x=pow(MV/MS,2.);
	double dvp=0.;
	if (MS<MV) return 0.;	
	else if (MS>MV && 2.*MV>MS) {
		if (dv==1.) dvp=7./12.-10.*pow(sw,2.)/9.+40.*pow(sw,4.)/27.;
		if (dv==2.) dvp=1.; 
		return dvp*pow(COUPL,2.)*(3./(64.*16.*pow(PI,3.)*pow(MV,4.)))*MS*R(x);
	}
	return dv*pow(COUPL,2)/(16.*PI)*kappa(MS,MV,MV)*kinsvv(MS,MV, MV)/pow(MS,3);
}

// FUNCTION WSVS
// =============

// PARAM: scalar masses, vector boson mass, coupling
// RETURN: decay width of a scalar into scalar/vector boson

double WSVS(double MS1, double MS2, double MV, double COUPL) 
{
	if (MS1==0.) return 0.;
	if (MS2+MV>MS1) return 0.;
	if (MV==0.) return (pow(COUPL,2)/(8.*PI))*MS1*(1-pow(MS2/MS1,4));
	return (pow(COUPL,2)/(16.*PI))/(pow(MV,2)*pow(MS1,3))*pow(kappa(MS1,MV,MS2),3);
}

// FUNCTION WSSS
// =============

// PARAM: d=1 if id. part. else d=2, scalar masses, coupling
// RETURN: decay width of a scalar into scalars

double WSSS(double d, double MS, double MS1, double MS2, double COUPL) 
{
	if (MS==0.) return 0.;
	if (MS1+MS2>MS) return 0.;
	return d/(32.*PI)*pow(COUPL,2)/pow(MS,3)*kappa(MS1,MS2,MS);
}

// FUNCTION MAXM23
// ==============
// evaluates the kinematic parameter
double MAXM23(double x, double Msq, double m1sq, double m2sq, double m3sq)
{
	double E2s=0.,E3s=0.;
	E2s = (x-m1sq+m2sq)/(2.*sqrt(x));
	E3s = (Msq-x-m3sq)/(2.*sqrt(x));
	
	return pow(E2s+E3s,2)-pow(sqrt(pow(E2s,2)-m2sq) - sqrt(pow(E3s,2)-m3sq),2);
}

// FUNCTION MINM23
// ==============
// evaluates the kinematic parameter
double MINM23(double x, double Msq, double m1sq, double m2sq, double m3sq)
{
	double E2s=0.,E3s=0.;
	E2s = (x-m1sq+m2sq)/(2.*sqrt(x));
	E3s = (Msq-x-m3sq)/(2.*sqrt(x));
	
	return pow(E2s+E3s,2)-pow(sqrt(pow(E2s,2)-m2sq) + sqrt(pow(E3s,2)-m3sq),2);
}


// FUNCTION XZZZ
// ==============

// evaluates the integral for the SZZZ decay
double XZZZ(double x, double y, double z, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double a=0,b=0;
	
	a = (x*y*(1.-x)*(1.-y))/(4.*pow(z,3));
	a = a + (-2.+2.*(x+y)-2.*pow(x+y,2)+3.*x*y*(x+y))/(4.*pow(z,2));
	a = a + (6.+4.*(x+y)+2.*pow(x+y,2)+5.*x*y)/(4.*z);
	a = a + 5./2.*z - (7.+3.*(x+y))/2.;
	
	b = 2.*ghl*ghl*((x-hl)*(y-hl)+hl*whl)/(pow((x-hl)*(y-hl)+hl*whl,2)+hl*whl*pow(x-y,2));
	b = b + 2.*ghh*ghh*((x-hh)*(y-hh)+hh*whh)/(pow((x-hh)*(y-hh)+hh*whh,2)+hh*whh*pow(x-y,2));
	b = b + ghl*ghh*((x-hl)*(y-hh)+sqrt(hh*whh*hl*whl))/(pow((x-hl)*(y-hh)+sqrt(hh*whh*hl*whl),2)+pow(sqrt(hl*whl)*(y-hh) - sqrt(hh*whh)*(x-hl),2) );
	b = b + ghl*ghh*((y-hl)*(x-hh)+sqrt(hh*whh*hl*whl))/(pow((y-hl)*(x-hh)+sqrt(hh*whh*hl*whl),2)+pow(sqrt(hl*whl)*(x-hh) - sqrt(hh*whh)*(y-hl),2) );	
	
	return a * b;
}

// FUNCTION XWWZ
// ==============

// evaluates the integral for the SWWZ decay
double XWWZ(double x,double z, double w, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double a=0,b=0,c=0;
	a = (pow(x/(2.*w),2)-x/w+3.)*(pow(1.-z-x,2)/z-4.*x);
	b = ghl*ghl/(pow(x-hl,2)+hl*whl)+ghh*ghh/(pow(x-hh,2)+hh*whh);
	b = b + 2.*ghl*ghh*((x-hl)*(x-hh)+sqrt(hl*hh*whl*whh))/((pow(x-hl,2)+hl*whl)*(pow(x-hh,2)+hh*whh));
	c = MAXM23(x,1.,z,w,w) - MINM23(x,1.,z,w,w);
	return a * b * c;
}

// FUNCTION XVFF
// ==============

// evaluates the integral for the SVFF decay

double XVFF(double x, double v, double f1, double f2, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double a=0,b=0,c=0;
	a = (pow(1.-v-x,2)/v-4.*x)*(x-pow(f1+f2,2));
	b = ghl*ghl/(pow(x-hl,2)+hl*whl)+ghh*ghh/(pow(x-hh,2)+hh*whh);
	b = b + 2.*ghl*ghh*((x-hl)*(x-hh)+sqrt(hl*hh*whl*whh))/((pow(x-hl,2)+hl*whl)*(pow(x-hh,2)+hh*whh));
	c = MAXM23(x,1.,v,pow(f1,2),pow(f2,2)) - MINM23(x,1.,v,pow(f1,2),pow(f2,2));
	return a * b * c;
}

// FUNCTION XSFF
// ==============

// evaluates the integral for the SSFF decay

double XSFF(double x, double s, double f1, double f2, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double a=0,b=0,c=0;
	a = x-pow(f1+f2,2);
	b = ghl*ghl/(pow(x-hl,2)+hl*whl)+ghh*ghh/(pow(x-hh,2)+hh*whh);
	b = b + 2.*ghl*ghh*((x-hl)*(x-hh)+sqrt(hl*hh*whl*whh))/((pow(x-hl,2)+hl*whl)*(pow(x-hh,2)+hh*whh));
	c = MAXM23(x,1.,s,pow(f1,2),pow(f2,2)) - MINM23(x,1.,s,pow(f1,2),pow(f2,2));
	return a * b * c;
}

// FUNCTION XSVV
// ==============

// evaluates the integral for the SSVV decay
double XSVV(double x, double s, double w, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double a=0,b=0,c=0;
	a = (pow(x/(2.*w),2)-x/w+3.);
	b = ghl*ghl/(pow(x-hl,2)+hl*whl)+ghh*ghh/(pow(x-hh,2)+hh*whh);
	b = b + 2.*ghl*ghh*((x-hl)*(x-hh)+sqrt(hl*hh*whl*whh))/((pow(x-hl,2)+hl*whl)*(pow(x-hh,2)+hh*whh));
	c = MAXM23(x,1.,s,w,w) - MINM23(x,1.,s,w,w);
	return a * b * c;
}

// FUNCTION XSVV2
// ==============

// evaluates the integral for the SVFF decay

double XSVV2(double x, double s, double w, double h, double wh, double gh)
{
	double a=0,b=0,c=0;
	a = (pow(1.-w-x,2)/w-4.*x)*(pow(x-w-s,2)/w-4.*s);
	b = gh*gh/(pow(x-h,2)+h*wh);
	c = MAXM23(x,1.,w,w,s) - MINM23(x,1.,w,w,s);
	return a * b * c;
}
// FUNCTION XSVV2INT
// =================

// evaluates the integral for the SVFF decay

double XSVV2INT(double x, double s, double w, double h, double wh, double gh)
{
	double a=0,b=0,c=0;
	a = (pow(1.-w-x,2)/w-4.*x)*(pow(x-w-s,2)/w-4.*s);
	b = gh*gh/(pow(x-h,2)+h*wh);
	c = MAXM23(x,1.,w,w,s) - MINM23(x,1.,w,w,s);
	return a * b * c;
}

// FUNCTION MCWWZ
// ==============

// evaluates the integral for the WWZ decay
double MCWWZ(double z, double w, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double xpos, amppos, amp, maxamp, intmax, intmin, xscale;
	
	intmax = pow(1.-sqrt(z),2);
	intmin = 4.*w;
	
	xscale = intmax-intmin;
	maxamp = 0;
	//searches for a good initial cap 
	xpos = intmin;
	double step = xscale/1000;
	double scale = 1.5;
	for(xpos=intmin;xpos<intmax;xpos+=step)
	{
		amp = XWWZ(xpos,z,w,hl,hh,whl,whh,ghl,ghh);
		if(amp>maxamp) maxamp=amp;
	}
	maxamp = maxamp*scale;
	
	//MC begins here
	int nstep = 10000,i=0,ng=0,nt=0;
	double lastintegral=0;
	double integral=0;
	int randp = 8;
	double desiredprecision = pow(10.,-4);
	
	for(i=0;i<nstep;i++)
	{
		xpos = getrand(randp)*xscale+intmin;
		amppos = getrand(randp)*maxamp;
		amp = XWWZ(xpos,z,w,hl,hh,whl,whh,ghl,ghh);
		if (amppos<amp) ng++;
		nt++;
		if(i==nstep-1)
		{
			integral = ng*1./(nt*1.);
			if((integral-lastintegral)/integral > desiredprecision)
			{
				i=0;
				lastintegral=integral;
			}
		}
		if(amp>maxamp)
		{
			maxamp = amp*scale;
			i=0; ng=0; nt=0; lastintegral=0;
		}
	}
	integral = ng*1./(nt*1.)*xscale*maxamp;
	return integral;
}

// FUNCTION MCZZZ
// ==============

// evaluates the integral of the interference terms in the ZZZ or WWW decays
double MCZZZ(double z, double w, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double xpos, ypos, amppos, amp, maxamp, xmax, xmin, xscale, yscale, ymax, ymin, area;
	
	xmax = pow(1.-sqrt(z),2);
	xmin = 4.*w;
	xscale = xmax-xmin;
	
	//gets area and creates an initial cap
	int nstep = 10000,i=0,ng=0,nt=0;
	double lastintegral=0;
	double integral=0;
	int randp = 10;
	double desiredprecision = pow(10.,-4);
	double scale = 1.5;
	maxamp=0.;
	
	for(i=0;i<nstep;i++)
	{
		xpos = getrand(randp)*xscale+xmin;
		ypos = getrand(randp);
		ymax = MAXM23(xpos,1.,z,w,w);
		ymin = MINM23(xpos,1.,z,w,w);
		if (ymin<ypos<ymax)
		{
			ng++;
			amppos = XZZZ(xpos,ypos,z,hl,hh,whl,whh,ghl,ghh);
			if (amppos>maxamp) maxamp = amppos*scale;
		}
		nt++;
		if(i==nstep-1)
		{
			integral = ng*1./(nt*1.);
			if((integral-lastintegral)/integral > desiredprecision)
			{
				i=0;
				lastintegral=integral;
			}
		}
	}
	area = ng*1./(nt*1.)*xscale*1.;
	
	//MC begins here
	lastintegral=0;
	integral=0;
	for(i=0;i<nstep;i++)
	{
		xpos = getrand(randp)*xscale+xmin;
		ymax = MAXM23(xpos,1.,z,w,w);
		ymin = MINM23(xpos,1.,z,w,w);
		yscale = ymax-ymin;
		ypos = getrand(randp)*yscale+ymin;
		amppos = getrand(randp)*maxamp;
		amp = XZZZ(xpos,ypos,z,hl,hh,whl,whh,ghl,ghh);
		if (amppos<amp) ng++;
		nt++;
		if(i==nstep-1)
		{
			integral = ng*1./(nt*1.);
			if((integral-lastintegral)/integral > desiredprecision)
			{
				i=0;
				lastintegral=integral;
			}
		}
		if(amp>maxamp)
		{
			maxamp = amp*scale;
			i=0; ng=0; nt=0; lastintegral=0;
		}
	}
	integral = ng*1./(nt*1.)*area*maxamp;
	return integral;
}

// FUNCTION MCVFF
// ==============

// evaluates the integral for the VFF decays
double MCVFF(double v, double f1, double f2, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double xpos, amppos, amp, maxamp, intmax, intmin, xscale;
	
	intmax = pow(1.-sqrt(v),2);
	intmin = pow(f1+fabs(f2),2);
	
	xscale = intmax-intmin;
	maxamp = 0;
	
	//searches for a good initial cap 
	xpos = intmin;
	double step = xscale/1000;
	double scale = 1.5;
	for(xpos=intmin;xpos<intmax;xpos+=step)
	{
		amp = XVFF(xpos,v,f1,f2,hl,hh,whl,whh,ghl,ghh);
		if(amp>maxamp) maxamp=amp;
	}
	maxamp = maxamp*scale;
	
	//MC begins here
	long nstep = 10000,i=0,ng=0,nt=0;
	double lastintegral=0;
	double integral=0;
	int randp = 8;
	double desiredprecision = pow(10.,-4);
	
	for(i=0;i<nstep;i++)
	{
		xpos = getrand(randp)*xscale+intmin;
		amppos = getrand(randp)*maxamp;
		amp = XVFF(xpos,v,f1,f2,hl,hh,whl,whh,ghl,ghh);
		if (amppos<amp) ng++;
		nt++;
		if(i==nstep-1)
		{
			integral = ng*1./(nt*1.);
			if((integral-lastintegral)/integral > desiredprecision)
			{
				i=0;
				lastintegral=integral;
			}
		}
		if(amp>maxamp)
		{
			maxamp = amp*scale;
			i=0; ng=0; nt=0; lastintegral=0;
		}
	}
	integral = ng*1./(nt*1.)*xscale*maxamp;
	return integral;
}

// FUNCTION MCSFF
// ==============

// evaluates the integral for the VFF decays
double MCSFF(double s, double f1, double f2, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double xpos, amppos, amp, maxamp, intmax, intmin, xscale;
	
	intmax = pow(1.-sqrt(s),2);
	intmin = pow(f1+fabs(f2),2);
	
	xscale = intmax-intmin;
	maxamp = 0;
	
	//searches for a good initial cap 
	xpos = intmin;
	double step = xscale/1000;
	double scale = 1.5;
	for(xpos=intmin;xpos<intmax;xpos+=step)
	{
		amp = XSFF(xpos,s,f1,f2,hl,hh,whl,whh,ghl,ghh);
		if(amp>maxamp) maxamp=amp;
	}
	maxamp = maxamp*scale;
	
	//MC begins here
	long nstep = 10000,i=0,ng=0,nt=0;
	double lastintegral=0;
	double integral=0;
	int randp = 8;
	double desiredprecision = pow(10.,-4);
	
	for(i=0;i<nstep;i++)
	{
		xpos = getrand(randp)*xscale+intmin;
		amppos = getrand(randp)*maxamp;
		amp = XSFF(xpos,s,f1,f2,hl,hh,whl,whh,ghl,ghh);
		if (amppos<amp) ng++;
		nt++;
		if(i==nstep-1)
		{
			integral = ng*1./(nt*1.);
			if((integral-lastintegral)/integral > desiredprecision)
			{
				i=0;
				lastintegral=integral;
			}
		}
		if(amp>maxamp)
		{
			maxamp = amp*scale;
			i=0; ng=0; nt=0; lastintegral=0;
		}
	}
	integral = ng*1./(nt*1.)*xscale*maxamp;
	return integral;
}

// FUNCTION MCSVV
// ==============

// evaluates the integral for the SVV decays
double MCSVV(double s, double w, double hl, double hh, double whl, double whh, double ghl, double ghh)
{
	double xpos, amppos, amp, maxamp, intmax, intmin, xscale;
	
	intmax = pow(1.-sqrt(s),2);
	intmin = 4.*w;
	
	xscale = intmax-intmin;
	maxamp = 0;
	
	//searches for a good initial cap 
	xpos = intmin;
	double step = xscale/1000;
	double scale = 1.5;
	for(xpos=intmin;xpos<intmax;xpos+=step)
	{
		amp = XSVV(xpos,s,w,hl,hh,whl,whh,ghl,ghh);
		if(amp>maxamp) maxamp=amp;
	}
	maxamp = maxamp*scale;
	
	//MC begins here
	long nstep = 10000,i=0,ng=0,nt=0;
	double lastintegral=0;
	double integral=0;
	int randp = 8;
	double desiredprecision = pow(10.,-4);
	
	for(i=0;i<nstep;i++)
	{
		xpos = getrand(randp)*xscale+intmin;
		amppos = getrand(randp)*maxamp;
		amp = XSVV(xpos,s,w,hl,hh,whl,whh,ghl,ghh);
		if (amppos<amp) ng++;
		nt++;
		if(i==nstep-1)
		{
			integral = ng*1./(nt*1.);
			if((integral-lastintegral)/integral > desiredprecision)
			{
				i=0;
				lastintegral=integral;
			}
		}
		if(amp>maxamp)
		{
			maxamp = amp*scale;
			i=0; ng=0; nt=0; lastintegral=0;
		}
	}
	integral = ng*1./(nt*1.)*xscale*maxamp;
	return integral;
}
// FUNCTION MCSVV2
// ===============

// evaluates the integral for the SVV decays
double MCSVV2(double s, double w, double h, double wh, double gh)
{
	double xpos, amppos, amp, maxamp, intmax, intmin, xscale;
	
	intmax = pow(1.-sqrt(w),2);
	intmin = pow(sqrt(w)+sqrt(s),2);
	
	xscale = intmax-intmin;
	maxamp = 0;
	
	//searches for a good initial cap 
	xpos = intmin;
	double step = xscale/1000;
	double scale = 1.5;
	for(xpos=intmin;xpos<intmax;xpos+=step)
	{
		amp = XSVV2(xpos,s,w,h,wh,gh);
		if(amp>maxamp) maxamp=amp;
	}
	maxamp = maxamp*scale;
	
	//MC begins here
	long nstep = 10000,i=0,ng=0,nt=0;
	double lastintegral=0;
	double integral=0;
	int randp = 8;
	double desiredprecision = pow(10.,-4);
	
	for(i=0;i<nstep;i++)
	{
		xpos = getrand(randp)*xscale+intmin;
		amppos = getrand(randp)*maxamp;
		amp = XSVV2(xpos,s,w,h,wh,gh);
		if (amppos<amp) ng++;
		nt++;
		if(i==nstep-1)
		{
			integral = ng*1./(nt*1.);
			if((integral-lastintegral)/integral > desiredprecision)
			{
				i=0;
				lastintegral=integral;
			}
		}
		if(amp>maxamp)
		{
			maxamp = amp*scale;
			i=0; ng=0; nt=0; lastintegral=0;
		}
	}
	integral = ng*1./(nt*1.)*xscale*maxamp;
	return integral;
}

// FUNCTION WSVVV
// ==============

// Three body decay - assumes that the A0 decays to a Z and an off-shell h or H
// which suseqently decays into WW or ZZ (d=1 or d=3 respectively)
// or the the A+- Decays to a W+- and then WW or ZZ (d=2)
// parameters are as follows:
// WWZ: d=1, MS=M_A0, MZ=Mz, MV=Mw, COUPL=g_{hAZ}*g_{hWW}, COUPL=g_{HAZ}*g_{HWW}, MSOL,WSOL=Mh,Wh, MSOH,WSOH=MH,WH
// ZZZ: d=3, MS=M_A0, MZ=Mz, MV=Mz, COUPL=g_{hAZ}*g_{hZZ}, COUPL=g_{HAZ}*g_{HZZ}, MSOL,WSOL=Mh,Wh, MSOH,WSOH=MH,WH
// WZZ: d=2, MS=M_A+-, MZ=Mw, MV=Mz, COUPL=g_{hA+W}*g_{hZZ}, COUPL=g_{HA+W}*g_{HZZ}, MSOL,WSOL=Mh,Wh, MSOH,WSOH=MH,WH
// WWW: d=2, MS=M_A+-, MZ=Mw, MV=Mw, COUPL=g_{hA+W}*g_{hZZ}, COUPL=g_{HA+W}*g_{HZZ}, MSOL,WSOL=Mh,Wh, MSOH,WSOH=MH,WH

double WSVVV(double d, double MS, double MZ, double MV, double MSOL, double WSOL, double MSOH, double WSOH, double COUPL, double COUPH)
{
	//conditions to quit
	if (MS>MSOL+MZ) return 0;
	if (MS<2.*MV+MZ) return 0;
	if (d!=1. && d!=2. && d!=3.) return 0;
	
	double coeff = MS/(256.*pow(PI,3));	
	double part1=0.,part2=0.;
	double z,w,hl,hh,whl,whh,ghl,ghh;
	z = pow(MZ/MS,2);
	w = pow(MV/MS,2);
	hl = pow(MSOL/MS,2);
	hh = pow(MSOH/MS,2);
	whl = pow(WSOL/MS,2);
	whh = pow(WSOH/MS,2);
	ghl = COUPL/MS;
	ghh = COUPH/MS;
	
	part1 = coeff*MCWWZ(z,w,hl,hh,whl,whh,ghl,ghh);
	//A0->WWZ
	if (d==1) return part1;
	//A+->WZZ
	if (d==2 && w!=z) return part1;
	
	part2 = coeff*MCZZZ(z,w,hl,hh,whl,whh,ghl,ghh);
	//A+->WWW
	if (d==2) return part1+part2/2.;
	//A+->ZZZ
	if (d==3) return (part1+part2)/2.;
	
	return 0.;
}

// FUNCTION WSVFF
// ==============

// Three body S->VS->VFF decay 
// Note: Scalar and pseudoscalars do not mix
// For PScalar set COUPH = 0;  vector or axial(0 or 1)

double WSVFF(double MS, double MV, double MF1, double MF2, double COLOR, double AV, double MSOL, double WSOL, double COUPL, double MSOH, double WSOH, double COUPH)
{
	//conditions to quit
	if (MS>MSOL+MV) return 0;
	if (MS<MF1+MF2+MV) return 0;
	if (AV!=0 && AV!=1) return 0;
	
	double coeff = COLOR*MS/(128.*pow(PI,3));	
	double integral=0.;
	double v,f1,f2,hl,hh,whl,whh,ghl,ghh;
	v = pow(MV/MS,2);
	f1 = MF1/MS;
	f2 = MF2/MS;
	hl = pow(MSOL/MS,2);
	hh = pow(MSOH/MS,2);
	whl = pow(WSOL/MS,2);
	whh = pow(WSOH/MS,2);
	ghl = COUPL;
	ghh = COUPH;
	AV = 1.-(2.*AV);
	
	integral = coeff*MCVFF(v,f1,AV*f2,hl,hh,whl,whh,ghl,ghh);
	
	return integral;
}

// FUNCTION WSSFF
// ==============
// Note: Code still in beta

// Three body S->SS->SFF decay 
// Note: Scalar and pseudoscalars do not mix
// For PScalar set COUPH = 0;  vector or axial(0 or 1)

double WSSFF(double MS, double MS1, double MF1, double MF2, double COLOR, double AV, double MSOL, double WSOL, double COUPL, double MSOH, double WSOH, double COUPH)
{
	//conditions to quit
	if (MS>MSOL+MS1) return 0;
	if (MS<MF1+MF2+MS1) return 0;
	if (AV!=0 && AV!=1) return 0;
	
	double coeff = COLOR*MS/(128.*pow(PI,3));	
	double integral=0.;
	double s,f1,f2,hl,hh,whl,whh,ghl,ghh;
	s = pow(MS1/MS,2);
	f1 = MF1/MS;
	f2 = MF2/MS;
	hl = pow(MSOL/MS,2);
	hh = pow(MSOH/MS,2);
	whl = pow(WSOL/MS,2);
	whh = pow(WSOH/MS,2);
	ghl = COUPL/MS;
	ghh = COUPH/MS;
	AV = 1.-(2.*AV);
	
	integral = coeff*MCSFF(s,f1,AV*f2,hl,hh,whl,whh,ghl,ghh);
	
	return integral;
}

// FUNCTION WSSVV
// ==============
// Note: Code still in beta

// Three body S->SS->SVV decay 

double WSSVV(double d, double MS, double MS1, double MV, double MSOL, double WSOL, double COUPL, double MSOH, double WSOH, double COUPH)
{
	//conditions to quit
	if (MS>MSOL+MS1) return 0;
	if (MS<MS1+2.*MV) return 0;
	
	double coeff = MS/(256.*pow(PI,3));	
	double integral=0.;
	double s,w,hl,hh,whl,whh,ghl,ghh;
	s = pow(MS1/MS,2);
	w = pow(MV/MS,2);
	hl = pow(MSOL/MS,2);
	hh = pow(MSOH/MS,2);
	whl = pow(WSOL/MS,2);
	whh = pow(WSOH/MS,2);
	ghl = COUPL/pow(MS,2);
	ghh = COUPH/pow(MS,2);
	
	integral = coeff*MCSVV(s,w,hl,hh,whl,whh,ghl,ghh)/d;
	
	return integral;
}
// FUNCTION WSSVV2
// ==============
// Note: Code still in beta

// Three body S->SV->SVV decay 
// W-> d=1, Z-> d=2

double WSSVV2(double d, double MS, double MS1, double MV, double MSO, double WSO, double COUP)
{
	//conditions to quit
	if (MS>MSO+MS1) return 0;
	if (MS<MS1+2.*MV) return 0;
	
	double coeff = MS/(256.*pow(PI,3));	
	double integral=0.;
	double s,w,h,wh,gh;
	s = pow(MS1/MS,2);
	w = pow(MV/MS,2);
	h = pow(MSO/MS,2);
	wh = pow(WSO/MS,2);
	gh = COUP;
	
	integral = coeff*MCSVV2(s,w,h,wh,gh);
	
	if(d==1) return integral;
	
	//integral += coeff*MCSVV2INT(s,w,h,wh,gh)/d;
	 return integral;
}


// FUNCTION SMTOPWID (from MGII)
// =============================

// c*************************************************************************
// c     THE TOTAL WEAK DECAY WIDTH OF THE TOP QUARK, INCLUDING
// c     THE EFFECTS OF BOTTOM MASS AND, IF IGW=1,  A FINITE W WIDTH.
// c     From James Stirling 6-10-94
// c
// c     RMT=TOP MASS
// c     RMW=W   MASS
// c     RMB=B   MASS
// c     RGW=W   WIDTH
// c     GW =WEAK COUPLING
// c
// c     RGT=TOP WIDTH
// c
// c*************************************************************************

double smtopwid(double RMT,double RMW, double RMB, double RGW, double GW) {
	
	double XW,XB,RGT;
	double complex XGW,RM,OM,Y1,Y0,Z1,Z0,D0,D1,A4,A3,A2,A1,B0,B1,B2,RINT,XGW4;
	
	XGW=GW/2./csqrt(2.);
	XB=RMB/RMT;
	XW=RMW/RMT;
	RM=cpow(XB,2.);
	
	OM=1.+RM-(pow(RMW,2.)+RMW*RGW*I)/pow(RMT,2.);
	Y1=OM+csqrt(OM*OM-4.*RM);
	Y0=OM-csqrt(OM*OM-4.*RM);
	Z1=2.;
	Z0=2.*csqrt(RM);
	
	D0=(-cpow(Y0,8.)+3.*cpow(Y0,7.)*RM+3.*cpow(Y0,7.)-8.*cpow(Y0,6.)*RM-12.*cpow(Y0,5.)*cpow(RM,
																							 2.)-12.*cpow(Y0,5.)*RM+96.*cpow(Y0,4.)*cpow(RM,2.)-48.*cpow(Y0,3.)*cpow(RM,3.)-48.*cpow(Y0,3.)*
		cpow(RM,2.)-128.*cpow(Y0,2.)*cpow(RM,3.)+192.*Y0*cpow(RM,4.)+192.*Y0*cpow(RM,3.)-256.*
		cpow(RM,4.))/(24.*cpow(Y0,4.)*(Y1-Y0));
	D1=(-cpow(Y1,8.)+3.*cpow(Y1,7.)*RM+3.*cpow(Y1,7.)-8.*cpow(Y1,6.)*RM-12.*cpow(Y1,5.)*cpow(RM,
																							 2.)-12.*cpow(Y1,5.)*RM+96.*cpow(Y1,4.)*cpow(RM,2.)-48.*cpow(Y1,3.)*cpow(RM,3.)-48.*cpow(Y1,3.)*
		cpow(RM,2.)-128.*cpow(Y1,2.)*cpow(RM,3.)+192.*Y1*cpow(RM,4.)+192.*Y1*cpow(RM,3.)-256.*
		cpow(RM,4.))/(24.*cpow(Y1,4.)*(Y1-Y0));
	A4=(32.*cpow(RM,4.)*(Y1-Y0))/(3.*Y1*Y0*(Y1-Y0));
	A3=(8.*cpow(RM,3.)*(-3.*cpow(Y1,2.)*Y0*RM-3.*cpow(Y1,2.)*Y0+4.*cpow(Y1,2.)*RM+3.*
						Y1*cpow(Y0,2.)*RM+3.*Y1*cpow(Y0,2.)-4.*cpow(Y0,2.)*RM))/(3.*cpow(Y1,2.)*cpow(Y0,2.)*(Y1
																											 -Y0));
	A2=(8.*cpow(RM,3.)*(2.*cpow(Y1,3.)*cpow(Y0,2.)-3.*cpow(Y1,3.)*Y0*RM-3.*cpow(Y1,3.)*Y0+4.
						*cpow(Y1,3.)*RM-2.*cpow(Y1,2.)*cpow(Y0,3.)+3.*Y1*cpow(Y0,3.)*RM+3.*Y1*cpow(Y0,3.)-4.*cpow(Y0
																												  ,3.)*RM))/(3.*cpow(Y1,3.)*cpow(Y0,3.)*(Y1-Y0));
	A1=(2.*cpow(RM,2.)*(3.*cpow(Y1,4.)*pow(Y0,3.)*RM+3.*cpow(Y1,4.)*cpow(Y0,3.)+8.*cpow(Y1,4.)*cpow(Y0
																									,2.)*RM-12.*cpow(Y1,4.)*Y0*cpow(RM,2.)-12.*cpow(Y1,4.)*Y0*RM+16.*cpow(Y1,4.)*cpow(RM,2.)
						-3.*cpow(Y1,3.)*cpow(Y0,4.)*RM-3.*cpow(Y1,3.)*cpow(Y0,4.)-8.*cpow(Y1,2.)*cpow(Y0,4.)*RM+12.*
						Y1*cpow(Y0,4.)*cpow(RM,2.)+12.*Y1*cpow(Y0,4.)*RM-16.*cpow(Y0,4.)*cpow(RM,2.)))/(3.*cpow(Y1,
																												4.)*cpow(Y0,4.)*(Y1-Y0));
	B0=(cpow(Y1,3.)-3.*cpow(Y1,2.)*RM-3.*cpow(Y1,2.)+8.*Y1*RM-cpow(Y0,3.)+3.*cpow(Y0,2.)*RM+
		3.*cpow(Y0,2.)-8.*Y0*RM)/(24.*(Y1-Y0));
	B1=(Y1+Y0-3.*RM-3.)/24.;
	B2=1./24.;
	
	RINT=D0*clog((Z1-Y0)/(Z0-Y0))
	-D1*clog((Y1-Z1)/(Y1-Z0))
	-A4/3.*(1./pow(Z1,3.)-1./pow(Z0,3.))
	-A3/2.*(1./pow(Z1,2.)-1./pow(Z0,2.))
	-A2   *(1./Z1   -1./Z0   )
	+A1*clog(Z1/Z0)
	+B0   *(Z1   -Z0   )
	+B1/2.*(pow(Z1,2.)-pow(Z0,2.))
	+B2/3.*(pow(Z1,3.)-pow(Z0,3.));
	
	XGW4=pow(XGW,4.);
	
	// TOTAL WIDTH INCLUDES FLAVOUR & COLOUR FACTORS
	RGT=pow(RMT,3.)/(RMW*RGW)*XGW4/(8.*pow(PI,3.))*cimag(RINT);
	RGT=9.*RGT;
	
	return RGT;
	
}

// FUNCTION HEAVY (from MGII)
// ==========================

double heavy(double x, double y, double z) { 
	
	return ( 1. - 0.5*(pow(y,2.)+pow(z,2.)) - 0.5*pow((pow(y,2.)-pow(z,2.)),2.)
			+ 3.*y*z*(pow(x,2.) - 1.)/(pow(x,2.) + 1.)        )
	* sqrt( pow((1.-pow(y,2.)-pow(z,2.)),2.) - 4. * pow(y,2.) * pow(z,2.) );
}

// FUNCTION SMZWID (from MGII)
// ===========================

double smzwid(double zmass, double lmass, double cmass, double bmass, double ee, double cw, double sw) {
	
	double ez  = ee/(sw*cw);
	double ey  = ee*(sw/cw);
	double sin2w=pow(sw,2.);
	
	double gzn1 = -ez*0.5;
	double gzn2 = 0.;
	double gzl1 = -ez*(-0.5 + sin2w);
	double gzl2 = -ey;
	double gzu1 = -ez*( 0.5 - sin2w*2./3.);
	double gzu2 = ey*2./3.;
	double gzd1 = -ez*(-0.5 + sin2w/3.);
	double gzd2 = -ey/3.;
	
	double decz = zmass / ( 24. * PI );
	double w_z_nn = decz * ( pow(gzn1,2.) +pow(gzn2,2.) );
	double w_z_ll = decz * ( pow(gzl1,2.) +pow(gzl2,2.) );
	decz = decz * 3.;
	double w_z_uu = decz * ( pow(gzu1,2.) +pow(gzu2,2.) );
	double w_z_dd = decz * ( pow(gzd1,2.) +pow(gzd2,2.) );
	double dum=(gzl2+gzl1)/(gzl2-gzl1);
	double w_z_tau = w_z_ll * heavy( dum, lmass/zmass, lmass/zmass );
	dum=(gzu2+gzu1)/(gzu2-gzu1);
	double w_z_cc = w_z_uu *  heavy( dum, cmass/zmass, cmass/zmass );
	dum=(gzd2+gzd1)/(gzd2-gzd1);
	double w_z_bb = w_z_dd *  heavy( dum, bmass/zmass, bmass/zmass );
	
	double zwidth =   3.*w_z_nn + 2.*w_z_ll + w_z_tau + 2.*w_z_dd + w_z_uu + w_z_cc + w_z_bb;
	
	return zwidth;
}

// FUNCTION SMWWID (from MGII)
// ===========================

double smwwid(double wmass, double lmass, double cmass, double ee, double cw, double sw) {
	
	double ez  = ee/(sw*cw);
	double ey  = ee*(sw/cw);
	double sin2w=pow(sw,2.);
	double gwf1 = -ee/sqrt(2.*sin2w);
	double gwf2 = 0;
	
	double decw = wmass / ( 24.0 * PI );
	double w_w_nl = decw * ( pow(gwf1,2.) + pow(gwf2,2.) );
	double dum = (gwf2+gwf1)/(gwf2-gwf1);
	double w_w_tau = w_w_nl * heavy( dum, lmass/wmass, 0. );
	double w_w_ud = w_w_nl * 3.;
	double w_w_cs = w_w_ud * heavy( dum, cmass/wmass, 0. );
	
	return 2.*w_w_nl + w_w_tau + w_w_ud + w_w_cs;
}



// FUNCTION MAIN
// =============

int main(int argc, char *argv[]) 
{	
	// FILE DECLARATIONS
	FILE* data;
	FILE* output;
	FILE* log;
	
	// INPUT VARIABLES WITH DEFAULT VALUES
	double Y1U=0,Y1C=0,Y1T=0,Y1D=0,Y1S=0,Y1B=0,Y1E=0,Y1Mu=0,Y1Ta=0;
	double Vud=1,Vus=0,Vub=0,Vcd=0,Vcs=1,Vcb=0,Vtd=0,Vts=0,Vtb=1;
	double MSMUMASS=0.106, MSTAMASS=1.777, MSSMASS=0.100, MSCMASS=1.2, MSBMASS=4.2, MSTMASS=174.;
	double MUMASS=0.106, TAMASS=1.777, SMASS=0.100, CMASS=1.2, BMASS=4.2, TMASS=174.;
	double alpha=1./128.91, G_Fermi=1.16637e-5, alpha_s=0.1172, ZMASS=91.1876;
	
	//Numbers 
	double zero = 0., one = 1., two = 2., three = 3., four = 4., eight = 8., half = 0.5;
	
	// 2HDM Params
	double MH1=120.,MH2=300.,MH3=250.,tb=1.,sa=0.8,sgamma=0.,MHC=0;
	double ca=0.,cb=0.,sb=1.,sab=1.,cab=1.,s2a=1.,s2b=1.;
	double sag=0.,cag=1.,cbg=0.,sbg=1.;
	
	//Warning params (defaults written)
	double LWARN = 80.;
	double TWARN = 6.;
	
	
	// ##############################################################################
	// ##
	// ##  INPUT READING
	// ##
	// ##############################################################################
	
	// OPEN FILES
	data=fopen(argv[1],"r");
	if (data==0) 
	{
		fprintf(stderr,"Cannot open LHA file %s for reading!\n",argv[1]);
		return -1;
	}
	output=fopen(argv[2],"w");
	if (output==0) 
	{
		fprintf(stderr,"Cannot open output file %s for writing!\n",argv[2]);
		return -1;
	}
	if(argc>=4) log=fopen(argv[3],"w"); else log=fopen("2hdm4tccalc.log","w");
	if (log==0) 
	{
		fprintf(stderr,"Cannot open log file for writing!\n");
		return -1;
	}
	
	// HEADER
	fprintf(output,"# 2HDM4TC Parameters\n");
    fprintf(output,"# Spectrum and decay, produced by 2HDM4TCCalc\n");
    fprintf(output,"# ********************************************\n");
	
	
	// Reading routine for LHA file. Based on M. Herquet's TwoHiggsCalc.
	
	// Dummy variables for reading
	char cbuf[MAXLG]; // Line buffer
	int i=0,i1=0,i2=0,i3=0; // Dum integers
	double val=0.; // Read value
	int blt=0; // Block tag
	int basis=0; // Basis tag
	
	// Read the input file
	while(fgets(cbuf,MAXLG,data)) 
	{
		// Everything to lowercase
		for(i=0;cbuf[i];i++)
			if(isalpha(cbuf[i]) && isupper(cbuf[i])) cbuf[i]=tolower(cbuf[i]);
		// Identify the block
		if(strncmp(cbuf,"block",5)==0) {
			char blnm[MAXBLNM];
			sscanf(cbuf+5,"%s",blnm);
			if(strcmp(blnm,"sminputs")==0) blt=1;
			else if(strcmp(blnm,"tcparams")==0) blt=2;
			else if(strcmp(blnm,"mgyukawa")==0) blt=3;
			else if(strcmp(blnm,"mass")==0) blt=4;
			else if(strcmp(blnm,"mgckm")==0) blt=11;
			else blt=0;
			continue;
		}
		// Ignore comments (or in general line starting without ' ') and unknown blocks
		if(cbuf[0]!=' ' || blt==0) continue;
		
		// Scan data line according to current block style and reject line with spurious format
		// Note: it is safe to write it like this because C does not evalute following expressions in a && list if the first one is false
		if(blt<10 && sscanf(cbuf,"%d%lf",&i1,&val)!=2) continue;
		if(blt>10 && sscanf(cbuf,"%d%d%lf",&i1,&i2,&val)!=3) continue;
		
		// Switch on the model tag
		switch(blt) 
		{
			// Read SM variables
			case 1: 
				if(i1==1) fprintf(log,"read alpha =% 16.8e\n",alpha=1/val);
				else if(i1==2) fprintf(log,"read G_Fermi =% 16.8e\n",G_Fermi=val);
				else if(i1==3) fprintf(log,"read alpha_s =% 16.8e\n",alpha_s=val);
				else if(i1==4) fprintf(log,"read ZMASS =% 16.8e\n",ZMASS=val);
				break;
			//Read Two Higgs Doublet Model Parameters
			case 2: 
				if(i1==1) fprintf(log,"read H1MASS =% 16.8e\n",MH1=val);
				else if(i1==2) fprintf(log,"read H2MASS =% 16.8e\n",MH2=val);
				else if(i1==3) fprintf(log,"read H3MASS =% 16.8e\n",MH3=val);
				else if(i1==4) fprintf(log,"read Tan(beta) =% 16.8e\n",tb=val);
				else if(i1==5) fprintf(log,"read Sin(alpha) =% 16.8e\n",sa=val);
				else if(i1==6) fprintf(log,"read Sin(gamma) =% 16.8e\n",sgamma=val);
				else if(i1==7) fprintf(log,"read lambda_i warning =% 16.8e\n", LWARN=val);
				else if(i1==8) fprintf(log,"read lambda_t warning =% 16.8e\n", TWARN=val);
				break;
			//Read Masses
			case 3: 
				if(i1==3) fprintf(log,"read MSSMASS =% 16.8e\n",MSSMASS=val);
				else if(i1==5) fprintf(log,"read MSBMASS =% 16.8e\n",MSBMASS=val);
				else if(i1==4) fprintf(log,"read MSCMASS =% 16.8e\n",MSCMASS=val);
				else if(i1==6) fprintf(log,"read MSTMASS =% 16.8e\n",MSTMASS=val);
				else if(i1==13) fprintf(log,"read MSMUMASS =% 16.8e\n",MSMUMASS=val);
				else if(i1==15) fprintf(log,"read MSTAMASS =% 16.8e\n",MSTAMASS=val); 
				break;
			case 4: 
				if(i1==3) fprintf(log,"read SMASS =% 16.8e\n",SMASS=val);
				else if(i1==5) fprintf(log,"read BMASS =% 16.8e\n",BMASS=val);
				else if(i1==4) fprintf(log,"read CMASS =% 16.8e\n",CMASS=val);
				else if(i1==6) fprintf(log,"read TMASS =% 16.8e\n",TMASS=val);
				else if(i1==13) fprintf(log,"read MUMASS=% 16.8e\n",MUMASS=val);
				else if(i1==15) fprintf(log,"read TAMASS =% 16.8e\n",TAMASS=val);
				break;
			//Read CKM
			case 11:	
				if(i1==1 && i2==1) fprintf(log,"read Vud =% 16.8e\n",Vud=val); 
				/*
				else if(i1==1 && i2==2) fprintf(log,"read Vus =% 16.8e\n",Vus=val); 
				else if(i1==1 && i2==3) fprintf(log,"read Vub =% 16.8e\n",Vub=val); 
				else if(i1==2 && i2==1) fprintf(log,"read Vcd =% 16.8e\n",Vcd=val); 
				else if(i1==2 && i2==2) fprintf(log,"read Vcs =% 16.8e\n",Vcs=val); 
				else if(i1==2 && i2==3) fprintf(log,"read Vcb =% 16.8e\n",Vcb=val); 
				else if(i1==3 && i2==1) fprintf(log,"read Vtd =% 16.8e\n",Vtd=val); 
				else if(i1==3 && i2==2) fprintf(log,"read Vts =% 16.8e\n",Vts=val); 
				else if(i1==3 && i2==3) fprintf(log,"read Vtb =% 16.8e\n",Vtb=val); */
				break;
		}
	}
	
	// Compute secondary SM variables
	double ee=sqrt(alpha*4.*PI), sw=sin(asin(2.*sqrt(PI*alpha/(sqrt(2.)*G_Fermi*pow(ZMASS,2.))))/2.); 
	double cw=sqrt(1.-pow(sw,2.)), WMASS=ZMASS*cw,v=2.*WMASS/(ee/sw), sin2w = pow(sw,2);
	double sc2 = sin2w*(one-sin2w),ee2=ee*ee, ez=ee/(sw*cw), gw=ee/sw;
	
	//same by custodial SU(2)
	MHC = MH3;

	//mixing angles
	double beta2HDM = atan(tb);
	double alpha2HDM = asin(sa);
	double gamma2HDM = asin(sgamma);
	sb = sin(beta2HDM);
	cb = cos(beta2HDM);
	sa = sin(alpha2HDM);
	ca = cos(alpha2HDM);
	sab = sin(alpha2HDM-beta2HDM);
	cab = cos(alpha2HDM-beta2HDM);
	s2a = sin(2*alpha2HDM);
	s2b = sin(2*beta2HDM);
	
	sag = sin(alpha2HDM+gamma2HDM);
	cag = cos(alpha2HDM+gamma2HDM);
	sbg = sin(beta2HDM+gamma2HDM);
	cbg = cos(beta2HDM+gamma2HDM);
	
	// Recopy input
	fprintf(output, "Block SMINPUTS      # Standard Model inputs\n");
	fprintf(output,  "     1      % 16.13E   # alpha em(MZ)(-1) SM MSbar\n",1/alpha);
	fprintf(output,  "     2      % 16.13E   # G Fermi\n",G_Fermi);
	fprintf(output,  "     3      % 16.13E   # alpha s(MZ) SM MSbar\n",alpha_s);
	fprintf(output,  "     4      % 16.13E   # MZ(Pole)\n",ZMASS);  
	fprintf(output, "Block MGCKM   # CKM Matrix\n");
	fprintf(output,  "    1     1      % 16.13E   # Vud\n",Vud); 
	/* currently not supported
	fprintf(output,  "    1     2      % 16.13E   # Vus\n",Vus);
	fprintf(output,  "    1     3      % 16.13E   # Vub\n",Vub);
	fprintf(output,  "    2     1      % 16.13E   # Vcd\n",Vcd);
	fprintf(output,  "    2     2      % 16.13E   # Vcs\n",Vcs);
	fprintf(output,  "    2     3      % 16.13E   # Vcb\n",Vcb);
	fprintf(output,  "    3     1      % 16.13E   # Vtd\n",Vtd);
	fprintf(output,  "    3     2      % 16.13E   # Vts\n",Vts);
	fprintf(output,  "    3     3      % 16.13E   # Vtb\n",Vtb); */
	fprintf(output, "Block MGYUKAWA   # Yukawa masses\n");
	fprintf(output,  "#    PDG          YMASS\n");
	fprintf(output,  "     3      % 16.13E   # Ms\n",MSSMASS);
	fprintf(output,  "     4      % 16.13E   # Mc\n",MSCMASS);
	fprintf(output,  "     5      % 16.13E   # Mb\n",MSBMASS);
	fprintf(output,  "     6      % 16.13E   # Mt\n",MSTMASS);
	fprintf(output,  "     13     % 16.13E   # Mmu\n",MSMUMASS);
	fprintf(output,  "     15     % 16.13E   # Mta\n",MSTAMASS);
	
	// Add Qnumbers block
	
	fprintf(output,  "#===========================================================\n");
	fprintf(output,  "# QUANTUM NUMBERS OF NEW STATE(S) (NON SM PDG CODE) IF ANY\n");
	fprintf(output,  "# (see below for masses and decay tables)\n");
	fprintf(output,  "# These blocks are automatically created by the MadGraph\n");
	fprintf(output,  "# qnumbers.pl script from the particles.dat model file\n");
	fprintf(output,  "#===========================================================\n");
	fprintf(output,  "# END of QNUMBERS blocks\n");
	fprintf(output,  "#===========================================================\n");
	
	//Masses
	fprintf(output, "Block MASS # Mass spectrum (kinematic masses)\n");
	fprintf(output, "# PDG Code                mass    particle\n");
	
	// SM Particles first
	// ------------------
	fprintf(output, "# SM PARTICLES\n");
	fprintf(output, "         3    %16.13E    # Ms\n",SMASS);
	fprintf(output, "         4    %16.13E    # Mc\n",CMASS);
	fprintf(output, "         5    %16.13E    # Mb\n",BMASS);
	fprintf(output, "         6    %16.13E    # Mt\n",TMASS);
	fprintf(output, "        13    %16.13E    # Mmu\n",MUMASS);
	fprintf(output, "        15    %16.13E    # Mtau\n",TAMASS);
	fprintf(output, "        23    %16.13E    # Mz\n",ZMASS);
	fprintf(output, "        24    %16.13E    # Mw\n",WMASS);
	fprintf(output, "# 2HDM ADD SCALARS\n");
	fprintf(output, "        25    %16.13E    # h mass\n",MH1);
	fprintf(output, "        35    %16.13E    # H mass\n",MH2);
	fprintf(output, "        36    %16.13E    # A mass\n",MH3);
	fprintf(output, "        37    %16.13E    # Charged Higgs\n",MHC);
	
	
	// Check fundamental lagrangian parameters / print warnings
	double lambda1 = two*(pow(ca*MH1,2)+pow(sa*MH2,2))/pow(v*sb,2);
	double lambda2 = two*(pow(ca*MH1,2)+pow(sa*MH2,2))/pow(v*cb,2);
	double lambda3 = fabs(sa*ca*(MH2*MH2-MH1*MH1)/(v*v*sb*cb)+two*pow(MH3/v,2));
	double lambda4 = two*pow(MH3/v,2);
	double lambdat = one/(sbg);
	int hadissue = 0;
	
    if(lambda1>LWARN)   
	{
		fprintf(log,"WARNING:  lambda_1 =% 11.4e",lambda1);
		fprintf(log, " which is greater than %4.2f\n",LWARN);
		fprintf(stdout,"WARNING:  lambda_1 =% 11.4e",lambda1);
		fprintf(stdout, " which is greater than %4.2f\n",LWARN);
		hadissue = 1;
	}
	if(lambda2>LWARN)   
	{
		fprintf(log,"WARNING:  lambda_2 =% 11.4e",lambda2);
		fprintf(log, " which is greater than %4.2f\n",LWARN);
		fprintf(stdout,"WARNING:  lambda_2 =% 11.4e",lambda2);
		fprintf(stdout, " which is greater than %4.2f\n",LWARN);
		hadissue = 1;
	}
	if(lambda3>LWARN)   
	{
		fprintf(log,"WARNING:  lambda_3 =% 11.4e",lambda3);
		fprintf(log, " which is greater than %4.2f\n",LWARN);
		fprintf(stdout,"WARNING:  lambda_3 =% 11.4e",lambda3);
		fprintf(stdout, " which is greater than %4.2f\n",LWARN);
		hadissue = 1;
	}
	if(lambda4>LWARN)   
	{
		fprintf(log,"WARNING:  lambda_4 =% 11.4e",lambda4);
		fprintf(log, " which is greater than %4.2f\n",LWARN);
		fprintf(stdout,"WARNING:  lambda_4 =% 11.4e",lambda4);
		fprintf(stdout, " which is greater than %4.2f\n",LWARN);
		hadissue = 1;
	}
	if(hadissue==1) fprintf(stdout,"Loop effects overtake tree effects near %4.2f\n",16.*pow(PI,2));
	hadissue = 0;
	if(lambdat>TWARN)   
	{
		fprintf(log,"WARNING:  lambda_t =% 11.4e",lambdat);
		fprintf(log, "which is greater than %4.2f\n",TWARN);
		fprintf(stdout,"WARNING:  lambda_t =% 11.4e",lambdat);
		fprintf(stdout, " which is greater than %4.2f\n",TWARN);
		hadissue = 1;
	}
	if(hadissue==1) fprintf(stdout,"Loop effects overtake tree effects near %4.2f\n",4.*PI);
	
	
	// ##############################################################################
	// ##
	// ##  WIDTH COMPUTATION
	// ##
	// ##############################################################################
	
	//couplings
	
	//  VVS couplings
	double gwwh1  = ee2/sin2w*half*v*sab;
	double gwwh2  = ee2/sin2w*half*v*cab;
	double gzzh1  = ee2/sc2*half*v*sab;
	double gzzh2  = ee2/sc2*half*v*cab;
	double gwwh3  = zero;
	double gzzh3  = zero;
	
	//  VSS couplings
		  
	double gzh1h3 = -half*ez*cab;
	double gzh2h3 = -half*ez*sab;
	double gzh1h2 = zero;
	
	//   charged higgs	 
	double gzhchc = half*gw;
	double gwhch1 = half*gw*cab;
	double gwhch2 = half*gw*sab;
	double gwhch3 = half*gw;
	
	//  SSS couplings

	double gh112 = (two*pow(MH1,2) + pow(MH2,2))*cab*s2a/(v*s2b);
	double gh122 = (pow(MH1,2) + two*pow(MH2,2))*sab*s2a/(v*s2b);
	double gh133 = ( (pow(MH1,2) - two*pow(MH3,2))*cos(alpha2HDM-three*beta2HDM) +(three*pow(MH1,2) + two*pow(MH3,2))*cos(alpha2HDM+beta2HDM) )/(two*v*s2b);
	double gh233 = ( (pow(MH2,2) - two*pow(MH3,2))*sin(alpha2HDM-three*beta2HDM) +(three*pow(MH2,2) + two*pow(MH3,2))*sin(alpha2HDM+beta2HDM) )/(two*v*s2b);
	double gh1hmhp = ( (pow(MH1,2) - two*pow(MH3,2))*cos(alpha2HDM-three*beta2HDM) +(three*pow(MH1,2) + two*pow(MH3,2))*cos(alpha2HDM+beta2HDM) )/(two*v*s2b);
	double gh2hmhp = ( (pow(MH2,2) - two*pow(MH3,2))*sin(alpha2HDM-three*beta2HDM) +(three*pow(MH2,2) + two*pow(MH3,2))*sin(alpha2HDM+beta2HDM) )/(two*v*s2b);
	double gh3hmhp = zero;
	
	// SM widths
	fprintf(output,"#     PDG               Width\n");
	fprintf(output,"DECAY   6    % 16.13E      # top width\n",smtopwid(TMASS,WMASS,BMASS,smwwid(WMASS,TAMASS,CMASS,ee,cw,sw),ee/sw));
	fprintf(output,"DECAY  23    % 16.13E      # z width \n",smzwid(ZMASS,TAMASS,CMASS,BMASS,ee,cw,sw));
	fprintf(output,"DECAY  24    % 16.13E      # w width (SM decays only)\n",smwwid(WMASS,TAMASS,CMASS,ee,cw,sw));
	
	// H1 DECAY -> 2 particles
	double H1FF = cag/sbg;
	
	double WH1MUMU = WSFF(one,MH1,MUMASS,MUMASS,H1FF*MSMUMASS/v,zero);
	double WH1SS = WSFF(three,MH1,SMASS,SMASS,H1FF*MSSMASS/v,zero);
	double WH1TATA = WSFF(one,MH1,TAMASS,TAMASS,H1FF*MSTAMASS/v,zero);
	double WH1CC = WSFF(three,MH1,CMASS,CMASS,H1FF*MSCMASS/v,zero);
	double WH1BB = WSFF(three,MH1,BMASS,BMASS,H1FF*MSBMASS/v,zero);
	double WH1TT = WSFF(three,MH1,TMASS,TMASS,H1FF*MSTMASS/v,zero);
	
	double WH1WHM = WSVS(MH1,MHC,WMASS,gwhch1);
	double WH1WHP = WSVS(MH1,MHC,WMASS,gwhch1);
	double WH1ZH3 = WSVS(MH1,MH3,ZMASS,gzh1h3);
	
	double WH1WW = WSVV(two,MH1,WMASS,gwwh1,sw);
	double WH1ZZ = WSVV(one,MH1,ZMASS,gzzh1,sw);
	double WH1HMHP = WSSS(two,MH1,MHC,MHC,gh1hmhp);
	double WH1H3H3 = WSSS(one,MH1,MH3,MH3,gh133);
	
	double WH1=WH1MUMU+WH1SS+WH1TATA+WH1CC+WH1BB+WH1TT;
	WH1 += WH1WW+WH1ZZ+WH1WHM+WH1WHP+WH1HMHP+WH1H3H3;
	
	if (WH1 == 0.) WH1=0.001;
	
	
	// H2 DECAY -> 2 particles
	
	double H2FF = sag/sbg;
	
	double WH2MUMU = WSFF(one,MH2,MUMASS,MUMASS,H2FF*MSMUMASS/v,zero);
	double WH2SS = WSFF(three,MH2,SMASS,SMASS,H2FF*MSSMASS/v,zero);
	double WH2TATA = WSFF(one,MH2,TAMASS,TAMASS,H2FF*MSTAMASS/v,zero);
	double WH2CC = WSFF(three,MH2,CMASS,CMASS,H2FF*MSCMASS/v,zero);
	double WH2BB = WSFF(three,MH2,BMASS,BMASS,H2FF*MSBMASS/v,zero);
	double WH2TT = WSFF(three,MH2,TMASS,TMASS,H2FF*MSTMASS/v,zero);
	
	double WH2WHM = WSVS(MH2,MHC,WMASS,gwhch2);
	double WH2WHP = WSVS(MH2,MHC,WMASS,gwhch2);
	
	double WH2WW = WSVV(two,MH2,WMASS,gwwh2,sw);
	double WH2ZZ = WSVV(one,MH2,ZMASS,gzzh2,sw);
	double WH2HMHP = WSSS(two,MH2,MHC,MHC,gh2hmhp);
	double WH2H3H3 = WSSS(one,MH2,MHC,MHC,gh233);
	double WH2ZH1 = WSVS(MH2,MH1,ZMASS,gzh1h2);
	double WH2ZH3 = WSVS(MH2,MH3,ZMASS,gzh2h3);
	double WH2H1H1 = WSSS(one,MH2,MH1,MH1,gh112);
	
	double WH2=WH2MUMU+WH2SS+WH2TATA+WH2CC+WH2BB+WH2TT;
	WH2 += WH2WW+WH2ZZ+WH2HMHP+WH2H3H3+WH2ZH1+WH2ZH3+WH2H1H1+WH2WHM+WH2WHP;
	
	if (WH2 == 0.) WH2=0.001;
	
	// H3 DECAY -> 2 particles
	
	double H3FF = cbg/sbg;
	
	double WH3MUMU = WSFF(one,MH3,MUMASS,MUMASS,H3FF*MSMUMASS/v,one);
	double WH3SS = WSFF(three,MH3,SMASS,SMASS,H3FF*MSSMASS/v,one);
	double WH3TATA = WSFF(one,MH3,TAMASS,TAMASS,H3FF*MSTAMASS/v,one);
	double WH3CC = WSFF(three,MH3,CMASS,CMASS,H3FF*MSCMASS/v,one);
	double WH3BB = WSFF(three,MH3,BMASS,BMASS,H3FF*MSBMASS/v,one);
	double WH3TT = WSFF(three,MH3,TMASS,TMASS,H3FF*MSTMASS/v,one);
	
	double WH3WHM = WSVS(MH3,MHC,WMASS,gwhch3);
	double WH3WHP = WSVS(MH3,MHC,WMASS,gwhch3);
	
	double WH3WW = WSVV(two,MH3,WMASS,gwwh3,sw);
	double WH3ZZ = WSVV(one,MH3,ZMASS,gzzh3,sw);
	double WH3HMHP = WSSS(two,MH3,MHC,MHC,gh3hmhp);
	double WH3ZH1 = WSVS(MH3,MH1,ZMASS,gzh1h3);
	double WH3ZH2 = WSVS(MH3,MH2,ZMASS,gzh2h3);
	
	double WH3=WH3MUMU+WH3SS+WH3TATA+WH3CC+WH3BB+WH3TT;
	WH3 += WH3WW+WH3ZZ+WH3HMHP+WH3ZH1+WH3ZH2+WH3WHM+WH3WHP;
	
	if (WH3 == 0.) WH3=0.001;
	
	
	// H+ DECAY -> 2 particles
	
	double HCFF = cbg/sbg;
	
	double WHPVEE = WSFF(one,MHC,zero,zero,zero,one);
	double WHPVMMU = WSFF(one,MHC,zero,MUMASS,HCFF*MSMUMASS/v,one);
	double WHPVTTA = WSFF(one,MHC,zero,TAMASS,HCFF*MSTAMASS/v,one);
	
	double WHPUS = WSFF(three,MHC,zero,SMASS,zero,one);
	double WHPCD = WSFF(three,MHC,CMASS,zero,zero,one);
	double WHPCS = WSFF(three,MHC,CMASS,SMASS,HCFF*(MSSMASS+MSCMASS)/v,one);
	double WHPTB = WSFF(three,MHC,TMASS,BMASS,HCFF*(MSBMASS+MSTMASS)/v,one);
	
	double WHPWH1 = WSVS(MHC,MH1,WMASS,gwhch1);
	double WHPWH2 = WSVS(MHC,MH2,WMASS,gwhch2);
	double WHPWH3 = WSVS(MHC,MH3,WMASS,gwhch3);
	
	double WHC=WHPVEE+WHPVMMU+WHPVTTA;
	WHC +=WHPUS+WHPCD+WHPCS+WHPTB+WHPWH1+WHPWH2+WHPWH3;
	
	if (WHC == 0.) WHC=0.001;
	
	//3-Body Widths (in beta)
	//Note, the widths inserted here may be a bit smaller than the actual value
	double WH3WWZ = WSVVV(one, MH3, ZMASS, WMASS, MH1, WH1, MH2, WH2, gwwh1*gzh1h3, gwwh2*gzh2h3);
	double WH3ZZZ = WSVVV(three, MH3, ZMASS, ZMASS, MH1, WH1, MH2, WH2, gzzh1*gzh1h3, gzzh2*gzh2h3);
	double WHPWZZ = WSVVV(two, MHC, WMASS, ZMASS, MH1, WH1, MH2, WH2, gzzh1*gwhch1, gzzh2*gwhch2);
	double WHPWWW = WSVVV(two, MHC, WMASS, WMASS, MH1, WH1, MH2, WH2, gwwh1*gwhch1, gwwh2*gwhch2);
	WH3 += WH3WWZ+WH3ZZZ;
	WHC += WHPWZZ+WHPWWW;
	
	double WH3WTB = WSVFF(MH3, WMASS, TMASS, BMASS, three, 1., MHC, WHC, gwhch3*HCFF*(MSBMASS+MSTMASS)/v, 0., 0., 0.);
	double WH3WBT = WH3WTB;
	double WH3ZTT = WSVFF(MH3, ZMASS, TMASS, TMASS, three, 0., MH1, WH1, gzh1h3*H1FF*MSTMASS/v, MH2, WH2, gzh2h3*H2FF*MSTMASS/v);
	double WH3ZBB = WSVFF(MH3, ZMASS, BMASS, BMASS, three, 0., MH1, WH1, gzh1h3*H1FF*MSBMASS/v, MH2, WH2, gzh2h3*H2FF*MSBMASS/v);
	WH3 += WH3WBT+WH3WTB+WH3ZTT+WH3ZBB;
	
	double WHPWTT = WSVFF(MHC, WMASS, TMASS, TMASS, three, 1., MH3, WH3, gwhch3*H3FF*MSTMASS/v, 0., 0., 0.);
	WHPWTT += WSVFF(MHC, WMASS, TMASS, TMASS, three, 0., MH1, WH1, gwhch1*H1FF*MSTMASS/v, MH2, WH2, gwhch2*H2FF*MSTMASS/v);
	double WHPWBB = WSVFF(MHC, WMASS, BMASS, BMASS, three, 1., MH3, WH3, gwhch3*H3FF*MSBMASS/v, 0., 0., 0.);
	WHPWBB += WSVFF(MHC, WMASS, BMASS, BMASS, three, 0., MH1, WH1, gwhch1*H1FF*MSBMASS/v, MH2, WH2, gwhch2*H2FF*MSBMASS/v);
	WHC += WHPWTT+WHPWBB;
	
	double WHPZTB = WSVFF(MHC, ZMASS, TMASS, BMASS, three, 1., MHC, WHC, gzhchc*HCFF*(MSBMASS+MSTMASS)/v, 0., 0., 0.);
	WHC += WHPZTB;
	
	double WH1ZTT = WSVFF(MH1, ZMASS, TMASS, TMASS, three, 1., MH3, WH3, gzh1h3*H3FF*MSTMASS/v, 0., 0., 0.);
    double WH1ZBB = WSVFF(MH1, ZMASS, BMASS, BMASS, three, 1., MH3, WH3, gzh1h3*H3FF*MSBMASS/v, 0., 0., 0.);
	double WH1WTB = WSVFF(MH1, WMASS, TMASS, TMASS, three, 1., MHC, WHC, gwhch1*HCFF*(MSBMASS+MSTMASS)/v, 0., 0., 0.);
	double WH1WBT = WH1WTB;
	WH1 += WH1ZTT+WH1ZBB+WH1WTB+WH1WBT;
	
	double WH2ZTT = WSVFF(MH2, ZMASS, TMASS, TMASS, three, 1., MH3, WH3, gzh2h3*H3FF*MSTMASS/v, 0., 0., 0.);
    double WH2ZBB = WSVFF(MH2, ZMASS, BMASS, BMASS, three, 1., MH3, WH3, gzh2h3*H3FF*MSBMASS/v, 0., 0., 0.);
	double WH2WTB = WSVFF(MH2, WMASS, TMASS, TMASS, three, 1., MHC, WHC, gwhch2*HCFF*(MSBMASS+MSTMASS)/v, 0., 0., 0.);
	double WH2WBT = WH2WTB;
	WH2 += WH2ZTT+WH2ZBB+WH2WTB+WH2WBT;
	
	// Functions still in beta
	double WH3H2TT = WSSFF(MH3, MH2, TMASS, TMASS, three, 1., MH3, WH3, gh233*H3FF*MSTMASS/v, 0., 0., 0.);
	double WH3H1TT = WSSFF(MH3, MH1, TMASS, TMASS, three, 1., MH3, WH3, gh133*H3FF*MSTMASS/v, 0., 0., 0.);
	WH3 += WH3H1TT+WH3H2TT;
	
	double WHPH2TB = WSSFF(MH3, MH2, TMASS, BMASS, three, 1., MHC, WHC, gh2hmhp*HCFF*(MSTMASS+MSBMASS)/v, 0., 0., 0.);
	double WHPH1TB = WSSFF(MH3, MH1, TMASS, BMASS, three, 1., MHC, WHC, gh1hmhp*HCFF*(MSTMASS+MSBMASS)/v, 0., 0., 0.);
	WHC += WHPH1TB+WHPH2TB;
	
	double WH2H3TT = WSSFF(MH2, MH3, TMASS, TMASS, three, 1., MH3, WH3, gh233*H3FF*MSTMASS/v, 0., 0., 0.);
	double WH2HCTB = WSSFF(MH2, MHC, TMASS, BMASS, three, 1., MHC, WHC, gh2hmhp*HCFF*(MSTMASS+MSBMASS)/v, 0., 0., 0.);
	double WH2HCBT = WH2HCTB;
	double WH1H3TT = WSSFF(MH1, MH3, TMASS, TMASS, three, 1., MH3, WH3, gh133*H3FF*MSTMASS/v, 0., 0., 0.);
	double WH1HCTB = WSSFF(MH1, MHC, TMASS, BMASS, three, 1., MHC, WHC, gh1hmhp*HCFF*MSTMASS/v, 0., 0., 0.);
	double WH1HCBT = WH1HCTB;
	WH1 += WH1H3TT+WH1HCTB+WH1HCBT;
	WH2 += WH2H3TT+WH2HCTB+WH2HCBT;
	
	double WH2H1TT = WSSFF(MH2, MH1, TMASS, TMASS, three, 0., MH1, WH1, gh112*H1FF*MSTMASS/v, MH2, WH2, gh122*H2FF*MSTMASS/v);
	double WH2H1WW = WSSVV(one, MH2, MH1, WMASS, MH1, WH1, gh112*gwwh1, MH2, WH2, gh122*gwwh2);
	WH2H1WW += WSSVV2(one, MH2, MH1, WMASS, MHC, WHC, gwhch2*gwhch1);
	double WH2H1ZZ = WSSVV(two, MH2, MH1, ZMASS, MH1, WH1, gh112*gzzh1, MH2, WH2, gh122*gzzh2);
	WH2H1ZZ += WSSVV2(two, MH2, MH1, ZMASS, MH3, WH3, gzh2h3*gzh1h3);
	WH2 += WH2H1TT+WH2H1WW+WH2H1ZZ;
	//
	
	//For BRIDGE comparison
	fprintf(log,"\n");
	fprintf(log,"W(H1->X) =% 16.8e\n",WH1);
	if(WH1WW!=0.)   fprintf(log,"BR(H1->WW) =% 16.8e\n",WH1WW/WH1);
	if(WH1ZZ!=0.)   fprintf(log,"BR(H1->ZZ) =% 16.8e\n",WH1ZZ/WH1);
	if(WH1WHM!=0.)  fprintf(log,"BR(H1->W+H-) =% 16.8e\n",WH1WHM/WH1);
	if(WH1WHP!=0.)  fprintf(log,"BR(H1->W-H+) =% 16.8e\n",WH1WHP/WH1);
	if(WH1ZH3!=0.)  fprintf(log,"BR(H1->ZH3) =% 16.8e\n",WH1ZH3/WH1);
	if(WH1HMHP!=0.) fprintf(log,"BR(H1->H+H-) =% 16.8e\n",WH1HMHP/WH1);
	if(WH1H3H3!=0.) fprintf(log,"BR(H1->H3H3) =% 16.8e\n",WH1H3H3/WH1);
	if(WH1TT!=0.)   fprintf(log,"BR(H1->tt~) =% 16.8e\n",WH1TT/WH1);
	if(WH1BB!=0.)   fprintf(log,"BR(H1->bb~) =% 16.8e\n",WH1BB/WH1);
	if(WH1CC!=0.)   fprintf(log,"BR(H1->cc~) =% 16.8e\n",WH1CC/WH1);
	if(WH1SS!=0.)   fprintf(log,"BR(H1->ss~) =% 16.8e\n",WH1SS/WH1);
	if(WH1TATA!=0.) fprintf(log,"BR(H1->ta+ta-) =% 16.8e\n",WH1TATA/WH1);
	if(WH1MUMU!=0.) fprintf(log,"BR(H1->mu+mu-) =% 16.8e\n",WH1MUMU/WH1);
	if(WH1ZTT!=0.)  fprintf(log,"BR(H1->Ztt~) =% 16.8e\n",WH1ZTT/WH1);
	if(WH1ZBB!=0.)  fprintf(log,"BR(H1->Zbb~) =% 16.8e\n",WH1ZBB/WH1);
	if(WH1WTB!=0.)  fprintf(log,"BR(H1->W+bt~) =% 16.8e\n",WH1WTB/WH1);
	if(WH1WBT!=0.)  fprintf(log,"BR(H1->W-tb~) =% 16.8e\n",WH1WBT/WH1);
	//
	if(WH1H3TT!=0.)  fprintf(log,"BR(H1->H3tt~) =% 16.8e\n",WH1H3TT/WH1);
	if(WH1HCTB!=0.)  fprintf(log,"BR(H1->HCtb~) =% 16.8e\n",WH1HCTB/WH1);
	if(WH1HCBT!=0.)  fprintf(log,"BR(H1->HCbt~) =% 16.8e\n",WH1HCBT/WH1);
	//
	
	fprintf(log,"\n");
	fprintf(log,"W(H2->X) =% 16.8e\n",WH2);
	if(WH2WW!=0.)   fprintf(log,"BR(H2->WW) =% 16.8e\n",WH2WW/WH2);
	if(WH2ZZ!=0.)   fprintf(log,"BR(H2->ZZ) =% 16.8e\n",WH2ZZ/WH2);
	if(WH2WHM!=0.)  fprintf(log,"BR(H2->W+H-) =% 16.8e\n",WH2WHM/WH2);
	if(WH2WHP!=0.)  fprintf(log,"BR(H2->W-H+) =% 16.8e\n",WH2WHP/WH2);
	if(WH2ZH3!=0.)  fprintf(log,"BR(H2->ZH3) =% 16.8e\n",WH2ZH3/WH2);
	if(WH2ZH1!=0.)  fprintf(log,"BR(H2->ZH1) =% 16.8e\n",WH2ZH1/WH2);
	if(WH2HMHP!=0.) fprintf(log,"BR(H2->H+H-) =% 16.8e\n",WH2HMHP/WH2);
	if(WH2H3H3!=0.) fprintf(log,"BR(H2->H3H3) =% 16.8e\n",WH2H3H3/WH2);
	if(WH2H1H1!=0.) fprintf(log,"BR(H2->H1H1) =% 16.8e\n",WH2H1H1/WH2);
	if(WH2TT!=0.)   fprintf(log,"BR(H2->tt~) =% 16.8e\n",WH2TT/WH2);
	if(WH2BB!=0.)   fprintf(log,"BR(H2->bb~) =% 16.8e\n",WH2BB/WH2);
	if(WH2CC!=0.)   fprintf(log,"BR(H2->cc~) =% 16.8e\n",WH2CC/WH2);
	if(WH2SS!=0.)   fprintf(log,"BR(H2->ss~) =% 16.8e\n",WH2SS/WH2);
	if(WH2TATA!=0.) fprintf(log,"BR(H2->ta+ta-) =% 16.8e\n",WH2TATA/WH2);
	if(WH2MUMU!=0.) fprintf(log,"BR(H2->mu+mu-) =% 16.8e\n",WH2MUMU/WH2);
	if(WH2ZTT!=0.)  fprintf(log,"BR(H1->Ztt~) =% 16.8e\n",WH2ZTT/WH2);
	if(WH2ZBB!=0.)  fprintf(log,"BR(H1->Zbb~) =% 16.8e\n",WH2ZBB/WH2);
	if(WH2WTB!=0.)  fprintf(log,"BR(H1->W+bt~) =% 16.8e\n",WH2WTB/WH2);
	if(WH2WBT!=0.)  fprintf(log,"BR(H1->W-tb~) =% 16.8e\n",WH2WBT/WH2);
	//
	 if(WH2H3TT!=0.)  fprintf(log,"BR(H2->H3tt~) =% 16.8e\n",WH2H3TT/WH2);
	if(WH2HCTB!=0.)  fprintf(log,"BR(H2->HCtb~) =% 16.8e\n",WH2HCTB/WH2);
	if(WH2HCBT!=0.)  fprintf(log,"BR(H2->HCbt~) =% 16.8e\n",WH2HCBT/WH2);
	if(WH2H1TT!=0.)  fprintf(log,"BR(H2->H1tt~) =% 16.8e\n",WH2H1TT/WH2);
	if(WH2H1WW!=0.)  fprintf(log,"BR(H2->H1WW) =% 16.8e\n",WH2H1WW/WH2);
	if(WH2H1ZZ!=0.)  fprintf(log,"BR(H2->H1ZZ) =% 16.8e\n",WH2H1ZZ/WH2);
	 //
	
	fprintf(log,"\n");
	fprintf(log,"W(H3->X) =% 16.8e\n",WH3);
	if(WH3WW!=0.)   fprintf(log,"BR(H3->WW) =% 16.8e\n",WH3WW/WH3);
	if(WH3ZZ!=0.)   fprintf(log,"BR(H3->ZZ) =% 16.8e\n",WH3ZZ/WH3);
	if(WH3WHM!=0.)  fprintf(log,"BR(H3->W+H-) =% 16.8e\n",WH3WHM/WH3);
	if(WH3WHP!=0.)  fprintf(log,"BR(H3->W-H+) =% 16.8e\n",WH3WHP/WH3);
	if(WH3HMHP!=0.) fprintf(log,"BR(H3->H+H-) =% 16.8e\n",WH3HMHP/WH3);
	if(WH3ZH1!=0.)  fprintf(log,"BR(H3->ZH1) =% 16.8e\n",WH3ZH1/WH3);
	if(WH3ZH2!=0.)  fprintf(log,"BR(H3->ZH2) =% 16.8e\n",WH3ZH2/WH3);
	if(WH3TT!=0.)   fprintf(log,"BR(H3->tt~) =% 16.8e\n",WH3TT/WH3);
	if(WH3BB!=0.)   fprintf(log,"BR(H3->bb~) =% 16.8e\n",WH3BB/WH3);
	if(WH3CC!=0.)   fprintf(log,"BR(H3->cc~) =% 16.8e\n",WH3CC/WH3);
	if(WH3SS!=0.)   fprintf(log,"BR(H3->ss~) =% 16.8e\n",WH3SS/WH3);
	if(WH3TATA!=0.) fprintf(log,"BR(H3->ta+ta-) =% 16.8e\n",WH3TATA/WH3);
	if(WH3MUMU!=0.) fprintf(log,"BR(H3->mu+mu-) =% 16.8e\n",WH3MUMU/WH3);
	if(WH3WWZ!=0.)  fprintf(log,"BR(H3->WWZ) =% 16.8e\n",WH3WWZ/WH3);
	if(WH3ZZZ!=0.)  fprintf(log,"BR(H3->ZZZ) =% 16.8e\n",WH3ZZZ/WH3);
	if(WH3ZTT!=0.)  fprintf(log,"BR(H3->Ztt~) =% 16.8e\n",WH3ZTT/WH3);
	if(WH3ZBB!=0.)  fprintf(log,"BR(H3->Zbb~) =% 16.8e\n",WH3ZBB/WH3);
	if(WH3WTB!=0.)  fprintf(log,"BR(H3->W+bt~) =% 16.8e\n",WH3WTB/WH3);
	if(WH3WBT!=0.)  fprintf(log,"BR(H3->W-tb~) =% 16.8e\n",WH3WBT/WH3);
	//
	if(WH3H1TT!=0.) fprintf(log,"BR(H3->H1tt~) =% 16.8e\n",WH3H1TT/WH3);
	if(WH3H2TT!=0.) fprintf(log,"BR(H3->H2tt~) =% 16.8e\n",WH3H2TT/WH3);
	//
	
	fprintf(log,"\n");
	fprintf(log,"W(HC->X) =% 16.8e\n",WHC);
	if(WHPWH1!=0.)  fprintf(log,"BR(HC->W+H1) =% 16.8e\n",WHPWH1/WHC);
	if(WHPWH2!=0.)  fprintf(log,"BR(HC->W+H2) =% 16.8e\n",WHPWH2/WHC);
	if(WHPWH3!=0.)  fprintf(log,"BR(HC->W+H3) =% 16.8e\n",WHPWH3/WHC);
	if(WHPTB!=0.)   fprintf(log,"BR(HC->bt~) =% 16.8e\n",WHPTB/WHC);
	if(WHPCS!=0.)   fprintf(log,"BR(HC->sc~) =% 16.8e\n",WHPCS/WHC);
	if(WHPVTTA!=0.) fprintf(log,"BR(HC->tavt) =% 16.8e\n",WHPVTTA/WHC);
	if(WHPVMMU!=0.) fprintf(log,"BR(HC->muvm) =% 16.8e\n",WHPVMMU/WHC);
	if(WHPWZZ!=0.)  fprintf(log,"BR(HC->WZZ) =% 16.8e\n",WHPWZZ/WHC);
	if(WHPWWW!=0.)  fprintf(log,"BR(HC->WWW) =% 16.8e\n",WHPWWW/WHC);
	if(WHPWTT!=0.)  fprintf(log,"BR(HC->Wtt~) =% 16.8e\n",WHPWTT/WHC);
	if(WHPWBB!=0.)  fprintf(log,"BR(HC->Wbb~) =% 16.8e\n",WHPWBB/WHC);
	if(WHPZTB!=0.)  fprintf(log,"BR(HC->Zbt~) =% 16.8e\n",WHPZTB/WHC);
	//
	if(WHPH1TB!=0.) fprintf(log,"BR(HC->H1bt~) =% 16.8e\n",WHPH1TB/WHC);
	if(WHPH2TB!=0.) fprintf(log,"BR(HC->H2bt~) =% 16.8e\n",WHPH2TB/WHC);
	//
	
	//write to file
	fprintf(output,"DECAY  25    % 16.13E      # H1 width\n",WH1);
	fprintf(output,"DECAY  35    % 16.13E      # H2 width\n",WH2);
	fprintf(output,"DECAY  36    % 16.13E      # H3 width\n",WH3);
	fprintf(output,"DECAY  37    % 16.13E      # H+ width\n",WHC);
	
	
	//MGUSER stuff
	fprintf(output,"BLOCK MGUSER\n");
	fprintf(output,"         1     %16.13e   # TanBeta ,v2/v1\n",tb);
	fprintf(output,"         2     %16.13e   # SinAlpha ,Sine of the h1/h2 mixing angle\n",sa);
	fprintf(output,"         3     %16.13e   # SinGamma ,Sine of the h1/h2 top angle\n",sgamma);
	// Close all files
	// ***************
	
	
	fprintf(stdout,"The generation of %s is complete\n",argv[2]);
	if(argc>=4) fprintf(stdout,"Rough branching ratio information stored in %s \n",argv[3]);
	else fprintf(stdout,"Rough branching ratio information stored in 2hdm4tccalc.log\n");
	fprintf(stdout,"The widths are generally reliable.\n");
	fprintf(stdout,"However, use of BRIDGE is recommended\n");
	fprintf(stdout,"for more precise branching ratios and widths\n");
	
	fclose(data);
	fclose(output);
	fclose(log);
	
	return 0;
	
}
