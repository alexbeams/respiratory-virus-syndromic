/* file sir.c */
#include <R.h>
#include <math.h>

static double parms[7];

#define eps parms[0]
#define phase parms[1]
#define beta parms[2]
#define gam parms[3]
#define rho parms[4]
#define T_pd parms[5]
#define Npop parms[6]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
	int N=7;
	odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
	if (ip[0] <1) error("nout should be at least 1");
	ydot[0] = -y[0]*beta*(1+eps*cos(2*M_PI*(*t-phase)/T_pd))/Npop*exp(y[1]) + rho*y[2];
	ydot[1] =  y[0]*beta*(1+eps*cos(2*M_PI*(*t-phase)/T_pd))/Npop - gam;
	ydot[2] = gam*exp(y[1]) - rho*y[2];
	yout[0] = cos((*t-phase)*2*M_PI/T_pd);
}


/* The Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu,
double *pd, int *nrowpd, double *yout, int *ip)
{
pd[0] 			= -beta*(1 + eps*cos(2*M_PI*(*t - phase)/T_pd))*exp(y[1])/Npop;
pd[(*nrowpd)] 		= -beta*(1 + eps*cos(2*M_PI*(*t - phase)/T_pd))*y[0]*exp(y[1])/Npop;
pd[2*(*nrowpd)] 	= rho;
pd[1] 			= beta*(1 + eps*cos(2*M_PI*(*t - phase)/T_pd))/Npop;
pd[(*nrowpd)+1] 	= 0; 
pd[2*(*nrowpd)+1] 	= 0;
pd[2] 			= 0;
pd[(*nrowpd)+2] 	= gam*exp(y[1]);
pd[2*(*nrowpd)+2] 	= -rho;
}
/* END file sir.c */
