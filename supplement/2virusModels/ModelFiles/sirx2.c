/* file sirx2.c */
#include <R.h>
#include <math.h>
static double parms[23];
#define eps parms[0]
#define phase1 parms[1]
#define phase2 parms[2]
#define beta1 parms[3]
#define beta1p parms[4]
#define beta1tilde parms[5]
#define beta2 parms[6]
#define beta2p parms[7]
#define beta2tilde parms[8]
#define gam1 parms[9]
#define gam2 parms[10]
#define gam1p parms[11]
#define gam2p parms[12]
#define gam1tilde parms[13]
#define gam2tilde parms[14]
#define rho1 parms[15]
#define rho2 parms[16]
#define rho1p parms[17]
#define rho2p parms[18]
#define rho1tilde parms[19]
#define rho2tilde parms[20]
#define T_pd parms[21]
#define Npop parms[22]
/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=23;
odeparms(&N, parms);
}
/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
double *yout, int *ip)
{
if (ip[0] <1) error("nout should be at least 1");
ydot[0] = -y[0]*(beta1)*(1+eps*sin(2*M_PI*(*t-phase1)/T_pd))/Npop*(y[1]+y[3]+y[7]) - y[0]*(beta2)*(1+eps*sin(2*M_PI*(*t-phase2)/T_pd))/Npop*(y[2]+y[3]+y[6]) + rho1*y[4] + rho2*y[5];
ydot[1] = y[0]*(beta1)*(1+eps*sin(2*M_PI*(*t-phase1)/T_pd))/Npop*(y[1]+y[3]+y[7]) - y[1]*(beta2p)*(1+eps*sin(2*M_PI*(*t-phase2)/T_pd))/Npop*(y[2]+y[3]+y[6])-gam1*y[1]+(rho2p)*y[7];
ydot[2] = y[0]*(beta2)*(1+eps*sin(2*M_PI*(*t-phase2)/T_pd))/Npop*(y[2]+y[3]+y[6]) - y[2]*(beta1p)*(1+eps*sin(2*M_PI*(*t-phase1)/T_pd))/Npop*(y[1]+y[3]+y[7])-gam2*y[2]+(rho1p)*y[6];
ydot[3] = y[1]*(beta2p)*(1+eps*sin(2*M_PI*(*t-phase2)/T_pd))/Npop*(y[2]+y[3]+y[6])+y[2]*(beta1p)*(1+eps*sin(2*M_PI*(*t-phase1)/T_pd))/Npop*(y[1]+y[3]+y[7])-(gam1p+gam2p)*y[3];
ydot[4] = gam1*y[1]-rho1*y[4]-(beta2tilde)*(1+eps*sin(2*M_PI*(*t-phase2)/T_pd))/Npop*y[4]*(y[2]+y[3]+y[6])+(rho2tilde)*y[8];
ydot[5] = gam2*y[2]-rho2*y[5]-(beta1tilde)*(1+eps*sin(2*M_PI*(*t-phase1)/T_pd))/Npop*y[5]*(y[1]+y[3]+y[7])+(rho1tilde)*y[8];
ydot[6] = (gam1p)*y[3]-(gam2tilde)*y[6]-(rho1p)*y[6]+(beta2tilde)*(1+eps*sin(2*M_PI*(*t-phase2)/T_pd))/Npop*y[4]*(y[2]+y[3]+y[6]);
ydot[7] = (gam2p)*y[3]-(gam1tilde)*y[7]-(rho2p)*y[7]+(beta1tilde)*(1+eps*sin(2*M_PI*(*t-phase1)/T_pd))/Npop*y[5]*(y[1]+y[3]+y[7]);
ydot[8] = (gam1tilde)*y[7]+(gam2tilde)*y[6]-(rho1tilde+rho2tilde)*y[8];
yout[0] = sin((*t-phase1)*2*M_PI/T_pd);
}
/* The Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu,
double *pd, int *nrowpd, double *yout, int *ip)
{
pd[0] = -beta1*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*(y[1] + y[3] + y[7])/Npop - beta2*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*(y[2] + y[3] + y[6])/Npop;
pd[(*nrowpd)] = -y[0]*beta1*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop;
pd[2*(*nrowpd)] = -y[0]*beta2*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop;
pd[3*(*nrowpd)] = -y[0]*beta1*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop - y[0]*beta2*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop;
pd[4*(*nrowpd)] = rho1;
pd[5*(*nrowpd)] = rho2;
pd[6*(*nrowpd)] = -y[0]*beta2*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop;
pd[7*(*nrowpd)] = -y[0]*beta1*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop;
pd[8*(*nrowpd)] = 0;
pd[1] = beta1*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*(y[1] + y[3] + y[7])/Npop;
pd[(*nrowpd)+1] = y[0]*beta1*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop - (beta2p)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*(y[2] + y[3] + y[6])/Npop - gam1;
pd[2*(*nrowpd)+1] = -(beta2p)*y[1]*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop;
pd[3*(*nrowpd)+1] = y[0]*beta1*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop - (beta2p)*y[1]*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop;
pd[4*(*nrowpd)+1] = 0;
pd[5*(*nrowpd)+1] = 0;
pd[6*(*nrowpd)+1] = -(beta2p)*y[1]*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop;
pd[7*(*nrowpd)+1] = y[0]*beta1*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop + rho2p;
pd[8*(*nrowpd)+1] = 0;
pd[2] = beta2*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*(y[2] + y[3] + y[6])/Npop;
pd[(*nrowpd)+2] = -y[2]*(beta1 + beta1p)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop;
pd[2*(*nrowpd)+2] = y[0]*beta2*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop - (beta1p)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*(y[1] + y[3] + y[7])/Npop - gam2;
pd[3*(*nrowpd)+2] = y[0]*beta2*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop - y[2]*(beta1p)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop;
pd[4*(*nrowpd)+2] = 0;
pd[5*(*nrowpd)+2] = 0;
pd[6*(*nrowpd)+2] = y[0]*beta2*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop + rho1p; 
pd[7*(*nrowpd)+2] = -y[2]*(beta1p)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop;
pd[8*(*nrowpd)+2] = 0;
pd[3] = 0;
pd[(*nrowpd)+3] = (beta2p)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*(y[2] + y[3] + y[6])/Npop + y[2]*(beta1p)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop;
pd[2*(*nrowpd)+3] = (beta2p)*y[1]*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop + (beta1p)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*(y[1] + y[3] + y[7])/Npop; 
pd[3*(*nrowpd)+3] =  (beta2p)*y[1]*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop + y[2]*(beta1p)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop - gam1p - gam2p; 
pd[4*(*nrowpd)+3] = 0;
pd[5*(*nrowpd)+3] = 0;
pd[6*(*nrowpd)+3] = (beta2p)*y[1]*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))/Npop;
pd[7*(*nrowpd)+3] = y[2]*(beta1p)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))/Npop;
pd[8*(*nrowpd)+3] = 0;
pd[4] = 0;
pd[(*nrowpd)+4] = gam1;
pd[2*(*nrowpd)+4] = -(beta2tilde)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*y[4]/Npop;
pd[3*(*nrowpd)+4] = -(beta2tilde)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*y[4]/Npop; 
pd[4*(*nrowpd)+4] = -rho1 - (beta2tilde)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*(y[2] + y[3] + y[6])/Npop;
pd[5*(*nrowpd)+4] = 0;
pd[6*(*nrowpd)+4] = -(beta2tilde)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*y[4]/Npop;
pd[7*(*nrowpd)+4] = 0;
pd[8*(*nrowpd)+4] = rho2tilde;
pd[5] = 0;   
pd[(*nrowpd)+5] = -(beta1tilde)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*y[5]/Npop;
pd[2*(*nrowpd)+5] = gam2;
pd[3*(*nrowpd)+5] = -(beta1tilde)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*y[5]/Npop;
pd[4*(*nrowpd)+5] = 0;
pd[5*(*nrowpd)+5] = -rho2 - (beta1tilde)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*(y[1] + y[3] + y[7])/Npop;
pd[6*(*nrowpd)+5] =  0;
pd[7*(*nrowpd)+5] = -(beta1tilde)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*y[5]/Npop;
pd[8*(*nrowpd)+5] = rho1tilde;
pd[6] = 0;
pd[(*nrowpd)+6] = 0;
pd[2*(*nrowpd)+6] = (beta2tilde)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*y[4]/Npop; 
pd[3*(*nrowpd)+6] = gam1p + (beta2tilde)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*y[4]/Npop; 
pd[4*(*nrowpd)+6] = (beta2tilde)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*(y[2] + y[3] + y[6])/Npop; 
pd[5*(*nrowpd)+6] = 0;
pd[6*(*nrowpd)+6] = - gam2tilde - rho1p + (beta2tilde)*(1 + eps*sin(2*M_PI*(*t - phase2)/T_pd))*y[4]/Npop;
pd[7*(*nrowpd)+6] = 0;
pd[8*(*nrowpd)+6] = 0;
pd[7] = 0;
pd[(*nrowpd)+7] = (beta1tilde)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*y[5]/Npop;
pd[2*(*nrowpd)+7] = 0;
pd[3*(*nrowpd)+7] = gam2p + (beta1tilde)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*y[5]/Npop;
pd[4*(*nrowpd)+7] = 0;
pd[5*(*nrowpd)+7] = (beta1tilde)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*(y[1] + y[3] + y[7])/Npop;
pd[6*(*nrowpd)+7] = 0;
pd[7*(*nrowpd)+7] =  - gam1tilde - rho2p + (beta1tilde)*(1 + eps*sin(2*M_PI*(*t - phase1)/T_pd))*y[5]/Npop; 
pd[8*(*nrowpd)+7] = 0;
pd[8] = 0;
pd[(*nrowpd)+8] = 0;
pd[2*(*nrowpd)+8] = 0;
pd[3*(*nrowpd)+8] = 0;
pd[4*(*nrowpd)+8] = 0;
pd[5*(*nrowpd)+8] = 0;
pd[6*(*nrowpd)+8] = gam2tilde;
pd[7*(*nrowpd)+8] = gam1tilde; 
pd[8*(*nrowpd)+8] = - rho1tilde - rho2tilde;
}
/* END file sirx2.c */
