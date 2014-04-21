#include <R.h>

/* A typical C trick to get readable names for the parameters is #define
   This method is simple and efficient, but there are, of course,
   other possibilities that use dynamic variables.
*/
static double parms[11];

#define C_hu_nh4 parms[0]
#define C_micb_hum parms[1]
#define C_pr_nh4 parms[2]
#define C_pr_micb parms[3]
#define C_micb_nh4 parms[4]
#define K_nh4mic parms[5]
#define K_no3mic parms[6]
#define K_nh4no3 parms[7]
#define K_no3ng parms[8]
#define K_nh4nh3 parms[9]
#define K_no3nh4 parms[10]


/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
    int N= 11 ;
    odeparms(&N, parms);
	}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
	double m,h,s,j,r,i_a,i_n,n,d,v,nd;
	double v15,n15,i_a_15,i_n_15,d15,nd15;
	double r15,m15,h15;
	double dNH4, dNO3, dMib, dPR, dHum, dTON;
	double dNH4_N15, dNO3_N15, dMib_N15, dHum_N15, dTON_N15, dPR_N15;
	double NH4, NO3, NH4_N15, NO3_N15;
	double Mib, PR, Hum, TON,TON_N15,PR_N15;
	double Mib_N15,Hum_N15;
	double dNH3, dN2O, dNH3_N15, dN2O_N15;
	double NH3, N2O, NH3_N15, N2O_N15;
//		
	
	NH4 = y[0];
	NO3 = y[1];
	
	NH4_N15 = y[2];
	NO3_N15 = y[3];
	
	TON_N15 = y[4];
	TON = y[5];
	
	Mib = y[6];
	PR = y[7];
	Hum = y[8];
	
	Mib_N15 = y[9];
	Hum_N15 = y[10];
	
	NH3 = y[11];
	N2O = y[12];
	NH3_N15 = y[13];
	N2O_N15 = y[14];
	
	PR_N15 = y[15];
//	
	m  = C_hu_nh4;
	h  = C_micb_hum;
	s  = C_pr_nh4;
	j  = C_pr_micb;
	r  = C_micb_nh4;

	i_a = K_nh4mic * NH4;
	i_n = K_no3mic * NO3;
	n  = K_nh4no3 * NH4;
					
	d  = K_no3ng * NO3;
	v  = K_nh4nh3 * NH4;
	nd = K_no3nh4 * NO3;
					
//#assume conversion from micb to humus will take very long time, then the h,s,j and m for 15N balance was omitted
					
	v15 = K_nh4nh3 * NH4_N15;
	d15 = K_no3ng * NO3_N15;
	n15 = K_nh4no3 * NH4_N15;
	nd15 = K_no3nh4 * NO3_N15;	
	i_a_15 = K_nh4mic * NH4_N15;
	i_n_15 = K_no3mic * NO3_N15;
	r15 = r;	
	h15 = h;
	m15 = m;

    if (ip[0] <1) error("nout should be at least 1");
	
	dNH4 = m + s + nd - n - i_a - v + r ;
	dNO3 = n - d - nd - i_n;
	dMib = i_a + i_n +j - h - r;
	dPR = - s - j;
	dHum = h - m;
	dTON = dMib + dPR + dHum;
	dNH3 = v;
	dN2O = d;
					
	dNH4_N15 = (nd15 -n15 - i_a_15 - v15 + r15 + m15 + s);
	dNO3_N15 = (n15 - i_n_15 - d15 - nd15);
	dMib_N15 = (i_a_15 + i_n_15 - r15 -h15 + j);
	dHum_N15 = h15 - m15;
	dPR_N15  = 0;
	dTON_N15 = (dMib_N15+dHum_N15+dPR_N15);
	dNH3_N15 = v15;
	dN2O_N15 = d15;

    ydot[0] = dNH4;
    ydot[1] = dNO3;
	ydot[2] = dNH4_N15;
	ydot[3] = dNO3_N15;
	ydot[4] = dTON_N15;
	
    ydot[5] = dTON;
	ydot[6] = dMib;
	ydot[7] = dPR;
	ydot[8] = dHum;
	
	ydot[9] = dMib_N15;
	ydot[10]= dHum_N15;
	
	ydot[11]= dNH3;
	ydot[12]= dN2O;
	ydot[13]= dNH3_N15;
	ydot[14]= dN2O_N15;
	ydot[15]= dPR_N15;

	yout[0] = NH4;
	yout[1] = NO3;
//	yout[2] = TON;
}

/* It is also possible to orivide a Jacobian (otionally) */
//void jlotka(int *neq, double *t, double *y, int *ml,
//           int *mu, double *pd, int *nrowpd, double *yout, int*ip) {
//   // ~~~~~ optional: code for the Jacobian 
//}


