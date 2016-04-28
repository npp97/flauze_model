#include <R.h>

/* A typical C trick to get readable names for the parameters is #define
   This method is simple and efficient, but there are, of course,
   other possibilities that use dynamic variables.
*/
static double parms[29];

#define C_NRec_2_NH4 parms[0] 
#define C_NH4_2_NRec parms[1]
#define C_NH4_2_Nlab parms[2]
#define C_NO3_2_NRec parms[3]
#define C_NRec_2_NO3 parms[4]

#define K_NLab_2_NH4 parms[5]
#define K_NH4_2_NO2nit parms[6]   
#define K_NH4_2_NH4ads parms[7]
#define K_NH4ads_2_NH4 parms[8]
#define K_NO3_2_NH4 parms[9]
#define K_NO2nit_2_NO3 parms[10]
#define K_NO2nit_2_N2Onit parms[11]

#define K_NO3_2_NO2den parms[12]
#define K_NO2den_2_N2Oden parms[13]
#define K_NO2den_2_N2Ocod parms[14]
#define K_NRec_2_NO2org parms[15]
#define K_NO2org_2_N2Oorg parms[16]
#define K_NO2org_2_N2Ocod parms[17]

#define K_N2Oorg_2_N2O parms[18]
#define K_N2Ocod_2_N2O parms[19]
#define K_N2Oden_2_N2O parms[20]
#define K_N2Onit_2_N2O parms[21]

#define K_N2Oorg_2_N2 parms[22]
#define K_N2Ocod_2_N2 parms[23]
#define K_N2Oden_2_N2 parms[24]
#define K_N2Onit_2_N2 parms[25]

#define f_N2Oorg parms[26]
#define f_N2Ocod parms[27]
#define f_N2Oden parms[28]

/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
  int N= 29 ;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
  //----------Full
  double NLab_2_NH4,NH4_2_NLab,NRec_2_NH4,NH4_2_NRec,NH4_2_NO2nit,NH4_2_NH4ads,NH4ads_2_NH4,NO3_2_NH4,NO2nit_2_NO3,NO2nit_2_N2Onit;
  double NO3_2_NRec,NRec_2_NO3,NO3_2_NO2den,NO2den_2_N2Oden,NO2den_2_N2Ocod,NRec_2_NO2org,NO2org_2_N2Oorg,NO2org_2_N2Ocod;
  double N2Oorg_2_N2O,N2Ocod_2_N2O,N2Oden_2_N2O,N2Onit_2_N2O,N2Oorg_2_N2,N2Ocod_2_N2,N2Oden_2_N2,N2Onit_2_N2;
  
  double NLab,NH4,NO3,NH4ads,NRec,NO2nit,NO2den,NO2org,N2O,N2Oorg,N2Ocod,N2Oden,N2Onit;
  double dNLab,dNH4,dNO3,dNH4ads,dNRec,dNO2nit,dNO2den,dNO2org,dN2O,dN2Oorg,dN2Ocod,dN2Oden,dN2Onit,dN2;  
  double f_N2Onit;
  
  //-----------N15
  double NLab_2_NH4_15,NH4_2_NLab_15,NRec_2_NH4_15,NH4_2_NRec_15,NH4_2_NO2nit_15,NH4_2_NH4ads_15;
  double NH4ads_2_NH4_15,NO3_2_NH4_15,NO2nit_2_NO3_15,NO2nit_2_N2Onit_15;
  double NO3_2_NRec_15,NRec_2_NO3_15,NO3_2_NO2den_15,NO2den_2_N2Oden_15,NO2den_2_N2Ocod_15,NRec_2_NO2org_15,NO2org_2_N2Oorg_15,NO2org_2_N2Ocod_15;
  double N2Oorg_2_N2O_15,N2Ocod_2_N2O_15,N2Oden_2_N2O_15,N2Onit_2_N2O_15,N2Oorg_2_N2_15,N2Ocod_2_N2_15,N2Oden_2_N2_15,N2Onit_2_N2_15;
  
  double NLab_15,NH4_15,NO3_15,NH4ads_15,NRec_15,NO2nit_15,NO2den_15,NO2org_15,N2O_15,N2Oorg_15,N2Ocod_15,N2Oden_15,N2Onit_15;
  double dNLab_15,dNH4_15,dNO3_15,dNH4ads_15,dNRec_15,dNO2nit_15,dNO2den_15,dNO2org_15,dN2O_15,dN2Oorg_15,dN2Ocod_15,dN2Oden_15,dN2Onit_15,dN2_15;  
  //

  NLab = y[0];
  NLab_15 = y[1];
  
  NH4  = y[2];
  NH4_15 = y[3];
  
  NO3 = y[4];
  NO3_15 = y[5];
  
  NH4ads =y[6];
  NH4ads_15 =y[7];
  
  NRec = y[8];
  NRec_15 = y[9];
  
  NO2nit = y[10];
  NO2nit_15 = y[11];
  NO2den = y[12];
  NO2den_15 = y[13];
  
  N2O = y[14];
  N2O_15 = y[15];
  
  N2Oden = f_N2Oden * N2O;
  N2Ocod = f_N2Ocod * N2O;
  N2Oorg = f_N2Oorg * N2O;
  f_N2Onit=1-f_N2Oden-f_N2Ocod-f_N2Oorg;
  N2Onit = f_N2Onit * N2O;
  
  N2Oden_15 = f_N2Oden * N2O_15;
  N2Ocod_15 = f_N2Ocod * N2O_15;
  N2Oorg_15 = f_N2Oorg * N2O_15;
//  f_N2Onit=1-f_N2Oden-f_N2Ocod-f_N2Oorg;
  N2Onit_15 = f_N2Onit * N2O_15;
  
    
//----full  
  NLab_2_NH4 = K_NLab_2_NH4*NLab;
  NH4_2_NLab = C_NH4_2_Nlab;
  
  dNLab = NH4_2_NLab - NLab_2_NH4; //1
  
  NRec_2_NH4 = C_NRec_2_NH4;
  NH4_2_NRec = C_NH4_2_NRec;
  NH4_2_NO2nit = K_NH4_2_NO2nit * NH4;
  NH4_2_NH4ads = K_NH4_2_NH4ads * NH4;
  NH4ads_2_NH4 = K_NH4ads_2_NH4 * NH4ads;
  NO3_2_NH4 = K_NO3_2_NH4 * NO3;
  
  dNH4 = NRec_2_NH4+NH4ads_2_NH4+NO3_2_NH4-NH4_2_NH4ads-NH4_2_NLab-NH4_2_NO2nit-NH4_2_NRec; //2
  
  dNH4ads = NH4_2_NH4ads-NH4ads_2_NH4; //3
  
  dNRec = NH4_2_NRec+NO3_2_NRec-NRec_2_NH4-NRec_2_NO2org-NRec_2_NO3; //4
  
  NO2nit_2_NO3 = K_NO2nit_2_NO3 * NO2nit;
  NO2nit_2_N2Onit = K_NO2nit_2_N2Onit * NO2nit;
  dNO2nit = NH4_2_NO2nit-NO2nit_2_NO3-NO2nit_2_N2Onit;//5
  
  NO3_2_NRec = C_NO3_2_NRec;
  NRec_2_NO3 = C_NRec_2_NO3;
  NO3_2_NO2den = K_NO3_2_NO2den * NO3;
  dNO3=NO2nit_2_NO3+NRec_2_NO3-NO3_2_NH4-NO3_2_NO2den+NO3_2_NRec; //6
  
  NO2den_2_N2Oden = K_NO2den_2_N2Oden * NO2den;
  NO2den_2_N2Ocod = K_NO2den_2_N2Ocod * NO2den;
  dNO2den = NO3_2_NO2den - NO2den_2_N2Ocod - NO2den_2_N2Oden; //7
  
  NRec_2_NO2org = K_NRec_2_NO2org * NRec;
  NO2org_2_N2Oorg = K_NO2org_2_N2Oorg * NO2org_2_N2Oorg;
  NO2org_2_N2Ocod = K_NO2org_2_N2Ocod * NO2org_2_N2Ocod;
  dNO2org = NRec_2_NO2org - NO2org_2_N2Ocod - NO2org_2_N2Oorg; //8
  
  N2Oorg_2_N2O = K_N2Oorg_2_N2O * N2Oorg;
  N2Ocod_2_N2O = K_N2Ocod_2_N2O * N2Ocod;
  N2Oden_2_N2O = K_N2Oden_2_N2O * N2Oden;
  N2Onit_2_N2O = K_N2Onit_2_N2O * N2Onit;
  
  dN2O = N2Oorg_2_N2O+N2Ocod_2_N2O+N2Oden_2_N2O+N2Onit_2_N2O; //9
  
  N2Oorg_2_N2 = K_N2Oorg_2_N2 * N2Oorg;
  N2Ocod_2_N2 = K_N2Ocod_2_N2 * N2Ocod;
  N2Oden_2_N2 = K_N2Oden_2_N2 * N2Oden;
  N2Onit_2_N2 = K_N2Onit_2_N2 * N2Onit;
  
  dN2Oorg = NO2org_2_N2Oorg - N2Oorg_2_N2O - N2Oorg_2_N2; //10
  dN2Ocod = NO2org_2_N2Ocod + NO2den_2_N2Ocod - N2Oden_2_N2O- N2Oden_2_N2; //11
  dN2Oden = NO2den_2_N2Oden - N2Oden_2_N2O - N2Oden_2_N2; //12
  dN2Onit = NO2nit_2_N2Onit - N2Onit_2_N2O - N2Onit_2_N2; //13
  
  dN2O = dN2Oorg + dN2Ocod + dN2Oden + dN2Onit; //14
  dN2 =N2Oorg_2_N2 + N2Oden_2_N2 + N2Ocod_2_N2 + N2Onit_2_N2; //15
  
 //------N15 
  NLab_2_NH4_15 = K_NLab_2_NH4*NLab_15;
  NH4_2_NLab_15 = C_NH4_2_Nlab;
  
  dNLab_15 = NH4_2_NLab_15 - NLab_2_NH4_15; //1
  
  NRec_2_NH4_15 = C_NRec_2_NH4;
  NH4_2_NRec_15 = C_NH4_2_NRec;
  NH4_2_NO2nit_15 = K_NH4_2_NO2nit * NH4_15;
  NH4_2_NH4ads_15 = K_NH4_2_NH4ads * NH4_15;
  NH4ads_2_NH4_15 = K_NH4ads_2_NH4 * NH4ads_15;
  NO3_2_NH4_15 = K_NO3_2_NH4 * NO3_15;
  
  dNH4_15 = NRec_2_NH4_15+NH4ads_2_NH4_15+NO3_2_NH4_15-NH4_2_NH4ads_15-NH4_2_NLab_15-NH4_2_NO2nit_15-NH4_2_NRec_15; //2
  
  dNH4ads_15 = NH4_2_NH4ads_15-NH4ads_2_NH4_15; //3
  
  dNRec_15 = NH4_2_NRec_15+NO3_2_NRec_15-NRec_2_NH4_15-NRec_2_NO2org_15-NRec_2_NO3_15; //4
  
  NO2nit_2_NO3_15 = K_NO2nit_2_NO3 * NO2nit_15;
  NO2nit_2_N2Onit_15 = K_NO2nit_2_N2Onit * NO2nit_15;
  dNO2nit_15 = NH4_2_NO2nit_15-NO2nit_2_NO3_15-NO2nit_2_N2Onit_15;//5
  
  NO3_2_NRec_15 = C_NO3_2_NRec;
  NRec_2_NO3_15 = C_NRec_2_NO3;
  NO3_2_NO2den_15 = K_NO3_2_NO2den * NO3_15;
  dNO3_15=NO2nit_2_NO3_15+NRec_2_NO3_15-NO3_2_NH4_15-NO3_2_NO2den_15+NO3_2_NRec_15; //6
  
  NO2den_2_N2Oden_15 = K_NO2den_2_N2Oden * NO2den_15;
  NO2den_2_N2Ocod_15 = K_NO2den_2_N2Ocod * NO2den_15;
  dNO2den_15 = NO3_2_NO2den_15 - NO2den_2_N2Ocod_15 - NO2den_2_N2Oden_15; //7
  
  NRec_2_NO2org_15 = K_NRec_2_NO2org * NRec_15;
  NO2org_2_N2Oorg_15 = K_NO2org_2_N2Oorg * NO2org_2_N2Oorg_15;
  NO2org_2_N2Ocod_15 = K_NO2org_2_N2Ocod * NO2org_2_N2Ocod_15;
  dNO2org_15 = NRec_2_NO2org_15 - NO2org_2_N2Ocod_15 - NO2org_2_N2Oorg_15; //8
  
  N2Oorg_2_N2O_15 = K_N2Oorg_2_N2O * N2Oorg_15;
  N2Ocod_2_N2O_15 = K_N2Ocod_2_N2O * N2Ocod_15;
  N2Oden_2_N2O_15 = K_N2Oden_2_N2O * N2Oden_15;
  N2Onit_2_N2O_15 = K_N2Onit_2_N2O * N2Onit_15;
  
  dN2O_15 = N2Oorg_2_N2O_15+N2Ocod_2_N2O_15+N2Oden_2_N2O_15+N2Onit_2_N2O_15; //9
  
  N2Oorg_2_N2_15 = K_N2Oorg_2_N2 * N2Oorg_15;
  N2Ocod_2_N2_15 = K_N2Ocod_2_N2 * N2Ocod_15;
  N2Oden_2_N2_15 = K_N2Oden_2_N2 * N2Oden_15;
  N2Onit_2_N2_15 = K_N2Onit_2_N2 * N2Onit_15;
  
  dN2Oorg_15 = NO2org_2_N2Oorg_15 - N2Oorg_2_N2O_15 - N2Oorg_2_N2_15; //10
  dN2Ocod_15 = NO2org_2_N2Ocod_15 + NO2den_2_N2Ocod_15 - N2Oden_2_N2O_15- N2Oden_2_N2_15; //11
  dN2Oden_15 = NO2den_2_N2Oden_15 - N2Oden_2_N2O_15 - N2Oden_2_N2_15; //12
  dN2Onit_15 = NO2nit_2_N2Onit_15 - N2Onit_2_N2O_15 - N2Onit_2_N2_15; //13
  
  dN2O_15 = dN2Oorg_15 + dN2Ocod_15 + dN2Oden_15 + dN2Onit_15; //14
  dN2_15 =N2Oorg_2_N2_15 + N2Oden_2_N2_15 + N2Ocod_2_N2_15 + N2Onit_2_N2_15; //15
  
  //	 double dNLab,dNH4,dNO3,dNH4ads,dNRec,dNO2nit,dNO2den,dNO2org,dN2O,dN2Oorg,dN2Ocod,dN2Oden,dN2Onit,dN2; 	
  ydot[0] = dNLab;
  ydot[1] = dNLab_15;
  ydot[2] = dNH4;
  ydot[3] = dNH4_15;
  ydot[4] = dNO3;
  ydot[5] = dNO3_15;
  
  ydot[6] = dNH4ads;
  ydot[7] = dNH4ads_15;
  ydot[8] = dNRec;
  ydot[9] = dNRec_15;
  
  ydot[10] =dNO2nit;
  ydot[11]= dNO2nit_15;
  
  ydot[12]= dNO2den;
  ydot[13]= dNO2den_15;
  
  ydot[14]= dNO2org;
  ydot[15]= dNO2org_15;
  
  ydot[16]= dN2O;
  ydot[17]= dN2O_15;
  
  ydot[18]= dN2Oorg;
  ydot[19]= dN2Oorg_15;
  
  ydot[20]= dN2Oden;
  ydot[21]= dN2Oden_15;
  
  ydot[22]= dN2Ocod;
  ydot[23]= dN2Ocod_15;
  
  ydot[24]= dN2Onit;
  ydot[25]= dN2Onit_15;
  
  
  yout[0] = NH4;
  yout[1] = NO3;
  yout[2] = NH4_15;
  yout[3] = NO3_15;
  yout[4] = N2O;
  yout[5] = N2O_15;
}