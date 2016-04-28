#include <R.h>

/* A typical C trick to get readable names for the parameters is #define
   This method is simple and efficient, but there are, of course,
   other possibilities that use dynamic variables.
*/
static double parms[9];

#define C_NRec_2_NH4 parms[0] //3 M_NREc
#define C_NH4_2_NLab parms[1] //2 I_NH4
#define C_NO3_2_NRec parms[2] //4 I_NO3
#define C_NRec_2_NO3 parms[3] //5 O_NRec
#define C_NO3sto_2_NO3 parms[4] //9 R_NO3s

#define K_NLab_2_NH4 parms[5] //1 M_Nlab
#define K_NH4_2_NO3 parms[6]  //6 O_NH4 
#define K_NH4ads_2_NH4 parms[7] //8 R_NH4a
#define K_NO3_2_NH4 parms[8] //7 D_NO3


/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
  int N= 9 ;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
  //----------Full
  double NLab_2_NH4,NH4_2_NLab,NRec_2_NH4,NRec_2_NO3,NO3_2_NRec,NH4ads_2_NH4,NO3_2_NH4,NH4_2_NO3,NO3sto_2_NO3;
  double NLab,NH4,NO3,NH4ads,NRec,NO3sto;
  double dNLab,dNH4,dNO3,dNH4ads,dNRec,dNO3sto;  
 
  //-----------N15
  double NLab_2_NH4_15,NH4_2_NLab_15,NRec_2_NH4_15,NRec_2_NO3_15,NO3_2_NRec_15,NH4ads_2_NH4_15,NO3_2_NH4_15,NH4_2_NO3_15,NO3sto_2_NO3_15;
  double NLab_15,NH4_15,NO3_15,NH4ads_15,NRec_15,NO3sto_15;
  double dNLab_15,dNH4_15,dNO3_15,dNH4ads_15,dNRec_15,dNO3sto_15;  
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
  
  NO3sto = y[10];
  NO3sto_15 =y[11];
  
  
//----full  
  NLab_2_NH4 = K_NLab_2_NH4*NLab;
  NH4_2_NLab = C_NH4_2_NLab;
  
  dNLab = NH4_2_NLab - NLab_2_NH4; //1
  
  NRec_2_NH4 = C_NRec_2_NH4;
  NRec_2_NO3 = C_NRec_2_NO3;
  NO3_2_NRec = C_NO3_2_NRec;
  dNRec = NO3_2_NRec - NRec_2_NH4 - NRec_2_NO3;
  
  NH4ads_2_NH4 = K_NH4ads_2_NH4 * NH4ads;
  NO3_2_NH4 = K_NO3_2_NH4 * NO3;
  NH4_2_NO3 = K_NH4_2_NO3 * NH4;
  dNH4 = NLab_2_NH4 + NRec_2_NH4 + NH4ads_2_NH4 +NO3_2_NH4 - NH4_2_NLab -NH4_2_NO3;
  
  dNH4ads = -NH4ads_2_NH4;
  
  dNO3 = NRec_2_NO3 + NO3sto_2_NO3 + NH4_2_NO3 - NO3_2_NH4 - NO3_2_NRec;
  dNO3sto = -NO3sto_2_NO3;
  //---N15
  
  NLab_2_NH4_15 = K_NLab_2_NH4*NLab_15;
  NH4_2_NLab_15 = C_NH4_2_NLab;
  
  dNLab_15 = NH4_2_NLab_15 - NLab_2_NH4_15; //1
  
  NRec_2_NH4_15 = C_NRec_2_NH4;
  NRec_2_NO3_15 = C_NRec_2_NO3;
  NO3_2_NRec_15 = C_NO3_2_NRec;
  dNRec_15 = NO3_2_NRec_15 - NRec_2_NH4_15 - NRec_2_NO3_15;
  
  NH4ads_2_NH4_15 = K_NH4ads_2_NH4 * NH4ads_15;
  NO3_2_NH4_15 = K_NO3_2_NH4 * NO3_15;
  NH4_2_NO3_15 = K_NH4_2_NO3 * NH4_15;
  dNH4_15 = NLab_2_NH4_15 + NRec_2_NH4_15 + NH4ads_2_NH4_15 +NO3_2_NH4_15 - NH4_2_NLab_15 -NH4_2_NO3_15;
  
  dNH4ads_15 = 0;
  
  dNO3_15 = NRec_2_NO3_15 + NO3sto_2_NO3_15 + NH4_2_NO3_15 - NO3_2_NH4_15 - NO3_2_NRec_15;
  dNO3sto_15 = 0;
  
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
  
  ydot[10] =dNO3sto;
  ydot[11]= dNO3sto_15;
  
  yout[0] = NH4;
  yout[1] = NO3;
  yout[2] = NH4_15;
  yout[3] = NO3_15;
}