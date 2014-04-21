%%flauz.m
function ydot = flauz(t,y0,parms)

%%Parameters
C_hu_nh4 = parms(1);
C_micb_hum = parms(2);
C_pr_nh4 = parms(3);
C_pr_micb = parms(4);
C_micb_nh4 = parms(5);
K_nh4mic = parms(6);
K_no3mic = parms(7);
K_nh4no3 = parms(8);
K_no3ng = parms(9);
K_nh4nh3 = parms(10);
K_no3nh4 = parms(11);

%%initilizal conditions
NH4 = y0(1);
NO3 = y0(2);
TON = y0(10);
	
NH4_N15 = y0(3);
NO3_N15 = y0(4);
	
Mib = y0(6);
PR = y0(7);
Hum = y0(8);

Mib_N15 = y0(9);
TON_N15 = y0(5);

%%calculation
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
dn = K_no3nh4 * NO3;					
%%//#assume conversion from micb to humus will take very long time, then the h,s,j and m for 15N balance was omitted
					
v15 = K_nh4nh3 * NH4_N15;
n15 = K_nh4no3 * NH4_N15;
i_a_15 = K_nh4mic * NH4_N15;
i_n_15 = K_no3mic * NO3_N15;
r15 = C_micb_nh4;	
d15 = K_no3ng * NO3_N15;
dn15= K_no3nh4 * NO3_N15;

dNH4 = m + s - n - i_a - v + dn;
dNO3 = n - d - i_n - dn;
dMib = i_a + i_n +j - h - r;
dPR = - s - j;
dHum = h - m;
dTON = dMib + dPR + dHum;
					
dNH4_N15 = (-n15 - i_a_15 - v15 + r15 + dn15);
dNO3_N15 =  (n15 - i_n_15 - d15 - dn15);
dMib_N15 = (i_a_15 + i_n_15 - r15);
dHum_N15 = h - m;
dTON_N15 = (dMib_N15 + dHum_N15);

%%ydot 
%
ydot = [dNH4;dNO3;dNH4_N15;dNO3_N15;dTON_N15;dTON;dMib;dPR;dHum;dMib_N15;dHum_N15;dTON_N15];
