library("FME")
library('compiler')
library('multicore')

dyn.load('flauz.dll')
loadcmp('flauz_subs.Rc')

# mdMCMC_MC<-function(ff=FRAUZcost1,pp=pp,obs=obs,niter=50000,y=y0,times,wvar0=0.05,updatecov=200,lower=lower,upper=upper,burninlength=10000){
	# set.seed(1000*rnorm(1))
# #	Fit <- modFit(f=ff,p=pp,y=y,times=times,obs=obs,lower=lower,upper=upper,method='SANN');
	# Fit <- modFit(f=ff,p=Fit$par,y=y,times=times,obs=obs,lower=lower,upper=upper,method='Pseudo');
	# pp=Fit$par*SS;
	# MCMC<-mdMCMC1(ff=ff,pp=pp,obs=obs,niter=niter,y=y,times,wvar0=wvar0,updatecov=updatecov,lower=lower,upper=upper,burninlength=burninlength)
	# return(MCMC)
# }

# mdMCMC_SS<-function(ff=FRAUZcost,pp=ppSS,obs=obs,niter=50000,y=y0,SS=SS,times=times,wvar0=0.1,updatecov=200,lower=lower,upper=upper,burninlength=10000){
	# set.seed(1000*rnorm(1))
	# pp=pp*SS;
# #	Fit <- modFit(f=ff,p=pp,y=y,times=times,obs=obs,lower=lower,upper=upper,method='SANN');
	# Fit <- modFit(f=ff,p=pp,y=y,times=times,obs=obs,lower=lower,upper=upper,method='Pseudo');
	# pp=Fit$par*SS;
	# MCMC<-mdMCMC(ff=ff,pp=pp,obs=obs,niter=niter,y=y,times=times,SS=SS,wvar0=wvar0,updatecov=updatecov,lower=lower,upper=upper,burninlength=burninlength)
	# return(MCMC)
# }

#######################################################################################################################
#------------------Inti STATE VARIABLES
y0 = c(NH4 = 44.45,
	   NO3 = 25.43,
	   NH4_N15 = 4.48,
	   NO3_N15 = 0.033,
	   TON_N15=0.049,
	   TON = 37.0,	
	   Mib = 0,
	   PR = 17,
	   Hum = 20,
	   Mib_N15 = 0.049,
	   Hum_N15 = 0.049,
	   NH3 = 0.000,
	   N2O = 0.000,
	   NH3_N15 = 0.000,
	   N2O_N15 = 0.000,
	   PR_N15 = 0
		   );
#--------------Parameters
pars0 = c(C_hu_nh4 = 0.51, 
		C_micb_hum= 0., 
		C_pr_nh4= 0. , 
		C_pr_micb= 0.0001,
		C_micb_nh4=0.1,
		K_nh4mic = 0.5,
		K_no3mic= 0.1,
		K_nh4no3= 0.4,
		K_no3ng= 0.0001,
		K_nh4nh3= 0.0001,
		K_no3nh4=0
		);

#----------Switch to on/off some processes,the values must be 0 or 1.
SS =  c(S_hu_nh4 = 1, 
		S_micb_hum= 1, 
		S_pr_nh4= 0, 
		S_pr_micb= 0,
		S_micb_nh4=1,
		S_nh4mic = 1,
		S_no3mic= 1,
		S_nh4no3= 1,
		S_no3ng= 0,
		S_nh4nh3= 0,
		S_no3nh4= 0
		);
 
#KEY question now become the START points!!!!!!!!!! 
#Read data into memory, variables must be some of the following variables:
#-----------------------------------------------------------------------------------------------------------------------------------------------
#Time, NH4, Err_NH4, NO3, Err_NO3, NH4_N15, Err_NH4_N15, NO3_N15, Err_NO3_N15, ON, Err_ON, ON_N15, Err_ON_N15, N2O,Err_N2O, N2O_N15, Err_N2O_N15
#-----------------------------------------------------------------------------------------------------------------------------------------------

obs<-read.csv('Data_1.csv')
names(obs)<-c('time', 'NH4', 'NO3', 'NH4_N15', 'NO3_N15','TON_N15')
obs[,'NH4_N15']<-obs[,'NH4']*obs[,'NH4_N15']/100
obs[,'NO3_N15']<-obs[,'NO3']*obs[,'NO3_N15']/100
obs[,'TON_N15']<-1307*obs[,'TON_N15']/100

times<- obs[,'time']

inrep = 300;
bestpar<-data.frame(matrix(-9999,ncol=12,nrow=inrep));
names(bestpar)<-c('CAT',names(pars0))
row.names<-rep(c('MCMC','MCMC_SS'),(inrep/2))
bestpar[,1]<-rep(c('MCMC','MCMC_SS'),(inrep/2))

### MCMC simulation
#----------------------------------------------------------------------------
pp<-pars0;
ppSS<-pars0*SS;
upper<-lower<-pars0;
lower[]= 0;
upper[]= 1;


ik = seq(1,inrep,6);

for (i in ik){
set.seed(1000*rnorm(1))

Fit <- modFit(f=FRAUZcost1,p=pp,y=y0,times=times,obs=obs,lower=lower,upper=upper,method='Pseudo');
pp = Fit$par
Fit <- modFit(f=FRAUZcost,p=pp,y=y0,SS=SS,times=times,obs=obs,lower=lower,upper=upper,method='Pseudo');
ppSS = Fit$par*SS;

print(system.time({
MC1<-parallel(mdMCMC_MC(ff=FRAUZcost1,pp=pp,obs=obs,niter=30000,y=y0,times=times,wvar0=0.05,updatecov=200,lower=lower,upper=upper,burninlength=10000),name='MC1');
MC2<-parallel(mdMCMC_MC(ff=FRAUZcost1,pp=pp,obs=obs,niter=30000,y=y0,times=times,wvar0=0.05,updatecov=200,lower=lower,upper=upper,burninlength=10000),name='MC2');
MC3<-parallel(mdMCMC_MC(ff=FRAUZcost1,pp=pp,obs=obs,niter=30000,y=y0,times=times,wvar0=0.05,updatecov=200,lower=lower,upper=upper,burninlength=10000),name='MC3');

# MC2<-parallel(MCMC<-mdMCMC1(ff=FRAUZcost1,pp=pp,obs=obs,niter=25000,y=y0,times,wvar0=0.05,updatecov=50,lower=lower,upper=upper,burninlength=5000),name='MC2');
# MC3<-parallel(MCMC<-mdMCMC1(ff=FRAUZcost1,pp=pp,obs=obs,niter=30000,y=y0,times,wvar0=0.05,updatecov=50,lower=lower,upper=upper,burninlength=5000),name='MC3');
MS1<-parallel(MCSS1<-mdMCMC_SS(ff=FRAUZcost,pp=ppSS,obs=obs,niter=30000,y=y0,SS=SS,times=times,wvar0=0.1,updatecov=200,lower=lower,upper=upper,burninlength=10000),name='MCS1');
MS2<-parallel(MCSS2<-mdMCMC_SS(ff=FRAUZcost,pp=ppSS,obs=obs,niter=30000,y=y0,SS=SS,times=times,wvar0=0.1,updatecov=200,lower=lower,upper=upper,burninlength=10000),name='MCS2');
MS3<-parallel(MCSS3<-mdMCMC_SS(ff=FRAUZcost,pp=ppSS,obs=obs,niter=30000,y=y0,SS=SS,times=times,wvar0=0.1,updatecov=200,lower=lower,upper=upper,burninlength=10000),name='MCS3');

fq<-collect(list(MC1,MC2,MC3,MS1,MS2,MS3),wait=TRUE);
#fq<-collect(list(MC,MC2,MC3),wait=TRUE);
}))
MCMC1<-fq[[1]];
MCMC2<-fq[[2]];
MCMC3<-fq[[3]];

MCMC_SS1<-fq[[4]];
MCMC_SS2<-fq[[5]];
MCMC_SS3<-fq[[6]];
# MCMC2<-fq[[2]];
# MCMC3<-fq[[3]];
### MCMC output
#------------------------------------------------------------------------------
#print('----MCMC ---------------------------------------------------------------')
#print(summary(MCMC))
# print(summary(MCMC2))
# print(summary(MCMC3))#
#print('----MCMC_SS ---------------------------------------------------------------')
#print(summary(MCMC_SS))
### MCMC best parameters
#------------------------------------------------------------------------------
bestpar[c(i:(i+5)),2:12]<-rbind(MCMC1$bestpar,MCMC2$bestpar,MCMC3$bestpar,MCMC_SS1$bestpar,MCMC_SS2$bestpar,MCMC_SS3$bestpar)
#row.names(bestpar)<-c('MCMC','MCMC_SS')
#print(bestpar)
#print(MCMC$bestpar)
# print(MCMC2$bestpar)
# print(MCMC3$bestpar)
#print(MCMC_SS$bestpar)
###Plot the processes
#-----------------------------------------------------------------------------
#pdf('cumuplot.pdf')
#cumuplot(as.mcmc(MCMC$pars))
#dev.off()
# ### Plot the state variables
# #-----------------------------------------------------------------------------
 # sim_data<-as.data.frame(FRAUZ(y0,MCMC_SS$bestpar,SS,times));
 # pv2<-valid_plotMCMCpdf(obs,sim_data,'MCMC_SS.pdf');
  # sim_data<-as.data.frame(FRAUZ1(y0,MCMC$bestpar,times));
  # pv2<-valid_plotMCMCpdf(obs,sim_data,'MCMC.pdf');
# # # ### Sensitivity ranges
# # #------------------------------------------------------------------------------
 # arg_g=c('NH4','NO3','NH4_N15','NO3_N15','TON_N15','Mib','Hum','NH3','N2O')
 # sR1 <- sensRange(func = FRAUZ1, parms = MCMC$bestpar, times=times, y=y0, parInput = MCMC$par, dist = "latin",sensvar=arg_g)
	# pdf('SR1.pdf')
	 # opar<-par();
	 # par(mfcol=c(4,4), mar=c(3,3,2,0.5))
	 # plot(summary(sR1), xlab = "time")
	 # par(opar);
	 # dev.off();
 # sR2 <- sensRange(func = FRAUZ, parms = MCMC_SS$bestpar, times=times, y=y0, SS=SS, parInput = MCMC_SS$par, dist = "latin",sensvar=arg_g)
	# pdf('SR2.pdf')
	 # opar<-par();
	 # par(mfcol=c(4,4), mar=c(3,3,2,0.5))
	 # plot(summary(sR2), xlab = "time")
	 # par(opar);
	 # dev.off();
	print(i) 
}	

save(bestpar,file='bestpar.RData') 
#-----------------------------------------------------------------------------
dyn.unload('flauz.dll')
