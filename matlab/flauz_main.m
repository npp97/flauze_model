clear model data parama options

data.ydata = [
%  Time (days)     [A] (mole/liter)
%time	NH4	    NO3	    NH4_N15	NO3_N15	TON_N15
0.0 	44.45	25.43	9.504	0.013	0.0025
3.0	    23.45	39.6	8.966	2.811	0.0691
5.0	    6.23	48.6	7.604	4.184	0.1322
9.0	    0.65	39.92	2.187	4.45	0.1842
13.0	0.93	34.33	1.992	4.236	0.2107
16.0	0.93	32.11	1.756	4.239	0.2369
21.0	1.00	31.25	1.819	4.193	0.2335
28.0	1.07	31.83	1.688	3.428	0.2534
 ];

data.ydata(:,4) = data.ydata(:,4).* data.ydata(:,2)./100;
data.ydata(:,5) = data.ydata(:,5).* data.ydata(:,3)./100;
data.ydata(:,6) = 1037.* data.ydata(:,6)./100;

%%
% Initial concentrations are saved in |data| to be used in sum of
% squares function.

NH4=44.45;NO3 = 25.43;TON = 1037.0;NH4_N15 = 4.22;
NO3_N15 = 0.0033;TON_N15=0.026;Mib = 0; PR = 500; Hum =537.0; Mib_N15=0;Hum_N15=0.026;
data.y0 = [NH4;NO3;NH4_N15;NO3_N15;TON_N15;TON;Mib;PR;Hum;Mib_N15;Hum_N15;TON_N15];

C_hu_nh4 = 0.01; C_micb_hum = 0.01; C_pr_nh4 = 0.01;C_pr_micb = 0.01;
C_micb_nh4 = 0.03;K_nh4mic = 0.06; K_no3mic= 0.06; K_nh4no3= 0.18;
K_no3ng = 0.01; K_nh4nh3 = 0.01;K_no3nh4 =0.001;

p00 = [C_hu_nh4;C_micb_hum;C_pr_nh4;C_pr_micb;C_micb_nh4;K_nh4mic;K_no3mic;K_nh4no3;K_no3ng;K_nh4nh3;K_no3nh4];

%%
% Refine the first guess for the parameters with |fminseacrh| and
% calculate residual variance as an estimate of the model error variance.
[p0,ss0] = fminsearch(@flauz_ss,p00,[],data);
mse = ss0/(length(data.ydata)-4);
p0=abs(p0);
%p0=p00;
%%
params = {
    {'C_hu_nh4', p0(1), 0, 1,0}
    {'C_micb_hum', p0(2),0, 1,0}
    {'C_pr_nh4', p0(3),0, 1,0}
    {'C_pr_micb',p0(4),0, 1,0}
    {'C_micb_nh4',p0(5),0, 1,0}
    {'K_nh4mic',p0(6),0, 1,0}
    {'K_no3mic',p0(7),0, 1,0}
    {'K_nh4no3',p0(8),0, 1,0}
    {'K_no3ng',p0(9),0, 1,0}
    {'K_nh4nh3',p0(10),0, 1,0}
    {'K_no3nh4',p0(11),0, 1,0}
    };

model.ssfun = @flauz_ss;
model.sigma2 = mse;
%options.method = 'am';

options.nsimu = 3000;
options.updatesigma = 1;

[results,chain,s2chain] = mcmcrun(model,data,params,options);
aresult = chainstats(chain,results);

p0 = aresult(:,1);
lo = aresult(:,1) - 2.* aresult(:,2);
uu = aresult(:,1) + 2.* aresult(:,2);
lo(:)=0;
tha = results.theta;
params = {
    {'C_hu_nh4', p0(1), lo(1), uu(1),tha(1)}
    {'C_micb_hum', p0(2),lo(2), uu(2),tha(2)}
    {'C_pr_nh4', p0(3),lo(3), uu(3),tha(3)}
    {'C_pr_micb',p0(4),lo(4),uu(4),tha(4)}
    {'C_micb_nh4',p0(5),lo(5), uu(5),tha(5)}
    {'K_nh4mic',p0(6),lo(6), uu(6),tha(6)}
    {'K_no3mic',p0(7),lo(7), uu(7),tha(7)}
    {'K_nh4no3',p0(8),lo(8), uu(8),tha(8)}
    {'K_no3ng',p0(9),lo(9), uu(9),tha(9)}
    {'K_nh4nh3',p0(10),lo(10), uu(10),tha(10)}
    {'K_no3nh4',p0(11),lo(11),uu(11),tha(11)}
    };
options.nsimu = 20000;
options.qcov = results.cov;
[results,chain,s2chain] = mcmcrun(model,data,params,options);
chainstats(chain,results)
%%
figure(1); clf
mcmcplot(chain,[],results,'chainpanel')
subplot(2,2,4)
mcmcplot(sqrt(s2chain),[],[],'dens',2)
title('error std')

%%
% Function |chainstats| lists some statistics, including the
% estimated Monte Carlo error of the estimates.

%%
figure(2); clf
[t,y] = ode45(@flauz,linspace(0,30),data.y0,[],mean(chain));
title('Data and fitted model')
subplot(2,3,1);plot(data.ydata(:,1),data.ydata(:,2),'s',t,y(:,1),'-','LineWidth',2);title('NH4');
subplot(2,3,2);plot(data.ydata(:,1),data.ydata(:,3),'s',t,y(:,2),'-','LineWidth',2);title('NO3');
subplot(2,3,3);plot(data.ydata(:,1),data.ydata(:,4),'s',t,y(:,3),'-','LineWidth',2);title('NH4-N15');
subplot(2,3,4);plot(data.ydata(:,1),data.ydata(:,5),'s',t,y(:,4),'-','LineWidth',2);title('NO3-N15');
subplot(2,3,5);plot(data.ydata(:,1),data.ydata(:,6),'s',t,y(:,5),'-','LineWidth',2);title('TON-N15');
subplot(2,3,6);plot(t,y(:,6:9),'-',t,y(:,10:11).*10,'LineWidth',2);title('Other');
legend('TON','Mib','PR','Hum','Mib_N15','Hum_N15',0);