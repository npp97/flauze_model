function ss = flauz_ss(param,data)
% sum-of-squares for Himmelblau 9.9
% logLike = 0;

time = data.ydata(:,1);
Aobs = data.ydata(:,2:6);
% [nrow,col] = size(Aobs);
y0   = data.y0;
[t,y] = ode45(@flauz,time,y0,[],param);
Amodel = y(:,1:5);
ss = sum(sum((Aobs-Amodel).^2));

% sigma = sqrt(ss/col);
% logLike= logLike+ss/(2*sigma.^2);
% logLike= logLike+col*log(sigma);
% logLike= logLike+col*sqrt(2*pi);
% ss = logLike;
  
