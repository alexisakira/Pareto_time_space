close all
clc;
clear

load testData

[Q,dist] = getQ(P,pi',0,[],Gw,wprimeMat,[],zeta);
wsDist = reshape(dist,[length(Gw) S])';
wDist = sum(wsDist,1);
wDistCDF = cumsum(wDist);

[Q_old,dist_cont] = getQ_cont(P,pi',0,[],Gw,wprimeMat,wprimeSlopeBar,zeta);
wsDist_cont = reshape(dist_cont,[length(Gw) S])';
wDist_cont = sum(wsDist_cont,1);
wDistCDF_cont = cumsum(wDist_cont);

[Q_temp,dist_temp] = getQ_temp(P,pi',0,[],Gw,wprimeMat,wprimeSlopeBar,zeta);
wsDist_temp = reshape(dist_temp,[length(Gw) S])';
wDist_temp = sum(wsDist_temp,1);
wDistCDF_temp = cumsum(wDist_temp);

figure
loglog(Gw,1-wDistCDF,Gw,1-wDistCDF_cont)
xlim([1 Gw(end-1)])

figure
loglog(Gw,(Gw.^zeta).*(1-wDistCDF),Gw,(Gw.^zeta).*(1-wDistCDF_cont))
xlim([1 Gw(end-1)])

figure
loglog(Gw,1-wDistCDF,Gw,1-wDistCDF_cont,Gw,1-wDistCDF_temp)
xlim([1 Gw(end-1)])

figure
loglog(Gw,(Gw.^zeta).*(1-wDistCDF),Gw,(Gw.^zeta).*(1-wDistCDF_cont),Gw,(Gw.^zeta).*(1-wDistCDF_temp))
xlim([1 Gw(end-1)])

figure
loglog(Gw,(Gw.^zeta).*(1-wDistCDF_temp))
xlim([1 Gw(end-1)])
