clear
close all
clc;

%% figure formatting

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter', 'latex')
   
set(0,'DefaultTextFontSize', 14)
set(0,'DefaultAxesFontSize', 14)
set(0,'DefaultLineLineWidth',2)

temp = get(gca,'ColorOrder');
c1 = temp(1,:);
c2 = temp(2,:);

close all

%% set parameters

% discretize Gaussian distribution by Gauss-Hermite quadrature
S = 7; % number of states
[x, w] = GaussHermite(S);
PR = w/sqrt(pi); % probability
x = sqrt(2)*x;

delta = 0.04; % discount rate
eta = 0.025; % death rate
D = 1/4; % length of one period
v = exp(-eta*D); % survival probability
beta = exp(-delta*D); % monthly discount factor
gamma = 2; % risk aversion
mu = 0.07; % expected return
sigma = 0.15;
R = exp((mu-sigma^2/2)*D + sigma*sqrt(D)*x); % realized asset return

rho = (beta*dot(PR,R.^(1-gamma)))^(1/gamma)
if rho < 1
    disp('MPC is positive');
else
    disp('MPC is zero');
end

alphaY = 3; % income Pareto exponent
yProm = 5; % years to get promoted on average
pProm = 1 - exp(-D/yProm);
PG = [1-pProm pProm]';
g = log((1-v+v*pProm)/(v*pProm))/alphaY; % income growth when promotion occurs
G = [1 exp(g)]';


N = 100;
ggrid = linspace(0.02,0.1,N)';

alpha_capital = zeros(N,1);
alpha_labor = zeros(N,1);

for n=1:N
    alpha_labor(n) = log((1-v+v*pProm)/(v*pProm))/ggrid(n);
    Gn = [1 exp(ggrid(n))]';
    Rtilde = kron(1./Gn,R);
    P = kron(PG,PR);
    lambda = @(z)(v*dot(P,(rho*Rtilde).^z));
    alpha_capital(n) = fzero(@(z)(log(lambda(z))),1);
end

figure
plot(100*ggrid,alpha_labor,100*ggrid,alpha_capital);
xlabel('Income growth rate when promoted (\%)')
ylabel('Pareto exponent')
legend('Income','Normalized wealth')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_IFP_alpha','-dpdf')

%% solve detrended income fluctuation problem

% compute wealth Pareto exponent
Rtilde = kron(1./G,R);
P = kron(PG,PR);
lambda = @(z)(v*dot(P,(rho*Rtilde).^z));
Ngrid = 100;
zGrid = linspace(0,alphaY,Ngrid);
fval = 0*zGrid;
for n=1:Ngrid
    fval(n) = lambda(zGrid(n));
end

figure
plot(zGrid,fval)

alpha = fzero(@(z)(log(lambda(z))),1);

% define detrended discount factor
Gtilde = kron(G,ones(S,1));
betatilde = beta*Gtilde.^(1-gamma);

MPC = 1 - (PR'*(beta.*R.^(1-gamma)))^(1/gamma); % asymptotic MPC;

% define asset grid

aMed = 10;
aMax = aMed*1e3;
aGrid = expGrid(0,aMax,aMed,Ngrid);

MaxIter = 1000;
tol = 1e-6;

tic
ca = getC(betatilde,gamma,P,Rtilde,aGrid,MaxIter,tol);
toc

aBar = 100;
caBar = interp1(aGrid,ca,aBar);
tangent = caBar + MPC*(aGrid - aBar);
MPCerror = ((ca(end) - ca(end-1))/(aGrid(end) - aGrid(end-1)))/MPC - 1

figure
plot(aGrid,ca); hold on
%plot(aGrid,tangent,'--','Color',c1)
xlim([0 aBar])
xlabel('Normalized wealth')
ylabel('Normalized consumption')
%legend('$\tilde{c}(a)$','Asymptotic slope')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_IFP_ca','-dpdf')

%% compute stationary wealth distribution
aNextMat = kron(Rtilde',(aGrid - ca)+ 1); % next periods' asset
[Q,pi] = getQ(1,P',1-v,1,aGrid,aNextMat,[],alpha); % use Pareto extrapolation
cCDF = 1 - cumsum(pi);
aGridTemp = linspace(aGrid(end),1e6,Ngrid);
cCDFTemp = cCDF(end-1)*(aGridTemp/aGrid(end-1)).^(-alpha);

figure
loglog(aGrid(1:end-1),cCDF(1:end-1)); hold on
loglog(aGridTemp,cCDFTemp,'--','Color',c1)
xlabel('Normalized wealth')
ylabel('Tail probability')
legend('Model','Extrapolated')
xlim([1e-2,1e6])

%% simulation

rng(1); % set seed for replicability

F = griddedInterpolant(aGrid,ca,'linear','linear'); % linear interpolation/extrapolation

I = 1e5; % number of agents
T = ceil(1e3/D); % number of periods

Rsim = exp((mu-sigma^2/2)*D + sigma*sqrt(D)*normrnd(0,1,[I,T])); % realized asset return
Vsim = binornd(1,v,[I,T]); % indicator of survival
Prom = binornd(1,pProm,[I,T]); % indicator of promotion
Gsim = (1-Prom) + exp(g)*Prom;
Rtildesim = Rsim./Gsim;

atilde = ones(I,1);
Y = exprnd(1/alphaY,[I,1]);
for t=1:T
    atilde = (1-Vsim(:,t)) + Vsim(:,t).*(Rtildesim(:,t).*(atilde - F(atilde)) + 1); % update normalized wealth
    Y = (1-Vsim(:,t)) + Vsim(:,t).*Gsim(:,t).*Y; % update income
end
a = Y.*atilde;

figure
histogram(log(Y))

figure
histogram(log(a))

atilde_sort = sort(atilde,'descend');
Y_sort = sort(Y,'descend');
Ygrid = exp(linspace(0,log(100),100));
cCDFY = Ygrid.^(-alphaY);

[CDF_Y,Yval] = ecdf(Y);
[CDF_atilde,atildeval] = ecdf(atilde);

figure
plot(Ygrid,cCDFY,'--','Color',c1); hold on
stairs(Yval,1-CDF_Y,'-','LineWidth',1,'Color',c1);
%loglog(Y_sort,[1:I]/I,'o','Color',c1);
plot(aGrid(1:end-1),cCDF(1:end-1),'--','Color',c2);
stairs(atildeval,1-CDF_atilde,'-','LineWidth',1,'Color',c2);
%loglog(atilde_sort,[1:I]/I,'o','Color',c2);
hold off
set(gca,'XScale','log','YScale','log');
xlabel('Size')
ylabel('Tail probability')
legend('$Y$ (theory)','$Y$ (simulation)','$a/Y$ (theory)','$a/Y$ (simulation)')
xlim([1,1e4])
ylim([1e-5,1])

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_IFP_sim','-dpdf')