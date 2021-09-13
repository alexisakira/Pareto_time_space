% create figures and tables in the paper

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
c3 = temp(3,:);

close all

Nmin = 1000; % minimum sample size to include in analysis
%% import data

rule = '5'; % 5% rule
%rule = '1'; % 1% rule

filename = 'New_wave_LIS.xlsx';
sheet = ['Labor capital Hill ' rule '% rule'];
sheetHoga = ['Test Labor Capital inequal ' rule '%'];

% labor income
[~,country] = xlsread(filename,sheet,'A4:A999');
year = xlsread(filename,sheet,'B4:B999');
alphaL = xlsread(filename,sheet,'C4:C999');
seL = xlsread(filename,sheet,'D4:D999');
NL = xlsread(filename,sheet,'F4:F999');

% capital income
alphaK = xlsread(filename,sheet,'G4:G999');
seK = xlsread(filename,sheet,'H4:H999');
NK = xlsread(filename,sheet,'J4:J999');

%% summary statistics

N_country = length(unique(country))
min(year)
max(year)

if min(NL) < Nmin
    warning('minimum sample size below Nmin; calculation inaccurate')
end
mean(alphaL)
med_lab = median(alphaL)
std(alphaL)

ind = find(NK >= Nmin);
NK_incl = length(ind)
mean(alphaK(ind))
med_cap = median(alphaK(ind))
std(alphaK(ind))

%% histogram

figure
histogram(alphaL,'normalization','pdf'); hold on
histogram(alphaK(ind),'normalization','pdf');
xline(med_lab,'k--');
xline(med_cap,'k--');
hold off
xlabel('Pareto exponent')
ylabel('Probability density')
legend('Labor income','Capital income','Median')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_alphaHist' rule],'-dpdf')

%% scatter plot

country_unique = unique(country);
C = length(country_unique); % number of countries

[R,P,RL,RU] = corrcoef(alphaL(ind),alphaK(ind));
rho = R(1,2);
rhoL = RL(1,2);
rhoU = RU(1,2);

aU = ceil(max([alphaL(ind); alphaK(ind)]));

% scatter plot of all data
figure
scatter(alphaL(ind),alphaK(ind)); hold on
plot([0 aU],[0 aU],'k-','Linewidth',1);
%str = ['$\rho= $' num2str(round(rho,2,'significant')) ' (95\%CI: ['...
%    num2str(round(rhoL,2,'significant')) ', ' num2str(round(rhoU,2,'significant')) '])'];
%T = text(0.1,aU-2,str);
%set(T, 'fontsize', 14, 'verticalalignment', 'middle', 'horizontalalignment', 'left');
text(0.1,aU-1,['$\rho= $' num2str(round(rho,2,'significant'))],'VerticalAlignment','top');
text(0.1,aU-2,['(95\%CI: [' num2str(round(rhoL,2,'significant')) ', ' num2str(round(rhoU,2,'significant')) '])'],...
    'VerticalAlignment','bottom');
axis equal
xlim([0 aU])
ylim([0 aU])
xlabel('Labor Pareto exponent')
ylabel('Capital Pareto exponent')
legend('Pareto exponent','$45^\circ$ line','Location','NW')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_alphaScatter' rule],'-dpdf')

%% scatter plot of 5% and 1% rules

alphaL5 = xlsread(filename,'Labor capital Hill 5% rule','C4:C999');
alphaL1 = xlsread(filename,'Labor capital Hill 1% rule','C4:C999');
alphaK5 = xlsread(filename,'Labor capital Hill 5% rule','G4:G999');
alphaK1 = xlsread(filename,'Labor capital Hill 1% rule','G4:G999');
indL = find(NL >= Nmin);
indK = find((NK >= Nmin)&(~isnan(alphaK5))&(~isnan(alphaK1)));

alphaL5 = alphaL5(indL);
alphaL1 = alphaL1(indL);
alphaK5 = alphaK5(indK);
alphaK1 = alphaK1(indK);

% labor income
[R,P,RL,RU] = corrcoef(alphaL5,alphaL1);
rhoL = R(1,2);

aU = ceil(max([alphaL5; alphaL1]));

figure
scatter(alphaL5,alphaL1); hold on
plot([0 aU],[0 aU],'k-','LineWidth',1);
text(0.1,aU-2,['$\rho= $' num2str(round(rhoL,2,'significant'))],'VerticalAlignment','top');
axis equal
xlabel('Pareto exponent for 5\% rule')
ylabel('Pareto exponent for 1\% rule')
legend('Pareto exponent','$45^\circ$ line','Location','NW')
xlim([0,aU])
ylim([0,aU])

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_alphaLcomp','-dpdf')

% capital income
[R,P,RL,RU] = corrcoef(alphaK5,alphaK1);
rhoK = R(1,2);

aU = ceil(max([alphaK5; alphaK1]));

figure
scatter(alphaK5,alphaK1); hold on
plot([0 aU],[0 aU],'k-','LineWidth',1);
text(0.1,aU-2,['$\rho= $' num2str(round(rhoK,2,'significant'))],'VerticalAlignment','bottom');
axis equal
xlabel('Pareto exponent for 5\% rule')
ylabel('Pareto exponent for 1\% rule')
legend('Pareto exponent','$45^\circ$ line','Location','NW')
xlim([0,aU])
ylim([0,aU])

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_alphaKcomp','-dpdf')

%% time series for median exponents across countries

Tobs = zeros(C,1);
alphaLt = alphaL(ind);
alphaKt = alphaK(ind);
countryt = country(ind);
yeart = year(ind);

for c=1:C
    indc = find(strcmp(countryt,country_unique{c}));
    Tobs(c) = length(indc);
end

Tmin = 8; % minimum number of observations
ind_incl = find(Tobs >= Tmin);
country_unique(ind_incl)
K = length(ind_incl);

time = [min(yeart):max(yeart)];
T = length(time);
AlphaL = zeros(K,T); % store interpolated Pareto exponents
AlphaK = zeros(K,T);
for k=1:K
    ind = find(strcmp(countryt,country_unique{ind_incl(k)}));
    if length(ind) == 1
        AlphaL(k,:) = alphaLt(ind);
        AlphaK(k,:) = alphaKt(ind);
    else % linear interpolation
        temp = interp1(yeart(ind),alphaLt(ind),time);
        if any(isnan(temp)) % needs extrapolation
            ind1 = find(~isnan(temp),1,'first');
            ind2 = find(~isnan(temp),1,'last');
            AlphaL(k,:) = temp;
            AlphaL(k,1:ind1) = temp(ind1);
            AlphaL(k,ind2:end) = temp(ind2);
        end
        temp = interp1(yeart(ind),alphaKt(ind),time);
        if any(isnan(temp)) % needs extrapolation
            ind1 = find(~isnan(temp),1,'first');
            ind2 = find(~isnan(temp),1,'last');
            AlphaK(k,:) = temp;
            AlphaK(k,1:ind1) = temp(ind1);
            AlphaK(k,ind2:end) = temp(ind2);
        end
    end
end

AlphaL_med = median(AlphaL,1); % median
AlphaK_med = median(AlphaK,1);

figure
plot(time,AlphaL_med,time,AlphaK_med);

xlabel('Year')
ylabel('Median Pareto exponent')
legend('Labor','Capital','Location','Best')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_alphaTS' rule],'-dpdf')

%% time series for selected countries

countryc = {'Canada','Germany','Switzerland','Taiwan','United Kingdom','United States'};
abbrev = {'CAN','GER','SWI','TW','UK','US'};
K = length(countryc);

for k = 1:K
    indc = find(strcmp(country,countryc{k})&(NK >= Nmin));
    yearc = year(indc);
    alphaLc = alphaL(indc);
    seLc = seL(indc);
    alphaKc = alphaK(indc);
    seKc = seK(indc);

    CIL = alphaLc + 1.96*seLc*[-1 1];
    CIK = alphaKc + 1.96*seKc*[-1 1];

    figure
    plot(yearc,alphaLc,'-','Color',c1); hold on
    plot(yearc,alphaKc,'-','Color',c2);
    plot(yearc,CIL,':','Color',c1);
    plot(yearc,CIK,':','Color',c2); hold off
    xlabel('Year')
    ylabel('Pareto exponent')
    legend('Labor','Capital','95\% CI','Location','Best')
    title(countryc{k})

    %save figure in pdf format
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,['fig_alphaTS_' abbrev{k} rule],'-dpdf')
end

%% Testing equality of capital and labor exponents

[~,countryHoga] = xlsread(filename,sheetHoga,'A4:A999');
resultsHoga = xlsread(filename,sheetHoga,'B4:I999');

% scatter plot of two labor Pareto exponents
alphaL_Hoga = resultsHoga(:,3);
ind = find(resultsHoga(:,end) >= Nmin);
NminLarge = 2000;
ind1 = find(resultsHoga(:,end) >= NminLarge);
ind2 = find((resultsHoga(:,end) >= Nmin)&(resultsHoga(:,end) < NminLarge));

[R,P,RL,RU] = corrcoef(alphaL(ind),alphaL_Hoga(ind));
rho = R(1,2);

aU = ceil(max([alphaL(ind); alphaL_Hoga(ind)]));

figure
scatter(alphaL(ind),alphaL_Hoga(ind)); hold on
plot([0 aU],[0 aU],'k-','LineWidth',1);
text(0.1,aU-2,['$\rho= $' num2str(round(rho,2,'significant'))],'VerticalAlignment','bottom');
axis equal
xlabel([rule '\% rule for full sample'])
ylabel([rule '\% rule for positive capital income'])
legend('Pareto exponent','$45^\circ$ line','Location','NW')
xlim([0,aU])
ylim([0,aU])

%{
figure
plot(alphaL(ind1),alphaL_Hoga(ind1),'o','Color',c1); hold on
plot(alphaL(ind2),alphaL_Hoga(ind2),'x','Color',c2);
plot([0 aU],[0 aU],'k-','LineWidth',1);
text(0.1,aU-2,['$\rho= $' num2str(round(rho,2,'significant'))],'VerticalAlignment','bottom');
axis equal
xlabel([rule '\% rule for full sample'])
ylabel([rule '\% rule for positive capital income'])
legend('$N\ge 5000$','$N< 5000$','$45^\circ$ line','Location','NW')
xlim([0,aU])
ylim([0,aU])
%}

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_alphaLScatter' rule],'-dpdf')

countryHoga = countryHoga(ind);
resultsHoga = resultsHoga(ind,:);

Ncap = size(resultsHoga,1)
crit = 55.44; % 5% critical value for Hoga test

TN = resultsHoga(:,2); % Hoga statistic
ind_reject = find(TN >= crit);
Nreject = length(ind_reject)
rejectPercent = round(100*Nreject/Ncap)

% count number of samples alphaL exceeds alphaK
sum((TN >= crit)&(resultsHoga(:,3)>=resultsHoga(:,4)))

figure
histogram(log(TN),'normalization','pdf'); hold on
xline(log(crit),'--');
text(log(crit)-0.1,0.3,'5\% critical value','VerticalAlignment','top','HorizontalAlignment','right');
xlabel('$\log T_N$')
ylabel('Probability density')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_TNHist' rule],'-dpdf')

%% appendix tables

% capital and labor income Pareto exponents

results = xlsread(filename,sheet,'B4:J999');

ind = find(results(:,end) < Nmin);
results(ind,end-3:end) = NaN;

t_alpha = table(country,results(:,1),round(results(:,2),2),round(results(:,3),2),...
    results(:,4),results(:,5),round(results(:,6),2),round(results(:,7),2),...
    results(:,8),results(:,9));

table2latex(t_alpha,'Tab_alpha');

% Hoga test

YesNo = cell(Ncap,1);
YesNo(:) = {'No'};
YesNo(ind_reject) = {'Yes'};

t_Hoga = table(countryHoga,resultsHoga(:,1),round(resultsHoga(:,3),2),...
    round(resultsHoga(:,5),2),round(resultsHoga(:,4),2),...
    round(resultsHoga(:,6),2),resultsHoga(:,7),resultsHoga(:,8),YesNo);

table2latex(t_Hoga,'Tab_Hoga');
