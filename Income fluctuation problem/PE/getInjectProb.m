%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getInjectProb
% (c) 2020 Emilien Gouin-Bonenfant and Alexis Akira Toda
% 
% Purpose: 
%       Compute injection probabilities from outside the grid using Pareto extrapolation
%
% Usage:
%       q = getInjectProb(xGrid,g,G,zeta)
%
% Inputs:
% xGrid - (1 x N) grid for state variable x (asset, wealth, etc.)
% g     - law of motion at largest grid point
% G     - slope of law of motion at largest grid point
%
% Optional:
% zeta  - Pareto exponent
%
% Output:
% q     - (1 x N) transition probabilities from largest grid point
%
% Version 1.0: February 4, 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function q = getInjectProb(xGrid,g,G,zeta)
%% some error checking

N = length(xGrid); % number of wealth grid points

if N < 2
    error('xGrid must contain at least two elements')
end

if any(diff(xGrid) <= 0)
    error('xGrid must be in strictly increasing order')
end

if xGrid(end-1) <= 0
    error('second largest grid point must be positive')
end

if nargin < 4 % Pareto exponent not provided
    zeta = Inf; % just use truncation method
end

if zeta <= 1
    error('zeta must exceed 1')
end

%% compute injection probabilities from outside the grid
a = xGrid(1:end-1); % lower end of integration
b = xGrid(2:end); % upper end of integration
alpha = max(0,a-g);
beta = max(0,b-g);

if (zeta == Inf)||(G <= 0)
    aveF = (beta-alpha)./(b-a); % average CDF between grid points
else % zeta > 1 and G > 0
    C = G*xGrid(end);
    %aveF = (beta-alpha -(C/(1-zeta))*((1+beta/C).^(1-zeta) - (1+alpha/C).^(1-zeta)))./(b-a);
    temp = xGrid(end-1)/xGrid(end);
    theta = (1-temp)*(zeta-1)/(temp^(1-zeta)-1); % correction term for quadrature
    %theta = (temp^zeta + theta)/2; % this one works well
    %theta = 1/(1+(zeta/2)*(1-temp)); % trapeziodal formula
    
    %wb = (1/(1/temp-1))*(1/(zeta-1)-zeta/(zeta-1)*temp^(zeta-1)+temp^zeta);
    aveF = (beta-alpha -theta*(C/(1-zeta))*((1+beta/C).^(1-zeta) - (1+alpha/C).^(1-zeta)))./(b-a);
end

q = [aveF(1) diff(aveF) 1-aveF(end)]; % injection probabilities

end

