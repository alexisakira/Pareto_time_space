%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getTopShares
% (c) 2019 Emilien Gouin-Bonenfant and Alexis Akira Toda
% 
% Purpose: 
%       Compute top wealth shares from wealth grid and stationary
%       distribution
%
% Usage:
%       topShare = getTopShares(topProb,wGrid,wDist,zeta)
%
% Inputs:
% topProb   - top probabilities to evaluate top shares. For example,
%           topProb = [0.001 0.01 0.1] evaluates top 0.1%, 1%, 10% shares
% wGrid     - wealth grid
% wDist     - wealth distribution on wGrid
% zeta      - Pareto exponent (set 0 if using truncation)
%
% Outputs:
% topShare  - top wealth share corresponding to topProb
%
% Version 1.1: June 16, 2019
%
% Version 1.2: April 20, 2020
% - Fixed bug in spline interpolation when the grid is fine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function topShare = getTopShares(topProb,wGrid,wDist,zeta)
%% some error checking

if nargin < 4
    zeta = 0; % use truncation
end

if (min(topProb) < 0)||(max(topProb) > 1)
    error('topProb must be vector of numbers between 0 and 1')
end

if any(diff(topProb) <= 0)
    error('topProb must be strictly increasing')
end

if any(diff(wGrid) <= 0)
    error('wGrid must be strictly increasing')
end

if length(wGrid) ~= length(wDist)
    error('length of wGrid and wDist must agree')
end

if size(wGrid,1) < size(wGrid,2)
    wGrid = wGrid'; % convert to column vector
end

if size(wDist,1) < size(wDist,2)
    wDist = wDist'; % convert to column vector
end

if size(topProb,1) < size(topProb,2)
    topProb = topProb'; % convert to column vector
end

tailProb = cumsum(flipud(wDist)); % tail probability
[~,ia,~] = unique(tailProb); % index of unique values

%% first, consider when using truncation
if zeta == 0
    aggW = dot(wDist,wGrid); % aggregate wealth
    topWealth = cumsum(flipud(wDist.*wGrid))/aggW; % top wealth shares on grid
    topShare = interp1([0;tailProb(ia)],[0;topWealth(ia)],topProb,'spline'); % spline interpolation
    return
end

%% next, consider Pareto extrapolation
if zeta <= 1
    error('zeta must exceed 1') % need zeta > 1 for finite mean
end

temp = wGrid;
temp(end) = (zeta/(zeta-1))*wGrid(end); % correct last grid point
aggW = dot(wDist,temp); % aggregate wealth using Pareto extrapolation
topWealth = cumsum(flipud(wDist.*temp))/aggW; % top wealth shares on grid

ind1 = find(topProb <= tailProb(1)); % index for which extrapolation is necessary
ind2 = find(topProb > tailProb(1)); % index for which extrapolation is unnecessary

topShare = 0*topProb;
topShare(ind1) = (zeta/(zeta-1))*wDist(end)^(1/zeta)*(wGrid(end)/aggW)*topProb(ind1).^(1-1/zeta);
% extrapolate top wealth shares using Pareto distribution
topShare(ind2) = interp1(tailProb(ia),topWealth(ia),topProb(ind2),'spline'); % spline interpolation

end