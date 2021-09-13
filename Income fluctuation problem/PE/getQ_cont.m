%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getQ
% (c) 2019 Emilien Gouin-Bonenfant and Alexis Akira Toda
% 
% Purpose: 
%       Compute transition probability matrix using Pareto extrapolation
%
% Usage:
%       [Q,pi] = getQ(PS,PJ,p,xGrid,x0,gstjn,Gstj,zeta)
%
% Inputs:
% PS    - (S x S) transition probability matrix of exogenous state
% PJ    - (S^2 x J) matrix of conditional probabilities of transitory state
%       if (1 x J), then assume distribution of j does not depend on (s,s')
%       if (S x J), then assume distribution of j depends only on s
% p     - birth/death probability (set 0 for infinitely-lived case)
% x0    - initial state variable of newborn agents (not used if p=0)
% xGrid - (1 x N) grid for state variable x (asset, wealth, etc.)
% gstjn - (S^2 x NJ) matrix of law of motion of x
%       if (S x NJ), then assume x does not depend on s'
%
% Optional:
% Gstj  - (S^2 x J) matrix of asymptotic slope of law of motion of x
%       if (S x J), then assume G does not depend on s'
% zeta  - Pareto exponent
%
% Output:
% Q     - (SN x SN) transition probability matrix of (s,x)
%
% Optinal:
% pi    - (SN x 1) stationary distribution of Q
%
% Version 1.2: February 4, 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [Q,pi] = getQ_cont(PS,PJ,p,x0,xGrid,gstjn,Gstj,zeta)
%% some error checking

S = size(PS,2); % number of exogenous states
J = size(PJ,2); % number of transitory states
N = length(xGrid); % number of wealth grid points

if N < 2
    error('xGrid must contain at least two elements')
end

if size(PS,1) ~= S
    error('PS must be a square matrix')
end

if size(PJ,1) == 1 % conditional distribution independent of states
    PJ = repmat(PJ,S^2,1);
elseif size(PJ,1) == S % conditional distribution depends only current state
    PJ = kron(PJ,ones(S,1));
elseif size(PJ,1) ~= S^2
    error('size of PS and PJ inconsistent')
end

if (p < 0)||(p >= 1)
    error('p must be in [0,1)')
end

if any(diff(xGrid) <= 0)
    error('xGrid must be in strictly increasing order')
end

xMax = max(xGrid);
if xMax <= 0
    error('largest grid point must be positive')
end

if size(gstjn,2) ~= N*J
    error('size of PJ, xGrid, and gstjn must be consistent')
end

if p>0 % there is birth/death
    if x0 < xGrid(1)
        error('it must be x0 >= min(xGrid)')
    elseif x0 >= xGrid(end)
        error('it must be x0 < max(xGrid)')
    end   
end

if size(gstjn,1) == S % law of motion does not depend on next state
    gstjn = kron(gstjn,ones(S,1)); % replicate rows to make it S^2 x NJ
end
if size(gstjn,1) ~= S^2
    error('size of PS and gstjn must be consistent')
end

%% define optional arguments if not provided

if (nargin < 7)||isempty(Gstj) % asymptotic slope not provided
    ind = N*[1:J]; % index of law of motion at largest grid point
    Gstj = (gstjn(:,ind) - gstjn(:,ind-1))/(xGrid(end) - xGrid(end-1));
    % compute slope from two largest grid points
end

if size(Gstj,1) == S % law of motion does not depend on next state
    Gstj = kron(Gstj,ones(S,1)); % replicate rows to make it S^2 x NJ
end
if size(Gstj,1) ~= S^2
    error('size of PS and Gstj must be consistent')
end

if any(Gstj(:) <= 0)
    error('asymptotic slope must be positive')
end

if nargin < 8 % Pareto exponent not provided
    zeta = getZeta(PS,PJ,p,Gstj);
end

if zeta <= 1
    warning('zeta must be larger than 1 for finite mean')
    Q = NaN(S*N);
    return
end

%% construct the combined transition probability matrix

Q = zeros(S*N); % transition probability matrix of (s,x)

for s = 1:S
    for t = 1:S
        pst = PS(s,t); % P(t | s)
        for j = 1:J
            pstj = pst*PJ(S*(s-1)+t,j); % P(j,t | s)
            % use nonstochastic simulation except largest grid point
            for n = 1:N-1 % actual grid points
                x = gstjn(S*(s-1)+t,N*(j-1)+n); % destination of x
                ind = find(xGrid < x, 1, 'last' );
                if isempty(ind) % x is below smallest grid point
                    Q(N*(s-1)+n,N*(t-1)+1) = Q(N*(s-1)+n,N*(t-1)+1) + pstj;
                elseif ind < N % x is strictly inside grid, use convex combination
                    theta = (x-xGrid(ind))/(xGrid(ind+1)-xGrid(ind));
                    Q(N*(s-1)+n,N*(t-1)+ind) = Q(N*(s-1)+n,N*(t-1)+ind) + (1-theta)*pstj;
                    Q(N*(s-1)+n,N*(t-1)+ind+1) = Q(N*(s-1)+n,N*(t-1)+ind+1) + theta*pstj;
                else % x is above largest grid point
                    Q(N*(s-1)+n,N*t) = Q(N*(s-1)+n,N*t) + pstj;
                end
            end
            % Pareto extrapolation
            
            g = gstjn(S*(s-1)+t,N*j); % intercept
            %G = (gstjn(S*(s-1)+t,N*j) - gstjn(S*(s-1)+t,N*j-1))/(xGrid(end) - xGrid(end-1)); % slope
            G = Gstj(S*(s-1)+t,j);
            Q(N*s,N*(t-1)+1:N*t) = Q(N*s,N*(t-1)+1:N*t) + pstj*getInjectProb(xGrid,g,G,zeta);
        end
    end
end

%% finally, adjust for birth/death
if p > 0
    Qx0 = zeros(N);
    ind0 = find(xGrid < x0, 1, 'last');
    if isempty(ind0) % x0 is below smallest grid point
        Qx0(:,1) = 1;
    elseif ind0 < N % x0 is strictly inside grid
        theta0 = (x0-xGrid(ind0))/(xGrid(ind0+1)-xGrid(ind0));
        Qx0(:,ind0) = 1-theta0; % convex combination
        Qx0(:,ind0+1) = theta0;
    else % x0 is above largest grid
        Qx0(:,end) = 1;
    end
    Q0 = kron(PS,Qx0);
    Q = (1-p)*Q + p*Q0;
end

if nargout > 1
    warning off
    [v,~] = eigs(Q',1,1);
    pi = v/sum(v);
end

end

