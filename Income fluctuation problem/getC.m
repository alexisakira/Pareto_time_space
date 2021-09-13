function ca = getC(beta,gamma,P,R,aGrid,MaxIter,tol)
% compute consumption function by policy function iteration

%% some error checking
if gamma <= 0
    error('gamma must be positive')
end
if any(P < 0)
    error('P must be nonnegative')
end
S = length(P);
if (length(beta) ~= 1)&&(length(beta) ~= S)
    error('length of beta must be 1 or same as P')
end
if any(beta <= 0)
    error('beta must be positive')
end
if length(R) ~= S
    error('length of R must be same as P')
end
if any(R <= 0)
    error('R must be positive')
end
if size(P,2) > size(P,1)
    P = P'; % convert to column vector
end
if size(R,2) > size(R,1)
    R = R'; % convert to column vector
end
if any(diff(aGrid <= 0))
    error('aGrid must be strictly increasing')
end
if aGrid(1) <= 0
    error('elements of aGrid must be positive')
end
if size(aGrid,1) > size(aGrid,2)
    aGrid = aGrid'; % convert to row vector
end
N = length(aGrid);
if nargin < 6
    MaxIter = 400;
end
if nargin < 7
    tol = 1e-6;
end

%% policy function iteration

MPC = 1 - (P'*(beta.*R.^(1-gamma)))^(1/gamma); % asymptotic MPC;

ca = min([aGrid; 1 + MPC*aGrid],[],1); % initialize consumption function
for i = 1:MaxIter
    aNextMat = R*(aGrid - ca) + 1; % S x N matrix of next period's asset
    cNextMat = zeros(S,N); % S x N matrix of next period's consumption
    for s = 1:S
        temp = interp1([0 aGrid],[0 ca],aNextMat(s,:),'linear','extrap'); % linear interpolation
        cNextMat(s,:) = max(temp,1e-10); % force it to be nonnegative
    end
    expect = P'*bsxfun(@times,beta.*R,cNextMat.^(-gamma));
    ca_new = min([expect.^(-1/gamma); aGrid],[],1); % update using Euler equation
    ca_new = cummax(ca_new); % force it to be an increasing function
    if max(abs(ca_new./ca - 1)) < tol
        disp('Converged!')
        break
    else
        ca = ca_new;
    end
end

end

