function [grid,s] = expGrid(a,b,c,N)
% construct N-point exponential grid on (a,b] with median grid point c

% some error checking
if (a >= b)||(c <= a)||(c >= (a+b)/2)
    error('it must be a<c<(a+b)/2')
end

s = (c^2-a*b)/(a+b-2*c); % shift parameter
temp = linspace(log(a+s),log(b+s),N+1); % even grid in log scale
grid = exp(temp(2:end))-s;

end

