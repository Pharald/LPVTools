function sys = repsys(sys,m,n)
% REPSYS   Replicate and tile for UPFRD objects
%   
% B = REPSYS(A,M,N) creates a larger UPFRD B consisting of an M-by-N tiling 
% of copies of A. B = REPSYS(A,[M N]) also creates an M-by-N tiling of A.
%
% See also: repsys.

% Argument checking
if nargin==2
    if isscalar(m)
        n = m;
    elseif isvector(m) && length(m)==2
        n = m(2);
        m = m(1); %JT was n=m(1)
    else
        error('Only the row and column dimensions can be replicated.');
    end
end
if isa(m,'upmat') || isa(n,'upmat');
    error('Tiling dimensions M and N must be non-negative integers.');
end

% Replicate SS Array
niv = sys.DomainPrivate.NumIV;
sys.DataPrivate = repsys(sys.DataPrivate,[m n ones(1,niv)]);
