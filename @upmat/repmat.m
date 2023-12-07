function mat = repmat(mat,m,n)
% REPMAT   Replicate and tile for UPMAT objects
%   
% B = REPMAT(A,M,N) creates a larger UPMAT B consisting of an M-by-N tiling 
% of copies of A. B = REPMAT(A,[M N]) also creates an M-by-N tiling of A.
%
% See also: repmat.

% Argument checking
% TODO PJS 4/3/2011: Revisit handling of PMATs for m/n and error messages.
% NOTE PJS 4/29/2011: If m and/or n is a PMAT then the dimensions of
%  the output mat will, in general, depend on the point in the domain.
%  This can't be handled within the current PMAT class and so varying
%  m and n should be disallowed.
if nargin==2
    if isscalar(m)
        n = m;
    elseif isvector(m) && length(m)==2
        n = m(2);
        m = m(1);
    else
        error('Only the row and column dimensions can be replicated.');
    end
end
if isa(m,'upmat') || isa(n,'upmat');
    error('Tiling dimensions M and N must be non-negative integers.');
end

% Replicate matrix data
niv = mat.DomainPrivate.NumIV;
mat.DataPrivate = repmat(mat.DataPrivate,[m n ones(1,niv)]);
