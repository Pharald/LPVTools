function mat = repmat(mat,m,n)
% REPMAT   Replicate and tile for PMAT objects
%   
% B = REPMAT(A,M,N) creates a larger PMAT B consisting of an M-by-N tiling 
% of copies of A. B = REPMAT(A,[M N]) also creates an M-by-N tiling of A.
%
% B = REPMAT(A,M) creates a B that is an M-by-M tiling of A.
%
% See also: repmat.

% Argument checking
% TODO PJS 4/3/2011: Revisit handling of PMATs for m/n and error messages.
% NOTE PJS 4/29/2011: If m and/or n is a PMAT then the dimensions of
%  the output mat will, in general, depend on the point in the domain.
%  This can't be handled within the current PMAT class and so varying
%  m and n should be disallowed.

if nargin==2
   % TODO 7/26/2012:
   % Check mat is a PMAT and m is a vector of integers
elseif nargin==3
   % TODO 7/26/2012:
   % Check mat is a PMAT and m/n are scalar integers
   m=[m n];   
end

niv = mat.Domain.NumIV;
Data = mat.Data;
m = m(:)';
Data = repmat(Data,[m(1:2) ones(1,niv) m(3:end)]);
mat = pmat(Data,mat.Domain);


