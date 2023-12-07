function C = norm(M,P)
% NORM   Matrix or vector norm of a PMAT object.
%
% NORM(M) is the norm of M at each point in the domain of M. 
%
% NORM(M,P) specifies the P-norm. For vectors, P can be any scalar double.
% For matrices, P can only be 1, 2, 'inf' or 'fro'.
%
% See also: norm, cond, rcond.

% TODO PJS 4/30/2011: Overload this function to allow the P to be a
% varying PMAT?  The case P='fro' cannot be handled with PMATs and
% would need to be handled separately.

% Check # of input/output arguments
nin = nargin;
narginchk(1, 2)

if nargin == 1
   P = 2;
end
% Don't allow P to be PMAT (for time being, see comment above.)
if ~isnumeric(P) && ~ischar(P)
    error('Input P can only be scalar double, ''inf'' or ''fro''')
end

szm = privatesize(M);
C = zeros([1 1 szm(3:end)]);
Data = M.DataPrivate;
for i=1:prod(szm(3:end))
   C(1,1,i) = norm(Data(:,:,i),P);
end
C = pmat(C,M.DomainPrivate);    

    
