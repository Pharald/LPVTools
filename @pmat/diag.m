function b = diag(a,k)
% DIAG  Diagonal PMATs and diagonals of a PMAT object.
%
% If A is a PMAT vector of length N then B=diag(A) is an N-by-N PMAT
% with the entries of A on the main diagonal of B. If A is an N-by-M PMAT
% (N,M > 1) then B=diag(A) is a column vector whose elements are the 
% entries on the main diagonal of A.
%
% If A is a PMAT vector of length N then B=diag(A,K) is PMAT matrix with
% the elements of A on the K^th diagonal of B. 
%
% If A is an N-by-M PMAT (N,M > 1) then B=diag(A,K) is a column vector
% whose elements are the entries of the K^th diagonal of A.
%
% See also: diag, triu, tril, blkdiag.

% Input / output error checking
nin = nargin;
narginchk(1, 2)

if nin==1
    k=0;
elseif isa(k,'pmat')
    % TODO PJS 4/3/2011: Revisit handling of PMAT inputs for k.
    % Current PMATs must have constant dimensions at each point in the
    % domain. Varying k would produce output that violates this.
    error('K-th diagonal input must be an integer scalar.');
    
%     if ~isscalar(k)
%         error('K-th diagonal input must be an integer scalar.');
%     end    
%     k = k.DataPrivate;
%     k1 = k(1);
%     if all( k(:) == k1 )
%         k=k1;
%     else
%         error(['For diag(a,k), k must have constant value at' ...
%             ' each point in the domain.']);
%     end
end

sza = privatesize(a);
sza = sza(1:2);
if sza(1)==1 || sza(2)==1
    % a is a vector --> b is a matrix
    
    % Use double/diag to find right dimensions and indices
    % TODO GJB 25Jul12 - check PMAT LENGTH is being used correctly
    la = length(a);
    bmat = diag(1:la,k);
    idx = find(bmat);
    
    % Assign values
    b = pmat( zeros(size(bmat)) );
    L.type = '()';
    L.subs = {idx};
    b = subsasgn(b,L,a);
else
    % a is a matrix --> b is a vector
    if k>=sza(2) || k <= -sza(1)
        b = pmat(zeros(0,1));
        return
    end
    
    % Use double/diag to find right indices
    amat = reshape( 1:(sza(1)*sza(2)), sza );
    idx = diag(amat,k);
    
    % Assign values
    L.type = '()';
    L.subs = {idx};
    b = subsref(a,L);
end



