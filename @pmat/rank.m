function out = rank(mat,tol)
% RANK   Matrix rank for PMAT objects.
%
% RANK(M) is an estimate of the matrix rank of M at each point in the
% domain of M. The rank is based on the number of singular values
% that exceed a default tolerance.  The RANK help for double matrices
% provides more details.
%
% RANK(M,TOL) uses the tolerance TOL instead of the default tolerance.
%
% See also: rank.


% Input / output error checking
nin = nargin;
narginchk(1, 2)

if nin==1
    szm = privatesize(mat);
    niv = mat.DomainPrivate.NumIV;
    Data = mat.DataPrivate;
    
    out = zeros([1 1 szm(3:2+niv)]);
    for i=1:prod(szm(3:end))
        out(:,:,i) = rank(Data(:,:,i));
    end
else
    [mat,tol] = domunion(mat,tol);
    szm = privatesize(mat);
    niv = mat.DomainPrivate.NumIV;
    matData = mat.DataPrivate;
    tolData = tol.DataPrivate;

    out = zeros([1 1 szm(3:2+niv)]);
    for i=1:prod(szm(3:end))
        out(:,:,i) = rank(matData(:,:,i),tolData(:,:,i));
    end
end
out = pmat(out,mat.DomainPrivate);

