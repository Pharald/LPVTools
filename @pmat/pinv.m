function out = pinv(mat,tol)
% PINV   Pseudoinverse of a PMAT object.
%
% B=PINV(M) is the pseudoinverse of M at each point in the domain of M.
% Singular values less than the default tolerance are treated as zero.
% The PINV help for double matrices provides more details.
%
% PINV(M,TOL) uses the tolerance TOL instead of the default tolerance.
%
% See also: pinv, rank.

% Input / output error checking
nin = nargin;
narginchk(1, 2)

if nin==1
    szm = privatesize(mat);
    niv = mat.DomainPrivate.NumIV;
    Data = mat.DataPrivate;
    
    out = zeros([szm(2) szm(1) szm(3:2+niv)]);
    for i=1:prod(szm(3:end))
        out(:,:,i) = pinv(Data(:,:,i));
    end
else
    [mat,tol] = domunion(mat,tol);
    szm = privatesize(mat);
    niv = mat.DomainPrivate.NumIV;
    matData = mat.DataPrivate;
    tolData = tol.DataPrivate;

    out = zeros([szm(2) szm(1) szm(3:2+niv)]);
    for i=1:prod(szm(3:end))
        out(:,:,i) = pinv(matData(:,:,i),tolData(:,:,i));
    end
end
out = pmat(out,mat.DomainPrivate);

