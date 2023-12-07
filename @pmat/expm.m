function out = expm(mat)
% EXPM   Matrix exponential for PMAT objects.
%
% B = EXPM(M) is the matrix exponential of M at each point in the domain
% of M.
%
% See also: expm, logm, sqrtm.

out = mat;
szm = privatesize(mat);
for i=1:prod(szm(3:end))
    out.DataPrivate(:,:,i) = expm(mat.DataPrivate(:,:,i));
end
