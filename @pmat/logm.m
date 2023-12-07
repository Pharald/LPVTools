function [varargout] = logm(mat)
% LOGM   Matrix logarithm for PMAT objects.
%
% LOGM(M) is the matrix logarithm of M at each point in the domain of M.
%
% See also: logm, expm, sqrtm.

% TODO PJS 4/1/2011: Implement funm?

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 1)
nargoutchk(0, 2)

out1 = mat;
szm = privatesize(mat);
npts = prod(szm(3:end));
if nout==1
    for i=1:npts
        out1.DataPrivate(:,:,i) = logm(mat.DataPrivate(:,:,i));
    end 
    varargout = {out1};
else
    out2 = zeros(1,1,npts);
    for i=1:npts
        [out1.DataPrivate(:,:,i),out2(:,:,i)] = logm(mat.DataPrivate(:,:,i));
    end 
    out2 = pmat(out2,out1.DomainPrivate);
    varargout = {out1,out2};
end





