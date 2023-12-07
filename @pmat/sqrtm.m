function [varargout] = sqrtm(mat)
% SQRTM   Matrix square root for PMAT objects.
%
% SQRTM(M) is the matrix square root of M at each point in the domain of M.
%
% See also: sqrtm, expm, logm.

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 1)
nargoutchk(0, 3)

out1 = mat;
szm = privatesize(mat);
npts = prod(szm(3:end));
if nout==1
    for i=1:npts
        out1.DataPrivate(:,:,i) = sqrtm(mat.DataPrivate(:,:,i));
    end 
    varargout = {out1};
elseif nout==2
    out2 = zeros(1,1,npts);
    for i=1:npts
        [out1.DataPrivate(:,:,i),out2(:,:,i)] = sqrtm(mat.DataPrivate(:,:,i));
    end 
    out2 = pmat(out2,out1.DomainPrivate);
    varargout = {out1,out2};
else
    out2 = zeros(1,1,npts);
    out3 = zeros(1,1,npts);
    for i=1:npts
        [out1.DataPrivate(:,:,i),out2(:,:,i),out3(:,:,i)] = ...
            sqrtm(mat.DataPrivate(:,:,i));
    end 
    out2 = pmat(out2,out1.DomainPrivate);
    out3 = pmat(out3,out1.DomainPrivate);
    varargout = {out1,out2,out3};
end




% ----- OLD CODE

% out = mat;
% szm = size(mat);
% if szm(1)==szm(2)
%     for i=1:prod(szm(3:end))
%         out.DataPrivate(:,:,i) = sqrtm(mat.DataPrivate(:,:,i));
%     end
% else
%     error('Matrix must be square');
% end

