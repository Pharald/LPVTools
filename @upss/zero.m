function [Z,GAIN] = zero(S)
% ZERO  Pointwise computes the transmission zeros of a UPSS object
%  
% Z = ZERO(SYS) returns the transmission zeros of the UPSS SYS at each
% point in the domain of SYS. Z is returned as a column vector PMAT 
% containing the transmission zeros at each point in the domain of SYS.
%
% [Z,GAIN] = ZERO(SYS) also returns the transfer function gain for
% SISO models SYS.
%
% See also: zero, damp, pole.

nout = nargout;

szS = privatesize(S);
SData = S.DataPrivate;
nx = size(SData(:,:,1).a,1);

Z = nan([nx 1 szS(3:end)]);
nzmax = 0;
npts = prod(szS(3:end));
if nout==1
    for i=1:npts
       Zi = zero( SData(:,:,i) );
       nz = length( Zi );
       Z(1:nz,:,i) = Zi;
       nzmax = max(nz,nzmax);       
    end
else
    GAIN = zeros([1 1 szS(3:end)]);
    for i=1:npts
       [Zi,GAIN(:,:,i)] = zero( SData(:,:,i) );
       nz = length( Zi );
       Z(1:nz,:,i) = Zi;
       nzmax = max(nz,nzmax);       
    end
    GAIN = pmat(GAIN,S.DomainPrivate);    
end
Z = reshape(Z(1:nzmax,:),[nzmax 1 szS(3:end)]);
Z = pmat(Z,S.DomainPrivate);

