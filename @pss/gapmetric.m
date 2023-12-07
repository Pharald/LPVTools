function [GAP,NUGAP] = gapmetric(P1,P2,tol)
% GAPMETRIC Pointwise gap and the Vinnicombe gap metric for PSS objects
% 
% [GAP,NUGAP] = gapmetric(P1,P2) calculates the gap metric between the
% systems P1 and P2 at each point in the parameter domain. See lti/gapmetric
% for details.
% 
% See also: norm, loopmargin, wcmargin.

% Check number of input arguments
% nargchk(2,3);

% Run gapmetric at each point in domain
% (Gapmetric does not work on SS arrays)
[P1ext,P2ext] = domunion(P1,P2);
P1Data = P1ext.DataPrivate;
P2Data = P2ext.DataPrivate;

Dom = P1ext.DomainPrivate;
szd = size(Dom);
numd = prod(szd);

GAP = zeros(numd,1);
NUGAP = zeros(numd,1);
for i=1:numd
    if nargin==2
        [GAP(i),NUGAP(i)] = gapmetric(P1Data(:,:,i),P2Data(:,:,i));
    else
        [GAP(i),NUGAP(i)] = gapmetric(P1Data(:,:,i),P2Data(:,:,i),tol);
    end
end
GAP = reshape(GAP,szd);
NUGAP = reshape(NUGAP,szd);

GAP = pmat(shiftdim(GAP,-2),Dom);
NUGAP = pmat(shiftdim(NUGAP,-2),Dom);
