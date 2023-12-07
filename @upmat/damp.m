function [Wn,Z,P] = damp(S)
% DAMP  Pointwise compute natural frequency and damping of UPMAT objects.
% 
% [Wn,Z] = DAMP(M) returns vector PMATs Wn and Z containing the natural 
% frequencies and damping factors of the UPMAT M at each point in the 
% domain of M.
%
% [Wn,Z,P] = DAMP(M) also returns the poles P of M.
%
% See also: damp, pole, zero.

szS = privatesize(S);
SData = S.DataPrivate;


% Compute nat freq and damping at each point in domain
npts = prod(szS(3:end));
if nout==2
    % Initialize outputs
    [Wn1,Z1]=damp(SData(:,:,1));
    Wn = zeros( [size(Wn1) szS(3:end)] );
    Wn(:,:,1) = Wn1;
    Z = zeros( [size(Z1) szS(3:end)] );
    Z(:,:,1) = Z1;

    for i=2:npts
       [Wn(:,:,i),Z(:,:,i)] = damp( SData(:,:,i) );
    end
    Wn = pmat(Wn,S.DomainPrivate);
    Z = pmat(Z,S.DomainPrivate);
else
    % Initialize outputs
    [Wn1,Z1,P1]=damp(SData(:,:,1));
    Wn = zeros( [size(Wn1) szS(3:end)] );
    Wn(:,:,1) = Wn1;
    Z = zeros( [size(Z1) szS(3:end)] );
    Z(:,:,1) = Z1;
    P = zeros( [size(P1) szS(3:end)] );
    P(:,:,1) = P1;
    
    for i=2:npts
        [Wn(:,:,i),Z(:,:,i),P(:,:,i)] = damp( SData(:,:,i) );
    end
    Wn = pmat(Wn,S.DomainPrivate);
    Z = pmat(Z,S.DomainPrivate);
    P = pmat(P,S.DomainPrivate);
end
