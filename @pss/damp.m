function [Wn,Z,P] = damp(S)
% DAMP  Pointwise compute natural frequency and damping of PSS objects.
% 
% [Wn,Z] = DAMP(SYS) returns vector PMATs Wn and Z containing the natural 
% frequencies and damping factors of the PSS SYS at each point in the 
% domain of SYS. For discrete-time models, the equivalent s-plane natural 
% frequency and damping ratio of an eigenvalue lambda are:
%     Wn = abs(log(lambda))/Ts ,   Z = -cos(angle(log(lambda))) 
% Wn and Z are empty vectors if the sample time Ts is undefined.
%
% [Wn,Z,P] = DAMP(SYS) also returns the poles P of SYS.
%
% See also: damp, pole, zero.

P = pole(S);

nout = nargout;
szS = privatesize(S);
SData = S.DataPrivate;

if nout==2
    % Initialize outputs
    [Wn1,Z1]=damp(SData(:,:,1));
    Wn = zeros( [size(Wn1) szS(3:end)] );
    Wn(:,:,1) = Wn1;
    Z = zeros( [size(Z1) szS(3:end)] );
    Z(:,:,1) = Z1;

    % Compute nat freq and damping at each point in domain
    npts = prod(szS(3:end));
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

    % Compute nat freq and damping at each point in domain
    npts = prod(szS(3:end));
    for i=2:npts
       [Wn(:,:,i),Z(:,:,i),P(:,:,i)] = damp( SData(:,:,i) );
    end
    Wn = pmat(Wn,S.DomainPrivate);
    Z = pmat(Z,S.DomainPrivate);
    P = pmat(P,S.DomainPrivate);
end
