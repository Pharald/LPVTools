function S = IQCtvreal(poleList,RateBound)
%  IQCtvreal  Creates IQC for rate-bounded, time-varying, real parameter
%     using swapping lemma
%  
%     S = IQCtvreal(POLELIST,RATEBOUND) creates IQC description struct for
%     a rate-bounded, real parameter DELTA, which satsifies |DELTA(t)|<=1,
%     and |(d/dt) DELTA(t)| <= RATEBOUND.  RATEBOUND should be a positive
%     scalar.   POLELIST is a vector of negative real numbers, used to 
%     parametrize first-order basis functions for the allowable dynamic IQC
%     multipliers.

S.IQCfunction = @LOCALIQCtvreal;
S.IQCparams.polelist = poleList;
S.IQCparams.RB = RateBound;
S.PsiFlag = true;

function [lmiout,IDcenter,Psi] = LOCALIQCtvreal(lmisys,ncopies,omeg,params)
% Assume tvreal is unit norm-bounded with symmetric rate bound RB
% polelist is a vector of real, negative poles used to construct Psi

polelist = params.polelist;
RB = params.RB;

% Grab LMI info
setlmis(lmisys)
ndec = decnbr(lmisys);
nlmis = lminbr(lmisys);

% Form state matrices for filters
% XXX We can reconsider the choice for A and B
A = diag(polelist);
nst = size(A,1);
B = ones(nst,ncopies); 
C = eye(nst);
D = zeros(nst,ncopies);

% Define matrix variables and filter for the norm bound constraint
[X1,ndec,X1dec] = lmivar( 3, symdec(nst+ncopies,ndec) );
[Y1,ndec,Y1dec] = lmivar( 3, skewdec(nst+ncopies,ndec) );
[M1,ndec,M1dec] = lmivar( 3, [X1dec Y1dec; Y1dec' -X1dec] );
% XXXX: Reversed order of lines below, since before these execute, NLMIS is the
% number of LMIs already in system.  We're trying to add a new one.
nlmis = nlmis+1;
lmiterm([-nlmis 1 1 X1],1,1);

Gfr = frd(ss(A,B,C,D),omeg);
sIAfr = frd(ss(A,eye(nst),C,0),omeg);
Psi1 = [Gfr zeros(nst,ncopies+nst);...
    eye(ncopies) zeros(ncopies,ncopies+nst); ...
    zeros(nst,ncopies) Gfr sIAfr; ...
    zeros(ncopies) eye(ncopies) zeros(ncopies,nst) ];

% Define matrix variables and filter for the rate bound constraint
[X2,ndec,X2dec] = lmivar( 3, symdec(nst,ndec) );
if nst>1
    [Y2,ndec,Y2dec] = lmivar( 3, skewdec(nst,ndec) );
    [M2,ndec,M2dec] = lmivar( 3, [X2dec Y2dec; Y2dec' -X2dec] );
else
    [M2,ndec,M2dec] = lmivar( 3, [X2dec 0; 0 -X2dec] );
end
% XXXX: Similar to above, although immediate line below was missing
nlmis = nlmis+1;
lmiterm([-nlmis 1 1 X2],1,1);

Psi2 = [ RB*Gfr zeros(nst,ncopies+nst); zeros(nst,2*ncopies) eye(nst) ];

% Pack up outputs
lmiout = getlmis;
IDcenter = [M1;M2];
Psi = {freqresp(Psi1,Psi1.Frequency) freqresp(Psi2,Psi2.Frequency)};

