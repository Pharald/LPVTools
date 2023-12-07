function S = IQCtvrealnrb()
S.IQCfunction = @LOCALIQCtvrealnrb;
S.IQCparams = [];
S.PsiFlag = true;


function [lmiout,IDcenter,Psi] = LOCALIQCtvrealnrb(lmisys,ncopies,omeg,params)
% Assume tvreal is unit norm-bounded with no rate bound 

% Grab LMI info
setlmis(lmisys)
ndec = decnbr(lmisys);
nlmis = lminbr(lmisys);

% Define matrix variables and filter for the norm bound constraint
[X1,ndec,X1dec] = lmivar( 3, symdec(ncopies,ndec) );
if ncopies > 1
    [Y1,ndec,Y1dec] = lmivar( 3, skewdec(ncopies,ndec) );
    [M1,ndec,M1dec] = lmivar( 3, [X1dec Y1dec; Y1dec' -X1dec] );
else
    [M1,ndec,M1dec] = lmivar( 3, [X1dec 0; 0 -X1dec] );
end
nlmis = nlmis+1;
lmiterm([-nlmis 1 1 X1],1,1);
Psi1 = eye(2*ncopies);

% Pack up outputs
lmiout = getlmis;
IDcenter = M1;
Psi = { freqresp(ss(Psi1),omeg) };

