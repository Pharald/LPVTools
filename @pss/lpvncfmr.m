function [Pred,info] = lpvncfmr(P,varargin)
%LPVNCFMR balanced truncation model reduction through contractive
% coprime factorization of a PSS.
% 
% [Pred,INFO] = LPVNCFMR(P,ORDER) finds a balanced contractive
% coprime factorization of the LPV system P (analogous to a normalized
% coprime factorization for LTI systems), and performs a balanced 
% truncation to remove those states that contribute the least to the 
% input-output mapping of the balanced LPV system. If P has M states, then 
% the reduced order system Pred will have ORDER states, with M-ORDER states 
% removed using the balanced truncation. INFO.hsv contains a vector of 
% singular values describing the input-output mapping of SYSB (comparable 
% to Hankel singular values for LTI systems).
%
% [Pred,INFO] = LPVNCFMR(P) computes the balanced contractive coprime 
% factorization of the LPV system P, and outputs it as Pred. 
% This is equivalent to the call [Pred,INFO] = LPVNCFMR(P,Nx) for a 
% system P with Nx states.
%
% See also	lpvbalancmr, lpvbalreal, lpvgram, modred. 


% TODO - AH 3/20/14
%      - Implement a call where the function loops over different
%        orders of reduction. This is needed to implement 'MaxError'
%        call. It is also needed to replicate functionality of standard
%        MATALB ncfmr which if called as e.g. ncfmr(G,[10:2:18]);
%        will return an array of reduced order models 
%        with 10,12,14,16, and 18 states respectively. 
% TODO - AH 3/20/14
%      - Implement 'MaxError' and 'Display' input arguments.
% TODO - AH 3/20/14
%      - Allow user to specify options for mincx and the optimization of 
%        the cost function J = trace(Wc*Wo). 
% TODO - AH/3/20/14 - See below
%      - Compute coprime factors and store in INFO.GL & INFO.GR

% Parse system data:
[A,B,C,D] = ssdata(P);
nX = size(A,1);

% Single balancing transformation
P0 = P;
P = lpvbalance(P0);

% Define stopping conditions:
maxiter = 10;
reltol = 1e-3;

if nargin == 1
    order = nX;
elseif nargin ==2
    order = varargin{1};
else
    if isa(varargin{1},'double')
        order = varargin{1};
        varargin(1) =[];
    end
    
    % XXX do more input parsing here to grab name/value pairs
    % This can be copied basically from lpvbalancmr except we can't
    % use weights here
end

% BALANCE MODEL:

% Implements coprime factor model reduction in Sec 7.5 of Woods Thesis
% (LMIs manipulated to allow for implentation of algorithm 7.7.1)
Wo = eye(nX);
Woi = eye(nX);
go = 1;
k = 1;
J = zeros(maxiter,1);
while go == 1
    
    % Compute Gramians
    [Wc,Wci] = LOCALncfgram(P,'c',Woi);
    [Wo,Woi] = LOCALncfgram(P,'o',Wci);
    
    % Compute cost function
    J(k) = trace(Wc*Wo);
    
    % Check stopping conditions
    if k == 1
        reldiff = 2*reltol;
    else
        reldiff = abs( (J(k)-J(k-1))/J(k) );
    end
    if (k>=maxiter) || (reldiff<reltol)
        go = 0;
    end
    k = k+1;
end

% Construct transformation using Moore's algorithm
[T,Ti,G] = LOCALbalT(Wo,Wc);

% Compute balanced original system:
Pbal = ss2ss(P,T);

% MODEL REDUCTION

% XXX do a loop on order, if order is a vector creating multiple reduced
% models for different levels of order reduction
% info.ErrorBound becomes a vector
% Pred becomes an array of pss
%
% XXX Compute coprime factors and store in INFO.GL & INFO.GR
elim = (order+1):nX;
Pred = modred(Pbal,elim,'Truncate');
info.hsv = G;




function [W,Wi] = LOCALncfgram(P,type,Wgt)

if isequal(type,'c')
    cflag = true;
elseif isequal(type,'o')
    cflag = false;
    Wci = Wgt;
end
Wgt = 0.5*(Wgt+Wgt');

% Get state space data
[A,B,C,D] = ssdata(P);
nX = size(A,1);
nU = size(B,2);
nY = size(C,1);

S = eye(nU)+D'*D;
R = eye(nY)+D*D';

% Create LMI variables
setlmis([]);
[Xidec,ndec,Xivar] = lmivar(1,[nX 1]);

% Create lmis at each point in the domain
Dom = P.DomainPrivate;
PIVName = P.DomainPrivate.IVName;
npts = prod(Dom.LIVData);
cnt = 1;
for i=1:npts
    % Get parameter value
    if isempty(Dom)
        pvaluec = cell(0,1);
    else
        pvaluec = num2cell(Dom(i));  % single index into RGRID gives value
    end
    
    % Evaluate state matrices at parameter
    % TODO PJS 5/16/2011: Replace with LPVSUBS? LPVSUBS returns data as a
    % double but it does an interpolation rather than a pure extraction.
    AV = double(lpvsplit(A,PIVName,pvaluec));
    BV = double(lpvsplit(B,PIVName,pvaluec));
    CV = double(lpvsplit(C,PIVName,pvaluec));
    DV = double(lpvsplit(D,PIVName,pvaluec));
    RV = double(lpvsplit(R,PIVName,pvaluec));
    SV = double(lpvsplit(S,PIVName,pvaluec));
    
    % Define gramian LMI
    if cflag
        % LMI: Xi*Ahat'+Ahat*Xi-B*inv(S)*B'+Xi*C'*inv(R)*C*Xi <0
        SiV = SV\eye(size(SV));
        AhatV = AV-BV*SiV*DV'*CV;
        lmiterm([cnt 1 1 Xidec],1,AhatV','s'); 
        lmiterm([cnt 1 1 0],-BV*SiV*BV'); 
        lmiterm([cnt 1 2 Xidec],1,CV');        
        lmiterm([cnt 2 2 0],-RV);        
    else
        % LMI: (Xi-Wci)*Atil+Atil'*(Xi-Wci) - C'*inv(R)*C +
        %            (Xi-Wci)*B*inv(S)*B'*(Xi-Wci)<0
        RiV = RV\eye(size(RV));
        AtilV = AV-BV*DV'*RiV*CV;
        lmiterm([cnt 1 1 Xidec],1,AtilV,'s'); 
        lmiterm([cnt 1 1 0],-CV'*RiV*CV-Wci*AtilV-AtilV'*Wci); 
        lmiterm([cnt 1 2 Xidec],1,BV);
        lmiterm([cnt 1 2 0],-Wci*BV);
        lmiterm([cnt 2 2 0],-SV);        
    end
    cnt = cnt+1;
    
end

% tolL*I < Xi
tolL = 0;
lmiterm([-cnt 1 1 Xidec],1,1);
lmiterm([cnt 1 1 0],tolL*eye(nX));
cnt = cnt+1;

% Xi < tolU*I 
% (This prevents X from becoming too large due to small HSV)
tolU = 1e10;
lmiterm([cnt 1 1 Xidec],1,1);
lmiterm([-cnt 1 1 0],tolU*eye(nX));
cnt = cnt+1;

% Objective: min -trace(Wgt*Xi) = min c'xvar
cobj = zeros(ndec,1);
cobj(diag(Xivar)) = diag(Wgt);
L = tril(Xivar,-1);
ldx = find(L);
cobj(L(ldx)) = 2*Wgt(ldx);
cobj = -cobj;

% Get LMI Options
% XXX - LMI Options? Solver Options?
LMIopt = [1e-4 200 0 0 1];

% Solve LMI
lmisys = getlmis;
[copt,xopt] = mincx(lmisys,cobj,LMIopt);
Wi = dec2mat(lmisys,xopt,Xidec);
W = inv(Wi);

% XXX The implementation of Moore's algorithm appears in several
% of the LPV model reduction functions. This should be implemented
% as a separate utility function.
function [T,Ti,G] = LOCALbalT(Wo,Wc)

% Balance System ( Moore's algorithm, using notation in section 9.5 in
% Robust Systems by Sanchez-Pena and Sznaier)
[Vc,Sc] = LOCALschursort(Wc);
Scrt = sqrt(diag(Sc));
T1 = lrscale(Vc',1./Scrt,[]); %T1 = inv(Scrt)*Vc'
T1i = lrscale(Vc,[],Scrt); %T1i = Vc*Scrt
Wotil = (T1')\(Wo/T1);
[Vo,So] = LOCALschursort(Wotil);
Sort = diag(So).^(1/4);
T2 = lrscale(Vo',Sort,[]); %T2 = Scrt*Vo'
T2i = lrscale(Vo,[],1./Sort); %T2i = Vo*inv(Scrt)

% Form final transformation
T = T2*T1;
Ti = T1i*T2i;

% Compute Hankel singular values:
G = sqrt(diag(So));




function [V,S] = LOCALschursort(W)
[V,S] = schur(W);
[Ssort,idx] = sort(diag(S),'descend');
S = diag(Ssort);
V = V(:,idx);







