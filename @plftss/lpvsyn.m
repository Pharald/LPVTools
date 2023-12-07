function [Kopt,gamopt,Info] = lpvsyn(P,nmeas,ncont,opt)
% LPVSYN  Parameter-dependent controller synthesis for PLFTSS
%
% [K,GAM,INFO] = LPVSYN(P,NMEAS,NCON) computes a parameter-varying
% controller K which minimizes the L2 norm of the interconnection 
% defined by lft(P,K). K is a PLFTSS with NMEAS inputs and NCON outputs. 
% GAM is the L2 norm of lft(P,K). INFO is a structure containing data from 
% the Linear Matrix Inequalities that are solved to obtain K.
%
% [K,GAM,INFO] = LPVSYN(P,NMEAS,NCON,OPT) allows the user to pass in
% a LPVSYNOPTIONS object. 
%
% The default algorithm for LPVSYN will solve the given synthesis problem
% twice. The first iteration attempts to find a solution that minimizes the
% induced L2 norm of lft(P,K). The second iteration will solve the 
% optimization problem again, with the caveat that any solution that is 
% found to lie within 15% of the optimal induced L2 norm of lft(P,K) from 
% the first iteration, is satisfactory. This formulation has been found to 
% yield controllers that are better numerically conditioned. The back-off 
% factor of 15% can be set in LPVSYNOPTIONS.
%
% See also: lpvsynOptions.



% PJS 12/31/2012: Initial Implementation
%
% References
%  1. Gain Scheduling via LFTs by Packard, Sys. and Control Letters, 1994.
%  2. A Convex Characterization of Gain-Scheduled Hinfty Controllers,
%       by Apkarian and Gahinet, IEEE TAC, 1995.
%  3. Erratum to "A Convex Characterization of Gain-Scheduled Hinfty
%       Controllers," by Apkarian and Gahinet, IEEE TAC, 1995, p1681.
%  4. Explicit Controller Formulas for LMI-based Hinfty Synthesis,
%       by Gahinet, Automatica, 1996.
%
% This function essentially implements the LMI conditions in Ref 2 for
% continuous-time LPV synthesis (Thm 5.1). Note Ref 3 corrects an error
% in LMI conditions given in Ref. The conditions in the
% reference are for a robust performance condition:
%      min gam subject to  Induced Gain <= gam for |theta(t)|<=1/gam
% This function slightly modifies the LMI conditions to solve the
% related worst-case gain condition with normalized parameters:
%      min gam subject to  Induced Gain <= gam for |theta(t)|<=1
%

%%
% Error Checking
narginchk(3,4);
if nargin == 3
    opt = lpvsynOptions;
end

%%
% Initialize Outputs
Kopt = [];
gamopt = [];
Info = [];

%% Parse the lpvsynOptions

Method = opt.Method;

%%
% Problem Dimensions
szP = iosize(P);
ny = nmeas;
nu = ncont;
ne = szP(1) -  ny;
nd = szP(2) - nu;

[M,DELTA,BLKSTRUCT,NORMUNC] = lftdata(P,[],'Parameters');
Nblk = length(BLKSTRUCT);
nri = zeros(Nblk,1);
for i=1:Nblk
    if ~strcmp( BLKSTRUCT(i).Type , 'tvreal' )
        error('Plant may only have tvreal blocks');
    else
        % XXX PJS: This assumes that BLKSTRUCT and NORMUNC have the
        % same ordering for the blocks (with repeats in NORMUNC)
        nri(i) = BLKSTRUCT(i).Occurrences;
    end
end
nr = sum(nri);
nx = order(M);

%% Balance the system
% % XXX Balance I/O channels associated with parameters?
% Min = M;
% M = ssbal(M);

%%
% Unpack Data: Follows notation in Eq 2.16 of Ref 2
[A,B,C,D] = ssdata(M);
it = 1:nr;
i1d = nr+1:nr+nd;
i2u = nr+nd+1:nr+nd+nu;
i1e = nr+1:nr+ne;
i2y = nr+ne+1:nr+ne+ny;

Bt = B(:,it);
B1 = B(:,i1d);
B2 = B(:,i2u);

Ct = C(it,:);
C1 = C(i1e,:);
C2 = C(i2y,:);

Dtt = D(it,it);
Dt1 = D(it,i1d);
Dt2 = D(it,i2u);

D1t = D(i1e,it);
D11 = D(i1e,i1d);
D12 = D(i1e,i2u);

D2t = D(i2y,it);
D21 = D(i2y,i1d);
D22 = D(i2y,i2u);

%%
% XXX PJS: Handle ill-conditioned D12 and D21
% (See comment 1 in Section 3.5 of Ref 4)


%%
% Compute optimal L2 gain bound
% Follows notation in Thm 5.1 of Ref 2 except that gam is the induced
% L2 gain bound assuming the parameter range is normalized to [-1,1].
NR = null([B2; Dt2; D12]');
NR = blkdiag(NR,eye(nr+nd));

NS = null([C2 D2t D21]);
NS = blkdiag(NS,eye(nr+ne));

B1hat = [Bt B1];
C1hat = [Ct;C1];
D11hat = [Dtt Dt1; D1t D11];


% Set up the LMIs
setlmis([]);

R = lmivar(1,[nx 1]);
S = lmivar(1,[nx 1]);
if nr>0
    J3 = lmivar(1,[nri(:) ones(Nblk,1)]);
    L3 = lmivar(1,[nri(:) ones(Nblk,1)]);
end
[gam,ndec] = lmivar(1,[1 1]);

if isequal(Method,'MaxFeas');
    % FV is an upper bound on R, S, J3 and L3.
    [FV,ndec] = lmivar(1,[1 1]);
end

% LMI in R (Eqn 5.2 in Ref 2 has an error. See Ref 3)
lmiterm([1 0 0 0],NR);
tmp = [eye(nx) zeros(nx,nr+ne)];
lmiterm([1 1 1 R],[A; C1hat],tmp,'s');
tmp = [zeros(nx,ne); zeros(nr,ne); eye(ne)];
lmiterm([1 1 1 gam],-tmp,tmp');
tmp = [zeros(nx,nr) B1; zeros(nr+ne,nr) [Dt1; D11]];
lmiterm([1 1 2 0],tmp);
tmp = [zeros(nr,nd); eye(nd)];
lmiterm([1 2 2 gam],-tmp,tmp');
if nr>0
    tmp = [zeros(nx,nr); eye(nr); zeros(ne,nr)];
    lmiterm([1 1 1 J3],-tmp,tmp');
    tmpl = [Bt; Dtt; D1t];
    tmpr = [eye(nr) zeros(nr,nd)];
    lmiterm([1 1 2 J3],tmpl,tmpr);
    tmp = [eye(nr); zeros(nd,nr)];
    lmiterm([1 2 2 J3],-tmp,tmp');
end

% LMI in S (Eqn 5.3 in Ref 2 has an error. See Ref 3)
lmiterm([2 0 0 0],NS);
tmp = [eye(nx) zeros(nx,nr+nd)];
lmiterm([2 1 1 S],[A'; B1hat'],tmp,'s');
tmp = [zeros(nx,nd); zeros(nr,nd); eye(nd)];
lmiterm([2 1 1 gam],-tmp,tmp');
tmp = [zeros(nx,nr) C1'; zeros(nr+nd,nr) [D1t'; D11']];
lmiterm([2 1 2 0],tmp);
tmp = [zeros(nr,ne); eye(ne)];
lmiterm([2 2 2 gam],-tmp,tmp');
if nr>0
    tmp = [zeros(nx,nr); eye(nr); zeros(nd,nr)];
    lmiterm([2 1 1 L3],-tmp,tmp');
    tmpl = [Ct'; Dtt'; Dt1'];
    tmpr = [eye(nr) zeros(nr,ne)];
    lmiterm([2 1 2 L3],tmpl,tmpr);
    tmp = [eye(nr); zeros(ne,nr)];
    lmiterm([2 2 2 L3],-tmp,tmp');
end

% (R,S) Coupling LMI (Eq 5.4)
lmiterm([-3 1 1 R],1,1);
lmiterm([-3 2 1 0],eye(nx));
lmiterm([-3 2 2 S],1,1);

% (L3,J3) Coupling LMI (Eq 5.5)
cnt = 4;
if nr>0
    lmiterm([-4 1 1 L3],1,1);
    lmiterm([-4 2 1 0],eye(nr));
    lmiterm([-4 2 2 J3],1,1);
    cnt = 5;
end


% Xlb*I < R < Xub*I
if opt.Xlb>0
    lmiterm([cnt 1 1 0],opt.Xlb*eye(nx));
    lmiterm([-cnt 1 1 R],1,1);
    cnt = cnt+1;
end
if isfinite(opt.Xub)
    lmiterm([-cnt 1 1 0],opt.Xub*eye(nx));    
    lmiterm([cnt 1 1 R],1,1);
    cnt = cnt+1;
end

% Ylb*I < S < Yub*I
if opt.Ylb>0
    lmiterm([cnt 1 1 0],opt.Ylb*eye(nx));
    lmiterm([-cnt 1 1 S],1,1);
    cnt = cnt+1;
end
if isfinite(opt.Yub)
    lmiterm([-cnt 1 1 0],opt.Yub*eye(nx));    
    lmiterm([cnt 1 1 S],1,1);
    cnt = cnt+1;
end

% Jlb*I < J3 < Jub*I
if nr> 0 && opt.Jlb>0
    lmiterm([cnt 1 1 0],opt.Jlb*eye(nr));
    lmiterm([-cnt 1 1 J3],1,1);
    cnt = cnt+1;
end
if nr> 0 && isfinite(opt.Jub)
    lmiterm([-cnt 1 1 0],opt.Jub*eye(nr));    
    lmiterm([cnt 1 1 J3],1,1);
    cnt = cnt+1;
end

% Llb*I < L3 < Lub*I
if nr> 0 && opt.Llb>0
    lmiterm([cnt 1 1 0],opt.Llb*eye(nr));
    lmiterm([-cnt 1 1 L3],1,1);
    cnt = cnt+1;
end
if nr> 0 && isfinite(opt.Lub)
    lmiterm([-cnt 1 1 0],opt.Lub*eye(nr));    
    lmiterm([cnt 1 1 L3],1,1);
    cnt = cnt+1;
end

% Gammalb*I < gam < Gammaub*I
if opt.Gammalb>0
    lmiterm([cnt 1 1 0],opt.Gammalb);
    lmiterm([-cnt 1 1 gam],1,1);
    cnt = cnt+1;
end
if isfinite(opt.Gammaub)
    lmiterm([-cnt 1 1 0],opt.Gammaub);    
    lmiterm([cnt 1 1 gam],1,1);
    cnt = cnt+1;
end


% Max feasability constraint = minimize bound on R, S, J3 and L3.
if isequal(Method,'MaxFeas')
    lmiterm([-cnt 1 1 FV],eye(nx),eye(nx));
    lmiterm([cnt 1 1 R],1,1);
    cnt = cnt+1;
    lmiterm([-cnt 1 1 FV],eye(nx),eye(nx));
    lmiterm([cnt 1 1 S],1,1);
    cnt = cnt+1;
    if nr>0
        lmiterm([-cnt 1 1 FV],eye(nr),eye(nr));
        lmiterm([cnt 1 1 J3],1,1);
        cnt = cnt+1;
        lmiterm([-cnt 1 1 FV],eye(nr),eye(nr));
        lmiterm([cnt 1 1 L3],1,1);
        cnt = cnt+1;        
    end
    
end

% Set Objective:  
cobj = zeros(ndec,1);
if isequal(Method,'MaxFeas')   
    % Method = 'MaxFeas':   minimize FV
    cobj(end) = 1;    
    
    % Add small penalty to gamma to force gamma near its lower bound
    % XXX What to choose for this penalty?
    cobj(end-1) = 1e-8;
else
    % Method = 'MinGamma' or 'BackOff'
    cobj(end) = 1;    
end

% Get LMI Options
% TODO PJS 5/29/2011: Default options for other solvers?
if ~isempty(opt.SolverOptions)
    LMIopt = opt.SolverOptions;
elseif isequal(opt.Solver,'lmilab')
    % Default settings for LMI Lab
    LMIopt = zeros(5,1);
    if isequal(Method,'MaxFeas')
        LMIopt(2) = 40;   % Max # of iters for MaxFeas problem
    else
        LMIopt(2) = 250;  % Max # of iters for non-rate bounded syn
    end
    LMIopt(5) = 1;        % Toggle display
else
    LMIopt = [];
end

% Get LMI Initial Condition
if ~isempty(opt.SolverInit)
    x0 = opt.SolverInit;
else
    x0 = [];
end


% Solve LMI
lmisys = getlmis;
% [gamopt,xopt]=mincx(lmisys,cobj,LMIopt);
FeasFlag = 1;
if isequal(opt.Solver,'lmilab')
    [copt,xopt] = mincx(lmisys,cobj,LMIopt,x0);
    if isempty(copt)
        FeasFlag = 0;
    end
elseif isequal(opt.Solver,'sedumi')
    % Convert to Sedumi format
    % TODO PJS 5/29/2011: Currently ignores x0
    [F0,Fi,blk] = lmitrans(lmisys);
    K.s = blk;
    [xdual,xopt,info]=sedumi(Fi,-cobj,F0,K,LMIopt);
    copt = cobj'*xopt;    
    if info.pinf==1 || info.dinf==1 || info.numerr~=0
        % TODO PJS 5/29/2011: Also check info.feasratio?
        FeasFlag = 0;
    end    
else
    % TODO PJS 5/20/2011: Implement other solvers with LMITRANS
    error('Specified solver is currently not available.');
end

% Handle Infeasible LMI Case
% TODO PJS 5/20/2011: What should we return in this case?
if ~FeasFlag
    Kopt = [];
    gamopt = inf;
    Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,'Ropt',[],...
                  'Sopt',[],'Lopt',[],'Jopt',[]);
    return;
end

gamopt=dec2mat(lmisys,xopt,gam);
Ropt=dec2mat(lmisys,xopt,R);
Sopt=dec2mat(lmisys,xopt,S);
if nr>0
    L3opt=dec2mat(lmisys,xopt,L3);
    J3opt=dec2mat(lmisys,xopt,J3);
else 
    L3opt = [];
    J3opt = [];
end

Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,...
               'Ropt',Ropt,'Sopt',Sopt,'Jopt',J3opt,'Lopt',L3opt);


if isequal(Method,'BackOff');
    % Solve Stage 2 LMI: Maximize the Min Eig of X/Y Coupling Constraint   
    opt2 = opt;
    opt2.Method = 'MaxFeas';
    opt2.Gammaub = opt.BackOffFactor*gamopt;
      
    % Construct feasible solution from optimal Stage 1 LMI answer
    FV0 = 1.1*max( eig( blkdiag(Ropt,Sopt,L3opt,J3opt) ) );
    x0 = [xopt; FV0];
    opt2.SolverInit = x0;
    
    % Solve Stage 2 Relaxed LMI
    Info1 = Info;
    gamopt1 = gamopt;    
    [Kopt,gamopt,Info2] = lpvsyn(P,nmeas,ncont,opt2);
    Info = struct('MinGamma',gamopt1,'Stage1Info',Info1,'Stage2Info',Info2);
    return
end


%%
% Controller Reconstruction
Kopt = klpv(M,nmeas,ncont,nri,Ropt,Sopt,J3opt,L3opt,gamopt);
Kopt = lft( Kopt, DELTA );
Kopt = feedback( Kopt, D22 );

%keyboard


return
% Step 1: Compute M,N such that MN' = I-RS
% nk is the order of the controller
[U,Sig,V] = svd( eye(nx) - Ropt*Sopt );
nk = sum( diag(Sig) > tol );
Srt = Sig(1:nk,1:nk).^0.5;
M = U(:,1:nk)*Srt;
N = V(:,1:nk)*Srt;

% Step 2: Compute Xcl
X22 = -N'*Ropt/(M');
Xcl = [Sopt N; N' X22];

% Step 3: Compute L
% nsi is the number of copies of parameter i used by the controller
L1 = cell(Nblk,1);
L2 = cell(Nblk,1);
ptr = 0;
nsi = zeros(Nblk,1);
for i=1:Nblk
    % Index for ith block
    idx = ptr+(1:nri(i));
    ptr = idx(end);
    
    % Compute M,N such that MN' = I-L3*J3 for i^th block
    [U,Sig,V] = svd( eye(nri(i)) - L3opt(idx,idx)*J3opt(idx,idx) );
    nsi(i) = sum( diag(Sig) > tol );
    Srt = Sig(1:nsi(i),1:nsi(i)).^0.5;
    M = U(:,1:nsi(i))*Srt;
    N = V(:,1:nsi(i))*Srt;
    
    % Solve for i^th block of L1 and L2
    L2{i} = N;
    L1{i} = -N'*J3opt(idx,idx)/(M');
end
if nr>0
    ns = sum(nsi);
    L1 = blkdiag( L1{:} );
    L2 = blkdiag( L2{:} );
    L = [L1 L2; L2' L3opt];
else
    ns = 0;
    L = [];
end

% Step 4: Solve Basic LMI  (Eqns 6.3, 6.6-6.8)
% Modify equations to account for possible differences in nsi and nri
calL = -blkdiag(L,gamopt*eye(nd));
calJ = -blkdiag(inv(L),gamopt*eye(ne));

A0 = blkdiag(A,zeros(nk));
B0 = [zeros(nx,ns) Bt B1; zeros(nk,ns+nr+nd)];
calB = [zeros(nx,nk) B2 zeros(nx,ns); eye(nk) zeros(nk,nu+ns)];
C0 = [zeros(ns,nx+nk); Ct zeros(nr,nk); C1 zeros(ne,nk)];
calD11 = blkdiag(zeros(ns),D11hat);
calD12 = [zeros(ns,nk+nu) eye(ns); zeros(nr,nk) Dt2 zeros(nr,ns); ...
    zeros(ne,nk) D12 zeros(ne,ns)];
calC = [zeros(nk,nx) eye(nk); C2 zeros(ny,nk); zeros(ns,nx+nk)];
calD21 = [zeros(nk,ns+nr+nd); zeros(ny,ns) D2t D21; eye(ns) zeros(ns,nr+nd)];

PSI = [A0'*Xcl+Xcl*A0 Xcl*B0 C0'; B0'*Xcl calL calD11'; C0 calD11 calJ];
Q = [calC calD21 zeros(nk+ny+ns,nr+ns+ne)];
PX = [calB' zeros(nk+nu+ns,nr+ns+nd) calD12'];
PX = PX*blkdiag(Xcl,eye(2*nr+2*ns+nd+ne));

OMEGA=basiclmi(PSI,PX,Q);
if isempty(OMEGA)
    Kopt = [];
    return
end

% Pack controller as a PSSLFT (Eq 6.2)
AK = OMEGA(1:nk,1:nk);
BK = OMEGA(1:nk,nk+1:end);
CK = OMEGA(nk+1:end,1:nk);
DK = OMEGA(nk+1:end,nk+1:end);
K = ss(AK,BK,CK,DK);

idx = [];
ptr = 0;
for i=1:Nblk
    % Index for ith block
    idx = [idx ptr+(1:nsi(i))];
    ptr = ptr+nri(i);
end
K = lft( K, DELTA(idx,idx) );

% Add back D22 term from plant
Kopt = feedback( K, D22 );



