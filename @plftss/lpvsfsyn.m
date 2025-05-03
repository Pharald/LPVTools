function [K,gamma,info] = lpvsfsyn(P,nu,Xb,opt)
% Parameter-dependent state feedback controller synthesis in LFT
% formualtion
% P is uss and is the generalised plant with weighted in/outputs.
% it must not include the state feedback signal to the controller (or dimensions will not work)
% nu = number of control inputs from controller to plant

%% To Do:
% - updating plant formulation to include disturbance inputs for more
% general case
% - lpvsynoptions not included: lb and ub on X, (Y, J, L), balance
% - solver options -> only lmilab supported
% - rate bounded vs non rate bounded

%%

narginchk(2,4);
nin = nargin;

if nin ==2
    opt = lpvsynOptions; % default settings
    Xb =[];
elseif nin == 3  && ~isa(Xb,'basis') % no basis functions specified -> will cause error
    opt = Xb;
    Xb = [];
end

if ~isa(P,'plftss')
    error(['Plant must be a plftss object'])
end

nx = order(P); % # of states

if isempty(Xb)
    Xb = basis(plftmat(eye(nx)),plftmat(zeros(nx,nx)));
end

npar = length(fieldnames(Xb.BasisFunction.Parameter)); % # parameters in basis function
% Check for non-rate bounded case
ratebndflg = true;
if npar==0
    ratebndflg = false;
end

% Assign LFT for basis function and partials
Gp = Xb.BasisFunction;
partialGp = [];
if ratebndflg
    partialGp = Xb.Partials;
end

Method = opt.Method;



% Dimensions
szP = size(P);
ne = szP(1);   % # of errors
ndu = szP(2);  % # of inputs of G = nd + nu
nd = ndu-nu;   % # of disturbances

simplifyopt = 'full'; % EB 31.07: Reduces number of occurences but is an approximation

% State-space data
[A, B, C, D] = ssdata(P);
B1 = B(:,1:nd);
B2 = B(:,nd+1:end);
D1 = D(:,1:nd);
D2 = D(:,nd+1:end);


% Determine if D2 has full column rank and scale D2 to
%    Q2*D2*R2INV = [0;I].
% Hence if the system is redefined with a unitary transformation
% on ERROR, etilde := Q2 e, and invertible transformation on CONTROL,
% u := R2INV utilde, in the new variables, D2Tilde = [0;I].  Note that
% the unitary transformation on e does not change ||e||, and the invertible
% transformation on u can be included in the overall controller.
%
% NOTE HP 02/05/25: This so far only works if D2 is not parameter
% dependent. 

nparD2 = length(fieldnames(D2.Parameter));

if nparD2 ~=0 
    error(['D2 must be a constant matrix'])
end

D2 = D2.Data.nominalvalue;
[q2,r2] = qr(D2);
rrk = double(rank(r2));
if (rrk ~= nu)
    error(' D12 DOES NOT HAVE FULL COLUMN RANK over IV')
end

q2 = q2(:,[nu+1:end 1:nu]);
r2inv = inv(r2(1:nu,:));  
    
B2 = B2*r2inv;
C = q2'*C;
D1 = q2'*D1;
D2 = [zeros(ne-nu,nu); eye(nu)];
    
C1 = C(1:ne-nu,:);
C2 = C(ne-nu+1:end,:);
D11 = D1(1:ne-nu,:);
D12 = D1(ne-nu+1:end,:);
    
Ahat = A-B2*C2; 
Bhat = B1-B2*D12; 

% Outer factor for full block S-procedure

np = size(Gp,1);
Qx = [Gp*Ahat' - partialGp, Gp*C1', zeros(np,ne);...
    Gp, zeros(np,2*ne-nu);...
    B2', zeros(nu,2*ne-nu);...
    zeros(ne-nu,nx), eye(ne-nu), zeros(ne-nu,ne);...
    zeros(nx,nx+ne), Bhat; ...
    eye(nx), zeros(nx,2*ne-nu); ...
    zeros(ne,nx+ne-nu),eye(ne);...
    zeros(ne-nu,nx), D11', zeros(ne-nu,ne)]; 
Qx = simplify(Qx,simplifyopt);


%%

% Initialize LMI and define variables
setlmis([])
% HP 02/05/25: check if lmis in ginv would not be better
[gam,ndec] = lmivar(1,[1 1]);

X = lmivar(1,[np 1]);    


if isequal(Method,'MaxFeas')
    % LBC is lower bound on X
    [LBC,ndec] = lmivar(1,[nx 0]);
end



%% LMI conditions

cnt = 0;

% first LMI condition
% MOVE THIS CHECK TO THE END and update Gamma value before running
% feasiblity problem for suboptimal
% if nnz(sqrt(Dmat'*Dmat) > eye(ne)) > 0
% warning('Dmat^T*Dmat^0.5 > I, first condition does not hold')
% end

% condition: gam^2*eye(ne) - Dmat'*Dmat >= 0];
% check at the end

if opt.Gammalb > 0
    cnt = cnt+1;
    lmiterm([cnt 1 1 0],opt.Gammalb);
    lmiterm([-cnt 1 1 gam],1,1);
else
    % 0 < gam
cnt = cnt + 1;
lmiterm([-cnt 1 1 gam],1,1);
end
if isfinite(opt.Gammaub)
    cnt = cnt +1;
    lmiterm([cnt 1 1 gam],1,1);
    lmiterm([-cnt 1 1 0],opt.Gammaub);    
end

% 0 < X_0
% cnt = cnt + 1;
% lmiterm([-cnt 1 1 X_0],1,1);
% if isequal(Method,'MaxFeas')
%     lmiterm([cnt 1 1 LBC],1,1);
% end

% 0 < Gp_0'*X_0*Gp_0
cnt = cnt + 1;
lmiterm([-cnt 1 1 X_0],Gp_0',Gp_0);
if isequal(Method,'MaxFeas')
    lmiterm([cnt 1 1 LBC],1,1);
end

% second LMI condition -> handled in fullBlockS
% R_1q < 0
cnt = cnt + 1;
[QQ,PiQx,~,cnt] = fullBlockS(Qx_pmat,n,cnt);

% QQ'*blkdiag(PiQx,X_0mat)*QQ <0
cnt = cnt + 1;
lmiterm([cnt 0 0 0],QQ);            % QQ'__QQ outer factor
lmiterm([cnt 1 1 PiQx],1,1);        % PiQx

% X_0mat
lmiterm([cnt 2 3 X_0],1,1);
lmiterm([cnt 4 4 gam],-eye(nu),1);
lmiterm([cnt 5 5 gam],-eye(ne),1);
lmiterm([cnt 6 7 0],eye(nx));
lmiterm([cnt 8 8 gam],-eye(ne),1);
lmiterm([cnt 8 9 0],eye(ne));

lmisys = getlmis;

% Get LMI Options
if ~isempty(opt.SolverOptions)
    LMIopt = opt.SolverOptions;
elseif isequal(opt.Solver,'lmilab')
    % Default settings for LMI Lab
    LMIopt = zeros(5,1);
    if isequal(Method,'MaxFeas')
        LMIopt(2) = 50;   % Max # of iters for MaxFeas problem
    elseif ratebndflg
        LMIopt(2) = 350;  % Max # of iters for rate bounded syn
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

 n = decnbr(lmisys); 
 c = zeros(n,1);
 c(1) = 1; % corresponding to first decision variable which is gamma

if isequal(Method,'Backoff')
   c(1) = -1; % maximise LBC 
end

if ~isequal(opt.Solver,'lmilab')
    error('Specified solver is currently not available.');
end

[copt,xopt] = mincx(lmisys,c,LMIopt,x0);

info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys);

if isempty(xopt)
    K = [];
    gamma = inf;
    return;
else
    gamma = dec2mat(lmisys,xopt,gam);
    Xp = Gp'*dec2mat(lmisys,xopt,X_0)*Gp;
    partialXp = partialGp'*dec2mat(lmisys,xopt,X_0)*Gp + Gp'*dec2mat(lmisys,xopt,X_0)*partialGp;

    if isequal(Method,'BackOff')
        opt2 = opt;
        opt2.Method = 'MaxFeas';
        opt2.Gammaub = opt.BackOffFactor*gamma;
        opt2.SolverInit = [xopt;0];
        info1 = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,'Xp',Xp,'partialXp',partialXp);
        gam1 = gamma;
        [K,gamma,info2] = lpvsfsyn(P,nu,Xb,opt2);
        info = struct('MinGamma',gam1,'Stage1Info',info1,'Stage2Info',info2);

    else
        
    % build controller
    Xpinv = Xp\eye(size(Xp,1));

    % feedback gain calculation
    K  = -(B2' + D112*Bhat'*gamma^(-2) - inv(D112*D111'*gamma^(-2) - eye(ne))*(C11*Xp*gamma^(-1) + D111*Bhat'*gamma^(-1)))*gamma*Xpinv - C21;
    % reconstruct as a plftss
    K = r2inv*plftmat(K,RB);
    end
end
end