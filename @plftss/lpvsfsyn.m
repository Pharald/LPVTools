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
elseif nin == 3  && ~isa(Xb,struct) % no basis functions specified -> will cause error
    opt = Xb;
    Xb.basis = [];
    Xb.partial = [];
    Xb.bSplit = [];
    error('Basis functions must be specified')
end

Method = opt.Method;

% State-space data
if isa(P,'plftss')
[At, Bt, Ct, Dt] = ssdata(P);
% data as umat objects
Asf = At.Data;
B = Bt.Data;
C = Ct.Data;
D = Dt.Data;
else % also compatible with P as a uss
 Asf = P.a;
 B = P.b;
 C = P.c;
 D = P.d;
end

% Dimensions
nx = size(Asf,1);
ndu = size(B,2);
ne = size(C,1) - nu;
nd = ndu-nu;

% partitioned
% Generalised plant does not include un-scaled in/outputs to controller

simplifyopt = 'full'; % EB 31.07: Reduces number of occurences 

nbasisX = size(Xb.basis,1); % # of basis functions X
if ~isa(Xb.basis,'double')
    Xbb = Xb.basis.Data;
    xbNames = fieldnames(Xbb.Uncertainty);
    nparx = size(xbNames,1); % # of params in basis X
    RB = Xb.basis.RateBounds;
%     RB{i,1} % parameter names
%     RB{i,2} % values, 2x1 double
else
    Xbb = Xb.basis;
    nparx = 0;
    RB = [];
end

if ~isa(Xb.partial,'double')
    Xbp = Xb.partial.Data;
else
    Xbp = Xb.partial;
    nparx = 0;
end


Gp = Xbb(1)*eye(nx);
partialGp = Xbp(1)*eye(nx);
for ii = 2:nbasisX
    Gp = [Gp; Xbb(ii)*eye(nx)];
    partialGp = [partialGp; Xbp(ii)*eye(nx)];
end

Gp_0 = Gp.nominalValue; % Gp is a umat

if nnz(strcmp(fieldnames(Xb),'bSplit')) ~= 0
    np = Xb.bSplit*nx;
else
    np = size(Gp,1);    
end

% Check for non-rate bounded case
% --> not implemented
ratebndflg = true;
if (nbasisX ==1 && nparx==0)
    ratebndflg = false;
end


% the 12 column is empty because there is only 1 weighted input
% EB: 04.25 this is not the general plant form including weighted
% disturbance input, it is the plant form specific to the structured
% controller synthesis
Csf11 = C(1:ne,:);
Csf21 = C(ne+1:ne+nu,:);
Bsf11 = B(:,1:ne);
Bsf12 = [];
Bsf2 = B(:,ne+1:ne+nu);
Dsf1111 = D(1:ne,1:ne);
Dsf1112 = [];
Dsf1121 = D(ne+1:ne+nu, 1:ne);
Dsf1122 = [];
Dsf2 = D(ne+1:ne+nu,nd+1:nd+nu);
Csf2 = eye(nx); % must be full state feedback

% checks
% if  Dsf2~= eye(nu) % check not compatible with umat
%         warning('general form of plant is not adhered to, solution may not work')
% end

%% LFT FORM
% General form following notation of theorem 3 in Theis, Pfifer 2020
Ahat = Asf- Bsf2*Csf21;
B2 = Bsf2;
D111 = [Dsf1111 Dsf1112];
D112 = [Dsf1121 Dsf1122];
Dmat = [D111; D112];
Bhat = [Bsf11 Bsf12]- Bsf2*D112;
C11 = Csf11;
C21 = Csf21;
% -----------------------

% LMI matrices
Qx_pmat = [Gp*Ahat' - partialGp, Gp*C11', zeros(sum(np),ne); Gp, zeros(sum(np),2*ne); B2', zeros(nu,2*ne); zeros(ne,nx), eye(ne), zeros(ne,ne); zeros(nx,nx+ne), Bhat; eye(nx), zeros(nx,2*ne); zeros(ne,nx+ne),eye(ne); zeros(ne,nx), D111', zeros(ne,ne)]; 
Qx_pmat = simplify(Qx_pmat,simplifyopt);

%%

setlmis([])

if isequal(Method,'MaxFeas')
    % lower bound on X
   [LBC,n] = lmivar(1,[nx 0]); 
end

gam = lmivar(1,[1 1]);              % decision variable 1

kx = size(np,1);
if kx > 1
    [X_0,n,~] = lmivar(1,[np ones(kx,1)]);  % blkdiag, npxnp
else
    [X_0,n,~] = lmivar(1,[np 1]);  % symmetric, npxnp
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
    K = plftss(K,RB);
    end
end
end