function [fact,Ml,Nl,L] = lpvlccf(P,Zb)
% LPVLCCF  Compute left contractive coprime factorization for an PLFTSS
%
% [FACT,Ml,Nl,L] = LPVLCCF(P,Zb) computes a contractive
% coprime factorization of the LPV system P (analogous to a normalized
% coprime factorization for LTI systems).Zb is a structure with fields 
% 'basis' and 'partial' containning PLFT objects describing the assumed parameter
% dependence of the solution matrix of the generalised filtering inequality.
% Zb.basis(i) and Zb.partial(i) describe the ith parameter dependency (basis function).
% If there is more than one parameter, an additional field Zb.bSplit is a
% nx1 vector describing the number of basis functions for n parameters.
% FACT = LPVLCCF(P,Zb) returns a state-space realization FACT of [Ml,Nl].
% [FACT,Ml,Nl,L] = LPVLCCF(P,Zb) also returns the coprime factors Ml, Nl and
% the gain L separately.
%
% The following assumptions are made in the current implementation:
% 1. The D matrix in the plant is not parameter dependent


nin = nargin;
narginchk(1, 2);
if nin==1
    Zb = [];
end

if ~isa(P,'plftss')
    error('P must be a plftss object')
end

nx = order(P);

% Parse system data:
if isempty(Zb)
    Zb = basis(plftmat(eye(nx)),plftmat(zeros(nx,nx)));
elseif ~isa(Zb,'basis')
    error('Zb must be a BASIS object')
end

npar = length(fieldnames(Zb.BasisFunction.Parameter)); % # parameters in basis function
% Check for non-rate bounded case
ratebndflg = true;
if npar==0
    ratebndflg = false;
end

% Assign LFT for basis function and partials
Hp = Zb.BasisFunction;
Hp0 = Hp.Data.nominalvalue;
partialHp = [];
if ratebndflg
    partialHp = Zb.Partials;
end


% Dimensions
szP = size(P);
ny = szP(1);    % # of measurements
nu = szP(2);    % # of inputs 
np = size(Hp,1);

simplifyopt = 'full'; % EB 31.07: Reduces number of occurences but is an approximation

% State-space data
[A, B, C, D] = ssdata(P);


% check for parameter dependence in D
nparD = length(fieldnames(D.Parameter));
if nparD ~=0 
    error('D must be a constant matrix')
end

R = eye(ny)+D*D'; R = R.Data.NominalValue;
S = eye(nu)+D'*D; S = S.Data.NominalValue;
Atil = A-B*inv(S)*D'*C;
% Ctil = C'*inv(R)*C;



%% Factorisation of Parameter dependence from LMI conditions --------------

% % define Qz
% G_matrix = blkdiag([Hp partialHp; zeros(sum(np),nstate), Hp], -eye(nu + ny));
% Q_matrix = [A, B; eye(nstate), zeros(nstate,nu); -C, zeros(ny,nu); zeros(nu,nstate) -eye(nu)];
% Qz_matrix = simplify(G_matrix*Q_matrix,simplifyopt);
% 
% Gw = [zeros(nstate,nstate) eye(nstate); eye(nstate) zeros(nstate,nstate); zeros(sum(np),nstate), Hp]; 

Qz = [Hp*Atil + partialHp, Hp*B; ...
      Hp, zeros(np,nu); ...
      C, zeros(ny,nu); ...
      zeros(nu,nx), eye(nu)];
Qz = simplify(Qz,simplifyopt);
                 

Gw = simplify([zeros(nx,nx) eye(nx); ... 
       eye(nx) zeros(nx,nx); ...
       zeros(np,nx), Hp]);


%% LMI Setup

setlmis([])

Z0 = lmivar(1,[np 1]);

[W,~,sW] = lmivar(1,[nx 1]);

%% LMI conditions
cnt = 1;

% % first LMI condition
% % 0 < Z_0
% lmiterm([-cnt 1 1 Z_0],1,1);  % Z_0
% cnt = cnt+1;

% 0 < Z(0)
lmiterm([-cnt 1 1 Z0],Hp0',Hp0);   % H_p'*Z_0*Hp_0
cnt = cnt+1;

% second LMI condition
% full block S procedure multipliers defined in fullBlockS
[QQ,PiQz,~,cnt] = fullBlockS(Qz,cnt);
cnt = cnt+1;

% QQ'*blkdiag(PiQz, Z_0mat)*QQ < 0
lmiterm([cnt 0 0 0],QQ);          % QQ'__ QQ outer factor
lmiterm([cnt 1 1 PiQz],1,1);      % PiQz
% Z_mat
lmiterm([cnt 2 3 Z0],1,1);
lmiterm([cnt 4 4 0],-inv(R));
lmiterm([cnt 5 5 0],-eye(nu));          % -eye(nu+ny)
cnt = cnt+1;

% inversion condition
[GW,PiGw,~,cnt] = fullBlockS(Gw,-cnt);
cnt = -cnt;
cnt = cnt + 1;

% 0 < GW'*blkdiag(PiGw, W_mat)*GW
lmiterm([-cnt 0 0 0],GW)          % GW'__GW outer factor
lmiterm([-cnt 1 1 PiGw],1,1);     % PiGw
% W_mat
lmiterm([-cnt 2 3 0],1);
lmiterm([-cnt 3 3 W],1,1);
lmiterm([-cnt 4 4 Z0],1,1);

% Get LMI Options
% TODO HP 10/11/2024: Shall we include options object in the code?

% Default settings for LMI Lab
LMIopt = zeros(5,1);
LMIopt(2) = 350;  % Max # of iters for rate bounded syn
LMIopt(5) = 1;        % Toggle display

% Solve LMI
lmisys = getlmis;
FeasFlag = 1;

% Create objective function:
ndec = decnbr(lmisys); 
cobj = zeros(ndec,1);
cobj(diag(sW)) = 1; % minimize Trace of W

[copt,xopt] = mincx(lmisys,cobj,LMIopt);
if isempty(copt)
    FeasFlag = 0;
end

% Handle Infeasible LMI Case
if ~FeasFlag
    fact = [];
    Ml =[];
    Nl = [];
    L =[];
    return;
end

%%
valZ0 = dec2mat(lmisys,xopt,Z0);
Zp = Hp'*valZ0*Hp;
% partialZp = partialHp'*valZ0*Hp + Hp'*valZ0*partialHp;

% Coprime Factorization State Matrices reconstruction
L = -(B*D'+Zp\C')/R; 

Afact = A+L*C;
Bfact = [L, B+L*D];
Cfact = R^(-0.5)*C;
Dfact = [R^(-0.5), R^(-0.5)*D];


Ml = ss(Afact,Bfact(:,1:ny),Cfact,Dfact(:,1:ny));
Nl = ss(Afact,Bfact(:,ny+1:end),Cfact,Dfact(:,ny+1:end));
fact = ss(Afact,Bfact,Cfact,Dfact);


end 


