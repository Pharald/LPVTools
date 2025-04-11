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
% 1. The D = 0 in the state space system P
% 2. All parameter variation is normalised in the range [-1 1]

nin = nargin;
narginchk(1, 2);
if nin==1
    Zb = [];
end


% Parse system data:
if isempty(Zb)
    Zb = basis(1,0);
elseif ~isa(Zb,'struct')
    error('Zb must be a structure with fields Zb.basis and Zb.partial')
end


% Single balancing transformation
% TODO HP 10/11/24: in the original code we did not balance.
% P0 = P;
% P = lpvbalance(P0);


%% =================================================================
% Transform system and basis to non-LPVTools format.
% [Pdata,RBz,BFz,Pz] = basis2data(P,Zb); % EB: 04.25 make lft basis
% definition compatible with pss basis objects

if nnz(strcmp(fieldnames(Zb),'bSplit')) == 0
    nparz = 0;
else
nparz = size(Zb.bSplit,1); % # of parameters in the basis functions
end
RBz = [-1 1]*ones(nparz,1); % EB 04.25: assuming normalised

% Get state space data of PMAT system P
% get state space data
if isa(P,'plftss')
[At, Bt, Ct, Dt] = ssdata(P);
% data as umat objects
A = At.Data;
B = Bt.Data;
C = Ct.Data;
D = Dt.Data;
else
 A = P.a;
 B = P.b;
 C = P.c;
 D = P.d;
end

ny = size(C,1);    % # of measurements
nu = size(B,2);    % # of inputs
nstate = size(A,1);
nbasis = size(Zb.basis,1); % # of basis functions
ybNames = fieldnames(Zb.basis.Uncertainty);
npary = size(ybNames,1); % # of params in basis Y
np = Zb.bSplit*nstate;

R = eye(ny); %eye(ny)+D*D'; % assuming D == 0, for operation R^(0.5) later, R cannot be an lft
S = eye(nu); %eye(nu)+D'*D; 

simplifyopt = 'full'; % EB 31.07: Reduces number of occurences

Hp = Zb.basis(1)*eye(nstate);
partialHp = Zb.partial(1)*eye(nstate);
for ii = 2:nbasis
    Hp = [Hp; Zb.basis(ii)*eye(nstate)];
    partialHp = [partialHp; Zb.partial(ii)*eye(nstate)];
end
Hp_0 = Hp.nominalValue;


%% Factorisation of Parameter dependence from LMI conditions --------------

% define Qz
G_matrix = blkdiag([Hp partialHp; zeros(sum(np),nstate), Hp], -eye(nu + ny));
Q_matrix = [A, B; eye(nstate), zeros(nstate,nu); -C, zeros(ny,nu); zeros(nu,nstate) -eye(nu)];
Qz_matrix = simplify(G_matrix*Q_matrix,simplifyopt);

Gw = [zeros(nstate,nstate) eye(nstate); eye(nstate) zeros(nstate,nstate); zeros(sum(np),nstate), Hp]; 

% Qz
% [QQ,ndeltaq] = partition_mat(Qz_matrix);

% Gw
% [GW,ndeltagw] = partition_mat(Gw);

%% LMI Setup

setlmis([])
k = size(np,1);
if  k > 1
    % blockdiagonal structure, basis functions for each parameter are
    % independent
    [Z_0,n,~] = lmivar(1,[np, ones(k,1)]); % blockdiagonal
else
    [Z_0,n,~] = lmivar(1,[np 1]); % symmetric, npxnp
end

[W,n,sW] = lmivar(1,[nstate 1]);

% defining multipliers
% sS_1w = skewdec(ndeltagw,n);
% [S_1w,n,~] = lmivar(3,sS_1w); % skew symmetric
% sR_1w = diag(n+1:ndeltagw+n);
% [R_1w,n,~] = lmivar(3,sR_1w); % diagonal
% PiGw = lmivar(3,[-sR_1w, sS_1w; sS_1w',sR_1w]);
% 
% sS_1q = skewdec(ndeltaq,n);
% [S_1q,n,sS_1q] = lmivar(3,sS_1q); % skew symmetric
% sR_1q = diag(n+1:ndeltaq+n);
% R_1q = lmivar(3,sR_1q); % diagonal
% PiQz = lmivar(3,[-sR_1q, sS_1q; sS_1q', sR_1q]);

%% LMI conditions
cnt = 1;

% first LMI condition
% 0 < Z_0
lmiterm([-cnt 1 1 Z_0],1,1);  % Z_0
cnt = cnt+1;

% 0 < Z(0)
lmiterm([-cnt 1 1 Z_0],Hp_0',Hp_0);   % H_p'*Z_0*Hp_0
cnt = cnt+1;

% second LMI condition
% full block S procedure multipliers defined in fullBlockS
% R_1q < 0;
% lmiterm([cnt 1 1 R_1q],1,1);  % R_1q
[QQ,PiQz,n,cnt] = fullBlockS(Qz_matrix,n,cnt);
cnt = cnt+1;

% QQ'*blkdiag(PiQz, Z_0mat)*QQ < 0
lmiterm([cnt 0 0 0],QQ);          % QQ'__ QQ outer factor
lmiterm([cnt 1 1 PiQz],1,1);      % PiQz
% Z_mat
lmiterm([cnt 2 3 Z_0],1,1);
lmiterm([cnt 4 4 0],-eye(nu+ny));          % -eye(nu+ny)
cnt = cnt+1;

% inversion condition
% 0 < R_1w
% lmiterm([-cnt 1 1 R_1w],1,1); % R_1w
[GW,PiGw,n,cnt] = fullBlockS(Gw,n,-cnt);
cnt = -cnt;
cnt = cnt + 1;

% 0 < GW'*blkdiag(PiGw, W_mat)*GW
lmiterm([-cnt 0 0 0],GW)          % GW'__GW outer factor
lmiterm([-cnt 1 1 PiGw],1,1);     % PiGw
% W_mat
lmiterm([-cnt 2 3 0],1);
lmiterm([-cnt 3 3 W],1,1);
lmiterm([-cnt 4 4 Z_0],1,1);

% % % ------------------------------------------------------------------------
% % positive definiteness of W
% % 0 < W(0)
% lmiterm([7 1 1 0],0);               % 0
% lmiterm([-7 1 1 W],Hp_0',Hp_0);     % H_p'*W_0*Hp_0
% % % ------------------------------------------------------------------------


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
valZ0 = dec2mat(lmisys,xopt,Z_0);
Zp = Hp'*valZ0*Hp;
partialZp = partialHp'*valZ0*Hp + Hp'*valZ0*partialHp;

% Coprime Factorization State Matrices reconstruction
L = -(B*D'+Zp\C')/R; 

Afact = A+L*C;
Bfact = [L, B+L*D];
Cfact = R^(-0.5)*C;
Dfact = [R^(-0.5), R^(-0.5)*D];

Ml = ss(Afact,Bfact(:,1:ny),Cfact,Dfact(:,1:ny));
Nl = ss(Afact,Bfact(:,ny+1:end),Cfact,Dfact(:,ny+1:end));
fact =ss(Afact,Bfact,Cfact,Dfact);

end 

function [RR1,ndelta] = partition_mat(G)

% ------------Partition R1 -----------------------------------------------
if isa(G,'uss') || isa(G,'umat')
G = simplify(G,'full');
    [G0,delta] = lftdata(G);

elseif isa(G,'plftss') || isa(G,'plftmat')
    G = simplify(G.Data,'full');
[G0,delta,~,~] = lftdata(G);
end

% lft(G0,deltar); 

n_in = size(G,2);
n_out = size(G,1);
ndelta = size(delta,1);

G011 = G0(1:ndelta,1:ndelta);
G021 = G0(ndelta+1:end,1:ndelta);
G022 = G0(ndelta+1:end,ndelta+1:end);
G012 = G0(1:ndelta,ndelta+1:end);

RR1 = [G011, G012; eye(ndelta), zeros(ndelta,n_in); G021, G022];
end
