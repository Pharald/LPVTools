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

% All state space matrices and Hp must be umat (or double) so that
% lftdata() correclty extracts the delta blocks after the outer factors are
% constructed

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
% both cases possible with plftmat object --------------------------------
% zbNames = fieldnames(Zb.basis.Uncertainty); % if params are ureal
zbNames = fieldnames(Zb.basis.Data.Uncertainty); % if params are tvreal
% ------------------------------------------------------------------------
nparz = size(zbNames,1); % # of params in basis Z
if nparz > 1
    np = Zb.bSplit*nstate;
else
    np = size(Zb.basis,1)*nstate;
end
RB = Zb.basis.RateBounds;
RBz = Zb.basis.RateBounds{2};

% not compatible with parameter varying D
if ~isa(D,'double') && ~isempty(fieldnames(D.Uncertainty))
    error('D matrix in plant cannot be parameter dependent')
end

R = eye(ny) + D*D';
S = eye(nu) + D'*D; 
Ab = A-B*inv(S)*D'*C;

simplifyopt = 'full'; % EB 31.07: Reduces number of occurences

Hp = Zb.basis.Data(1)*eye(nstate);

if ~isa(Zb.partial,'double')
    flgplft = 1;
    partialHp = Zb.partial.Data(1)*eye(nstate);
else
    flgplft = 0;
    partialHp = Zb.partial(1)*eye(nstate);
end

% stack basis functions
for ii = 2:nbasis
    Hp = [Hp; Zb.basis.Data(ii)*eye(nstate)];
    if flgplft
        partialHp = [partialHp; Zb.partial.Data(ii)*eye(nstate)];
    else
        partialHp = [partialHp; Zb.partial(ii)*eye(nstate)];
    end
end

partialHp = partialHp*RBz(2);
Hp_0 = Hp.NominalValue;


%% Factorisation of Parameter dependence from LMI conditions --------------

% % define Qz
% G_matrix = blkdiag([Hp partialHp; zeros(sum(np),nstate), Hp], -eye(nu + ny));
% Q_matrix = [A, B; eye(nstate), zeros(nstate,nu); -C, zeros(ny,nu); zeros(nu,nstate) -eye(nu)];
% Qz_matrix = simplify(G_matrix*Q_matrix,simplifyopt);
% 
% Gw = [zeros(nstate,nstate) eye(nstate); eye(nstate) zeros(nstate,nstate); zeros(sum(np),nstate), Hp]; 

Qz_matrix = simplify([Hp*Ab + partialHp, Hp*B; ...
                     Hp, zeros(size(Hp,1),size(B,2)); ...
                     -C, zeros(size(C,1),size(B,2)); ...
                  1/2*(R\C), zeros(size(C,1),size(B,2)); ...
                 zeros(size(B,2),size(C,2)), -eye(size(B,2)); ...
                 zeros(size(D,1),size(Ab,2)), -D],simplifyopt);

Gw = simplify([zeros(nstate,nstate) eye(nstate); ... 
       eye(nstate) zeros(nstate,nstate); ...
       zeros(np,nstate), Hp]);


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

%% LMI conditions
cnt = 1;

% % first LMI condition
% % 0 < Z_0
% lmiterm([-cnt 1 1 Z_0],1,1);  % Z_0
% cnt = cnt+1;

% 0 < Z(0)
lmiterm([-cnt 1 1 Z_0],Hp_0',Hp_0);   % H_p'*Z_0*Hp_0
cnt = cnt+1;

% second LMI condition
% full block S procedure multipliers defined in fullBlockS
[QQ,PiQz,n,cnt] = fullBlockS(Qz_matrix,n,cnt);
cnt = cnt+1;

% QQ'*blkdiag(PiQz, Z_0mat)*QQ < 0
lmiterm([cnt 0 0 0],QQ);          % QQ'__ QQ outer factor
lmiterm([cnt 1 1 PiQz],1,1);      % PiQz
% Z_mat
lmiterm([cnt 2 3 Z_0],1,1);
lmiterm([cnt 4 5 0],eye(ny));
lmiterm([cnt 6 6 0],-eye(nu+ny));          % -eye(nu+ny)
cnt = cnt+1;

% inversion condition
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
if isa(R,'double')
    Cfact = R^(-0.5)*C;
    Dfact = [R^(-0.5), R^(-0.5)*D];
else

    Cfact = lftinvsqrt(R)*C;
    Dfact = [lftinvsqrt(R), lftinvsqrt(R)*D];

end

Ml = ss(Afact,Bfact(:,1:ny),Cfact,Dfact(:,1:ny));
Nl = ss(Afact,Bfact(:,ny+1:end),Cfact,Dfact(:,ny+1:end));
fact = ss(Afact,Bfact,Cfact,Dfact);

% --- PLFTSS ---
Ml = plftss(Ml,RB);
Nl = plftss(Nl,RB);
fact = plftss(fact,RB);

end 


function y = lftinvsqrt(R)

% R^(-0.5)
% "quick fix" using symbolic toolbox

% sqrtR = R^(0.5);
[Rb,deltar,blkstruct] = lftdata(R);
Uname = fieldnames(R.Uncertainty);

nr = size(Rb,1);
nd = size(deltar,1);
np = size(blkstruct,1); % # of uncertain params

for iv = 1:np
    sId = blkstruct(iv).Occurrences; % size of param in delta block 
end

urealVec = cell([1 np]);

if nd~=0
    Rb11 = Rb(1:nd,1:nd);
    Rb12 = Rb(1:nd,nd+1:end);
    Rb21 = Rb(nd+1:end,1:nd);
    Rb22 = Rb(nd+1:end,nd+1:end);

    delta = sym('p1')*eye(sId(1));
    deltaVec = 'p1';
    urealVec{1} = R.Uncertainty.(Uname{1});

    if np > 1 
        for ii = 2:nb
            name = append('p',char(string(ii)));
            deltaVec = [deltaVec, name];
            delta = blkdiag(delta,sym(name)*eye(sId(ii)));
        end
    end


Rsym = Rb22 + Rb21*delta*inv(eye(nd) - Rb11*delta)*Rb12;
ysym = inv(sqrtm(Rsym));

y = subs(ysym,deltaVec,urealVec);


% ----------------------------------
% yb22 = sqrt(Rb22);
% 
% y = inv(sqrtR);

else
    y = R.NominalValue^(-0.5);
end

end
