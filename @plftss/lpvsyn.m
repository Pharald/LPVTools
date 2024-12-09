
function [Kopt,gamopt,Info] = lpvsyn(P,nmeas,ncont,Xb,Yb,opt)
% LPVSYN  Parameter-dependent controller synthesis for PLFTSS
%
%
% See also: lpvsynOptions.

% sysb = P; % no balancing occurs
narginchk(5,6);
nin = nargin;


if nin<=4
    % EB 06.12.24: How is no basis function case handled?
%     if nin==3
%         % K = lpvsyn(sys,nmeas,ncont)
%         opt = lpvsynOptions;
%     else
%         % K = lpvsyn(sys,nmeas,ncont,opt)
%         opt = Xb;
%     end
%     Yb = [];
%     Xb = Yb;
elseif nin==5
    % K = lpvsyn(sys,nmeas,ncont,Xb,Yb)
    opt = lpvsynOptions;
elseif nin~=6
    error('Incorrect number of input arguments.')
end

Method = opt.Method;

% Dimensions
szP = iosize(P);
nu = ncont;
ne = nmeas;
nd = nu; % disturbance is summative on input

%%

sysb = P; % no balancing

params = fieldnames(sysb.Parameter); % 09.12.24: Error sysb.Parameter no longer exists?
for ii = 1:length(params)
    domRange(ii,:) = sysb.Parameter.(params{ii}).range;
    domRB(ii,:) = sysb.Parameter.(params{ii}).ratebounds;
    domData{ii,:} = linspace(domRange(ii,1),domRange(ii,2),30); % 30 grid points per param
end
rgridPlaceholder = rgrid(params,domData',domRB);

% Orthogonalize general OLIC interconnection
[sysb,TL,TR,FT] = orthog4syn_plftss(sysb,rgridPlaceholder,nmeas,ncont);
%%

[M,DELTA,BLKSTRUCT,NORMUNC] = lftdata(sysb,[],'Parameters');
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

simplifyopt = 'full'; % EB 31.07: Reduces number of occurences in lft blocks

%% Basis functions
% define basis functions for feedback
nbasisX = size(Xb.basis,1); % # of basis functions X
xbNames = fieldnames(Xb.basis.Data.Uncertainty);
nparx = size(xbNames,1); % # of params in basis X
Gp = Xb.basis(1)*eye(nx);
partialGp = Xb.partial(1)*eye(nx);
for ii = 2:nbasisX
    Gp = [Gp; Xb.basis(ii)*eye(nx)];
    partialGp = [partialGp; Xb.partial(ii)*eye(nx)];
end

Gp_0 = Gp.Data.nominalValue;
np_g = size(Gp,1);

% define basis functions for filter
nbasisY = size(Yb.basis,1); % # of basis functions Y
ybNames = fieldnames(Yb.basis.Data.Uncertainty);
npary = size(ybNames,1); % # of params in basis Y
Hp = Yb.basis(1)*eye(nx);
partialHp = Yb.partial(1)*eye(nx);

for ii = 2:nbasisY
    Hp = [Hp; Yb.basis(ii)*eye(nx)];
    partialHp = [partialHp; Yb.partial(ii)*eye(nx)];
end

Hp_0 = Hp.Data.nominalValue;
np_h = size(Hp,1);

%%

% State-space data
[A,B,C,D] = ssdata(sysb);

% partitioned
% General form following notation in thesis of Wu
C11 = C(1:ne,:);
C12 = C(ne+1:ne+nu,:);
C2 = C(ne+nu+1:ne+nu+ne,:); 

B11 = B(:,1:nd);
B12 = B(:,nd+1:nd+ne);
Bof2 = B(:,ne+nd+1:ne+nd+nu);

D1111 = D(1:ne,1:nd);
D1112 = D(1:ne,nd+1:nd+ne);
D1121 = D(ne+1:ne+nu, 1:nd);
D1122 = D(ne+1:ne+nu,nd+1:nd+ne);

B1 = [B11 B12];
B2 = Bof2;
C1 = [C11; C12];
D111dot = [D1111 D1112];
D112dot = [D1121 D1122];
Dmat = [D111dot; D112dot];
D11dot1 = [D1111; D1121];
D11dot2 = [D1112; D1122];
D11 = [D111dot; D112dot];
D12 = D(1:ne+nu, nd+ne+1:end);
D21 = D(ne+nu+1:end, 1:nd+ne);
D22 = D(ne+nu+1:end, nd+ne+1:end);

Ahat = A- B2*C12;
Abar = A-B12*C2;

Bhat = B1 - B2*D112dot;
Cbar = C1 - D11dot2*C2;


%% LMI matrices
RpX_G_col1 = [eye(sum(np_g)); zeros(sum(np_g))]*Gp; % need two repeats of Gp as they are diagonal
RpX_G_col2 = [zeros(sum(np_g)); eye(sum(np_g))]*Gp; 
RpX_G = [RpX_G_col1,RpX_G_col2];
RpX_mult = [Ahat' , C11', zeros(nx,nd+ne); eye(nx), zeros(nx,2*ne + nd)];
RpX_rest =   [B2', zeros(nu,2*ne+nd); zeros(ne,nx), eye(ne), zeros(ne,ne+nd); zeros(nx,nx+ne), Bhat; eye(nx), zeros(nx,2*ne+nd);zeros(ne+nd,nx+ne),eye(ne+nd); zeros(ne+nd,nx), D111dot', zeros(ne+nd,ne+nd)];
RpX_Gpart = [-partialGp, zeros(sum(np_g),nd+2*ne); zeros(sum(np_g) + 4*ne+2*nd+2*nx, nx+2*ne+nd)];

RpXmat = [RpX_G*RpX_mult; RpX_rest] + RpX_Gpart;

RpY_H_col1 = [eye(sum(np_h)); zeros(sum(np_h))]*Hp;
RpY_H_col2 = [zeros(sum(np_h)); eye(sum(np_h))]*Hp; 
RpY_H = [RpY_H_col1,RpY_H_col2];
RpY_mult = [Abar, B11, zeros(nx,nu+ne); eye(nx), zeros(nx,nd +ne + nu)];
RpY_rest =  [C2 zeros(ne,nd +ne + nu); zeros(nd, nx), eye(nd), zeros(nd,ne+nu); zeros(nx,nx+ne), Cbar'; eye(nx), zeros(nx,nd+ne+nu); zeros(ne+nu,nx+nd), eye(ne+nu);zeros(ne+nu,nx), D11dot1, zeros(ne+nu,ne+nu)];
RpY_Hpart = [partialHp, zeros(sum(np_h),nd+nu+ne); zeros(sum(np_h) + nd + 2*nx + 3*ne + 2*nu, nx+nd+nu+ne)];

RpYmat = [RpY_H*RpY_mult; RpY_rest] + RpY_Hpart;

GHxy = [Gp zeros(sum(np_g),nx); zeros(nx,nx), eye(nx); eye(nx), zeros(nx,nx); zeros(sum(np_h),nx), Hp];


% cast outer factors into LFTs
[RRY,ndeltary] = partition_mat(RpYmat);
[RRX,ndeltarx] = partition_mat(RpXmat);
[RRXY,ndeltaxy] = partition_mat(GHxy);

nxy = size(RRXY,2);

% Create LMI variables
setlmis([])

kx = size(np_g,1);
if kx > 1
    [X_0,~,sX_0] = lmivar(1,[np_g ones(kx,1)]);  % blkdiag, npxnp
else
    [X_0,~,sX_0] = lmivar(1,[np_g 1]);  % symmetric, npxnp
end

ky = size(np_h,1);
if ky > 1
[Y_0,ndec,sY_0] = lmivar(1,[np_h ones(ky,1)]);  % blkdiag, npxnp
else
[Y_0,ndec,sY_0] = lmivar(1,[np_h 1]);  % symmetric, npxnp
end

% mult 1
sSx = skewdec(ndeltarx,ndec);
[Sx,ndec,sSx] = lmivar(3,sSx);         % skew symmetric
sRx = diag(ndec+1:ndeltarx+ndec);
[Rx,ndec,~] = lmivar(3,sRx);                 % diagonal
PiRX = lmivar(3,[-sRx, sSx; sSx', sRx]);

% mult 2
sSy = skewdec(ndeltary,ndec);
[Sy,ndec,sSy] = lmivar(3,sSy);         % skew symmetric
sRy = diag(ndec+1:ndeltary+ndec);
[Ry,ndec,~] = lmivar(3,sRy);                 % diagonal
PiRY = lmivar(3,[-sRy, sSy; sSy', sRy]);

% mult 3
sSxy = skewdec(ndeltaxy,ndec);
[Sxy,ndec,sSxy] = lmivar(3,sSxy);      % skew symmetric
sRxy = diag(ndec+1:ndeltaxy+ndec);
Rxy = lmivar(3,sRxy);               % diagonal
PiXY = lmivar(3,[-sRxy, sSxy; sSxy', sRxy]);


[gam,ndec] = lmivar(1,[1 1]);

if isequal(opt.Method,'MaxFeas')
    [FV,ndec] = lmivar(1,[1 1]);
% [FVone,ndec,sFVone] = lmivar(1,[1 1]);
% [FV,ndec,sFV] = lmivar(3,sFVone*eye(nxy))
end

% LMI conditions
cnt = 0;

if nnz(sqrt(Dmat.data.nominal'*Dmat.data.nominal) > eye(ne + nu)) > 0
    % Dmat is plftmat but should have no blocks so equal to nominal
warning('Dmat^T*Dmat^0.5 > I, first condition does not hold')
end
% 

% 0 < gam
cnt = cnt + 1;
lmiterm([-cnt 1 1 gam],1,1);

% Gammalb*I < gam < Gammaub*I
if opt.Gammalb>0
    cnt = cnt+1;
    lmiterm([cnt 1 1 0],opt.Gammalb);
    lmiterm([-cnt 1 1 gam],1,1);
end
if isfinite(opt.Gammaub)
    cnt = cnt+1;
    lmiterm([-cnt 1 1 0],opt.Gammaub);    
    lmiterm([cnt 1 1 gam],1,1);
end


% X and Y positive definite
% 0 < Y_0
cnt = cnt + 1;
lmiterm([-cnt 1 1 Y_0],1,1);

% 0 < Hp_0'*Y_0*Hp_0
cnt = cnt + 1;
lmiterm([-cnt 1 1 Y_0],Hp_0',Hp_0);

% 0 < X_0
cnt = cnt + 1;
lmiterm([-cnt 1 1 X_0],1,1);

% 0 < Gp_0'*X_0*Gp_0
cnt = cnt + 1;
lmiterm([-cnt 1 1 X_0],Gp_0',Gp_0);

% 0 < Rxy
cnt = cnt + 1;
lmiterm([-cnt 1 1 Rxy],1,1);

%     0 < RRXY'*blkdiag(PiXY,XY_0mat)*RRXY
cnt = cnt + 1;
lmiterm([-cnt 0 0 0],RRXY);  % RRXY'__RRXY outer factor
lmiterm([-cnt 1 1 PiXY],1,1);
lmiterm([-cnt 2 2 X_0],1,1);
lmiterm([-cnt 3 4 0],eye(nx));
lmiterm([-cnt 5 5 Y_0],1,1);
% if isequal(Method,'MaxFeas')
% lmiterm([cnt 0 0 FV],eye(nxy));
% end
if isequal(Method,'MaxFeas')
   cnt = cnt+1;
    lmiterm([-cnt 1 1 FV],eye(np_h),eye(np_h));
    lmiterm([cnt 1 1 X_0],1,1);
    cnt = cnt+1;
    lmiterm([-cnt 1 1 FV],eye(np_g),eye(np_g));
    lmiterm([cnt 1 1 Y_0],1,1);
end

% Rx < 0
cnt = cnt + 1;
lmiterm([cnt 1 1 Rx],1,1);

% RRX'*blkdiag(PiRX,X_0mat)*RRX < 0
cnt = cnt + 1;
lmiterm([cnt 0 0 0],RRX);  % RRX'__RRX outer factor
lmiterm([cnt 1 1 PiRX],1,1); % PiRX multiplier
% X_0mat
lmiterm([cnt 2 3 X_0],1,1);
lmiterm([cnt 4 4 gam],-eye(nu),1);
lmiterm([cnt 5 5 gam],-eye(ne),1);
lmiterm([cnt 6 7 0],eye(nx));
lmiterm([cnt 8 8  gam],-eye(ne + nd),1);
lmiterm([cnt 8 9 0],eye(ne + nd));

% Ry < 0
cnt = cnt + 1;
lmiterm([cnt 1 1 Ry],1,1);

% RRY'*blkdiag(PiRY,Y_0mat)*RRY < 0
cnt = cnt + 1;
lmiterm([cnt 0 0 0],RRY);  % RRY'__RRY outer factor
lmiterm([cnt 1 1 PiRY],1,1); % PiRY multiplier
% Y_0mat
lmiterm([cnt 2 3 Y_0],1,1);
lmiterm([cnt 4 4 gam],-eye(ne),1);
lmiterm([cnt 5 5 gam],-eye(nd),1);
lmiterm([cnt 6 7 0],eye(nx));
lmiterm([cnt 8 8 gam],-eye(ne + nu),1);
lmiterm([cnt 8 9 0],eye(ne + nu));



% Set Objective
 cobj= zeros(ndec,1);

 if isequal(opt.Method,'MaxFeas')
     % Method = 'MaxFeas':   min bound on X_0 Y_0 
    cobj(end) = 1;
 else
     % Method = 'BackOff' or 'MinGamma'
     cobj(end) = 1;
 end

% Get LMI Options
% TODO PJS 5/29/2011: Default options for other solvers?
if ~isempty(opt.SolverOptions)
    LMIopt = opt.SolverOptions;
elseif isequal(opt.Solver,'lmilab')
    % Default settings for LMI Lab
    LMIopt = zeros(5,1);
    %LMIopt(1) = 1/(gmax-gmin); % Tol setting in old code
    if isequal(Method,'MaxFeas')
        LMIopt(2) = 50;   % Max # of iters for MaxFeas problem
%     elseif RateBndFlag
%         %LMIopt(2) = 600;  % Setting in old LPVOFSYN1
%         LMIopt(2) = 350;  % Max # of iters for rate bounded syn
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
FeasFlag = 1;
if isequal(Method,'PoleCon')
    % Pole constraint is a GEVP
    % Currently only implemented using LMILAB/GEVP
    nlfc = nmod;
    [copt,xopt] = gevp(lmisys,nlfc,LMIopt);
    
    % TODO: Add initial conditions
    %[copt,copt] = gevp(lmisys,nlfc,LMIopt,t0,x0);
elseif isequal(opt.Solver,'lmilab')
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
    K = [];
    gamopt = inf;
    Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys);
    return;
end

% Get optimal X/Y/Gamma variables
% if ~isequal(Method,'PoleCon')
    gamopt = dec2mat(lmisys,xopt,gam);
% end

   Yp = Hp'*dec2mat(lmisys,xopt,Y_0)*Hp;
   partialYp = partialHp'*dec2mat(lmisys,xopt,Y_0)*Hp + Hp'*dec2mat(lmisys,xopt,Y_0)*partialHp;
   Xp = Gp'*dec2mat(lmisys,xopt,X_0)*Gp;
   partialXp = partialGp'*dec2mat(lmisys,xopt,X_0)*Gp + Gp'*dec2mat(lmisys,xopt,X_0)*partialGp;

    Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys);

if isequal(opt.Method,'BackOff')
    % Solve Stage 2 LMI: Maximize the Min Eig of X/Y Coupling Constraint   
    opt2 = opt;
    opt2.Method = 'MaxFeas';
    opt2.Gammaub = opt.BackOffFactor*gamopt;
      
    % Construct feasible solution from optimal Stage 1 LMI answer
    % TO DO
%     x0 = [xopt];
    opt2.SolverInit = []; %x0;
    
    % Solve Stage 2 Relaxed LMI
    Info1 = Info;
    gamopt1 = gamopt;    
   
    [Kopt,gamopt,Info2] = lpvsyn(P,nmeas,ncont,opt2);
%  lmisys = setmvar(lmisys,gam,opt2.Gammaub);
% [tfeas,xfeas] = feasp(lmisys,[0 1000 0 10 0]);

    Info = struct('MinGamma',gamopt1,'Stage1Info',Info1,'Stage2Info',Info2);
    return
else
%     Construct controller (see thesis of Wu)

% 02.12.2024 EB: there is a controller construction function klpv()
% is this the same reconstruction?

 % the following controller formulation uses the Y and X defined in the
 % thesis of Wu, the above optimisation rearranged to have non-inverse
 % gamma terms -> this has to be reversed for the next step
     YP = 1/gamopt*Yp;
     XP = 1/gamopt*Xp;
     partialXP = 1/gamopt*partialXp;

g2 = gamopt^(2);
g_2 = 1/g2;

XPinv = XP\eye(size(XP,1));
partialXPinv = -XP\(partialXP/XP);

Q = YP-g_2*XPinv;

Om = - D1122 - D1121*((g2*eye(size(D1111,2))\eye(size(D1111,2))) - D1111'* D1111)*D1111'*D1112;

A_ = A + B2*Om*C2;
B1_ = B1 + B2*Om*D21;
C1_ = C1 + D12*Om*C2;
D11_ = D11 + D12*Om*D21;

Dh = (eye(size(D11_,1)) - g_2*(D11_*D11_'))\eye(size(D11_,1));
Dt = (eye(size(D11_,2)) - g_2*(D11_'*D11_))\eye(size(D11_,2));

    F = (-D12'*Dh*D12)\((B2+B1_*D11_'*Dh*D12*g2)'/XP + D12'*Dh*C1_);
    L = -(YP\(C2+D21*Dt*D11_'*C1_*g_2)' + B1_*Dt*D21')/(D21*Dt*D21');

    Af = A_ + B2*F;
    Cf = C1_ + D12*F;

    Afx = XP\Af;
    left = XP\B1_ + Cf'*D11_;
    
    H = -(Afx+Afx'+partialXPinv+Cf'*Cf+g_2*left*Dt*left');

    q = YP - g_2*XPinv;
    qiy = q\YP;

        m1 = H + F'*(B2'/XP + D12'*Cf);
        m2 = (q*(-qiy*L*D21 - B1_) + g_2*F'*D12'*D11_)*Dt*left';
        m = m1 + m2;  

        %% 02.12.24 EB: a way to reduce number of occurences than using simplify
        
        K.a = Af + qiy*L*C2 - g_2*(q\m);
        K.b = -qiy*L;
        K.c = F;
        K.d = Om;
        Kopt = plftss(K.a,K.b,K.c,K.d);
     
        Info = struct('xopt',xopt,'lmisys',lmisys,...
               'Xopt',Xp,'Yopt',Yp);

end
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

% check partition
if size(G011) ~= [ndelta ndelta] | size(G021) ~= [n_out ndelta] | size(G012) ~= [ndelta n_in] | size(G022) ~= [n_out n_in]
    error('dimensions of partitions incorrect')
end

RR1 = [G011, G012; eye(ndelta), zeros(ndelta,n_in); G021, G022];
end





