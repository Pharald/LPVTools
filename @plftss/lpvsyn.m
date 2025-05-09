function [K,Gamma,Info] = lpvsyn(sys,nmeas,ncont,Xb,Yb,opt)
% LPVSYN  Parameter-dependent controller synthesis for PLFTSS
%
% P is a plftss
% Xb and Yb are structures with fields basis containing the basis functions
% and partial containting the partial derivatives of the basis functions.
% The basis and partial fields must be plft or double.
%
% lpvsyn(P,nmeas,ncont,Xb,Yb,opt)
% lpvsyn(P,nmeas,ncont,opt) % no basis NOT SUPPORTED
% lpvsyn(P,nmeas,ncont,Xb,Yb) % default options
% lpvsyn(P,nmeas,ncont) % default options no basis functions NOT SUPPORTED

% See also: lpvsynOptions.


% Parse Inputs
nin = nargin;
narginchk(3, 6)
nout = nargout;
if nin<=4
    if nin==3
        % K = lpvsyn(sys,nmeas,ncont)
        opt = lpvsynOptions;
    else
        % K = lpvsyn(sys,nmeas,ncont,opt)
        opt = Xb;
    end
    Yb = [];
    Xb = Yb;
elseif nin==5
    % K = lpvsyn(sys,nmeas,ncont,Xb,Yb)
    opt = lpvsynOptions;
elseif nin~=6
    error('Incorrect number of input arguments.')
end
Method = opt.Method;

nx = order(sys); % # of states


if isempty(Yb)
    Yb = basis(plftmat(eye(nx)),plftmat(zeros(nx,nx)));
elseif ~isa(Yb,'basis')
    error('Yb must be a BASIS object')
end
if isempty(Xb)
    Xb = basis(plftmat(eye(nx)),plftmat(zeros(nx,nx)));
elseif ~isa(Xb,'basis')
    error('Xb must be a BASIS object')
end

% Check for non-rate bounded case
nparx = length(fieldnames(Xb.BasisFunction.Parameter)); % # parameters in basis function
npary = length(fieldnames(Yb.BasisFunction.Parameter)); % # parameters in basis function

ratebndflg = true;
if nparx==0 && npary==0
    ratebndflg = false; 
end

% Dimensions
nd = size(sys,2) - ncont;
ne = size(sys,1) - nmeas;

sysb = sys; % no balancing

% Orthogonalize general OLIC interconnection
[sysb,TL,TR,FT] = orthog4syn(sysb,nmeas,ncont);
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

simplifyopt = 'full'; % EB 31.07.24: Reduces number of occurences in lft blocks

%% Basis functions
% define basis functions for feedback
% Assign LFT for basis function and partials
Gp = Xb.BasisFunction;
Gp0 = Gp.Data.nominalvalue;
nXp = size(Gp,1);
partialGp = zeros(nXp,nx);
if ratebndflg
    partialGp = Xb.Partials;
end

Hp = Yb.BasisFunction;
Hp0 = Hp.Data.nominalvalue;
nYp = size(Hp,1);
partialHp = zeros(nYp,nx);
if ratebndflg
    partialHp = Yb.Partials;
end


%%

% State-space data
[a,b,c,d] = ssdata(sysb);
ne1 = ne-ncont;
ne2 = ncont;
nd1 = nd-nmeas;
nd2 = nmeas;

d11 = d(1:ne,1:nd);
d11dot1 = d(1:ne,1:nd1);
d11dot2 = d(1:ne,nd1+1:nd);
d111dot = d(1:ne1,1:nd);
d112dot = d(ne1+1:ne,1:nd);
d1111 = d(1:ne1,1:nd1);
d1112 = d(1:ne1,nd1+1:nd);
d1121 = d(ne1+1:ne,1:nd1);
d1122 = d(ne1+1:ne,nd1+1:nd);
d12 = [zeros(ne1,ne2); eye(ne2,ne2)];
d21 = [zeros(nd2,nd1) eye(nd2,nd2)];

b11 = b(:,1:nd1);
b12 = b(:,nd1+1:nd);
b1 = [b11 b12];
b2 = b(:,nd+1:end);
c11 = c(1:ne1,:);
c12 = c(ne1+1:ne,:);
c1 = [c11;c12];
c2 = c(ne+1:end,:);

Ahat = a - b2*c12;
Bhat = b1 - b2*d112dot;
Atil = a - b12*c2;
Ctil = c1 - d11dot2*c2;



% LMI Outer factors

Qx = [Gp*Ahat' - partialGp, Gp*c11', zeros(nXp,nd);...
    Gp, zeros(nXp,ne1+nd);...
    b2', zeros(ncont,ne1+nd);...
    zeros(ne1,nx), eye(ne1), zeros(ne1,nd);...
    zeros(nx,nx+ne1), Bhat; ...
    eye(nx), zeros(nx,ne1+nd); ...
    zeros(nd,nx+ne1),eye(nd);...
    zeros(nd,nx), d111dot', zeros(nd,nd)]; 
Qx = simplify(Qx,simplifyopt);


Qy = [Hp*Atil + partialHp, Hp*b11, zeros(nYp,ne);...
    Hp, zeros(nXp,nd1+ne);...
    c2, zeros(nmeas,nd1+ne);...
    zeros(nd1,nx), eye(nd1), zeros(nd1,ne);...
    zeros(nx,nx+nd1), Ctil'; ...
    eye(nx), zeros(nx,nd1+ne); ...
    zeros(ne,nx+nd1),eye(ne);...
    zeros(ne,nx), d11dot1, zeros(ne,ne)]; 
Qy = simplify(Qy,simplifyopt);


Qxy = [Gp zeros(nXp,nx); zeros(nx,nx), eye(nx); eye(nx), zeros(nx,nx); zeros(nYp,nx), Hp];
Qxy = simplify(Qxy,simplifyopt);


% Create LMI variables
setlmis([])

if isequal(Method,'PoleCon')
    error('PoleCon Method is currently not supported for plftss objects')
end

if isequal(Method,'MaxFeas')
    % LBC is lower bound on X
    [LBC,ndec] = lmivar(1,[nx 0]);
end

[ginv,ndec] = lmivar(1,[1 1]);

X = lmivar(1,[nXp 1]);    
Y = lmivar(1,[nYp 1]);    


cnt = 1;


% get outer factors and multipliers
% LMIs relating to multipliers also defined in fullBlockS
[RRX,PiRX,ndec,cnt] = fullBlockS(Qx,cnt);
cnt = cnt+1;

[RRY,PiRY,ndec,cnt] = fullBlockS(Qy,cnt);
cnt = cnt+1;

[RRXY,PiXY,ndec,cnt,XYinfo] = fullBlockS(Qxy,-cnt); % XY xondition is 0 <
cnt = -cnt;
cnt = cnt+1;

% LMI conditions

% % X and Y positive definite
% X_0 > 0
lmiterm([-cnt 1 1 X],1,1);
cnt = cnt+1;

% Y_0 > 0
lmiterm([-cnt 1 1 Y],1,1);
cnt = cnt+1;

% ypdlow*I < Hp_0'*Y_0*Hp_0 < ypdupp*I
if opt.Ylb > 0
    ypdlow = opt.Ylb;
else
    ypdlow = 1e-6;
end
lmiterm([cnt 1 1 0],ypdlow*eye(nx));
lmiterm([-cnt 1 1 Y],Hp0',Hp0);
cnt = cnt + 1;

if isfinite(opt.Yub)
    ypdupp = opt.Yub;
else
    ypdupp = 1e6;
end
lmiterm([-cnt 1 1 0],ypdupp*eye(nx));
lmiterm([cnt 1 1 Y],Hp0',Hp0);
cnt = cnt+1;



% xpdlow*I < Gp_0'*X_0*Gp_0 < xpdupp*I
if opt.Xlb > 0
    xpdlow = opt.Xlb;
else
    xpdlow = 1e-6;
end
if isequal(Method,'MaxFeas')
    lmiterm([cnt 1 1 LBC],1,1);
else
    lmiterm([cnt 1 1 0],xpdlow*eye(nx));
end
lmiterm([-cnt 1 1 X],Gp0',Gp0);
cnt = cnt + 1;

if isfinite(opt.Xub)
    xpdupp = opt.Xub;
else
    xpdupp = 1e6;
end
lmiterm([-cnt 1 1 0],xpdupp*eye(nx));
lmiterm([cnt 1 1 X],Gp0',Gp0);
cnt = cnt+1;


% RRX'*blkdiag(PiRX,X_0mat)*RRX < 0
lmiterm([cnt 0 0 0],RRX);  % RRX'__RRX outer factor
lmiterm([cnt 1 1 PiRX],1,1); % PiRX multiplier
% X_0mat
lmiterm([cnt 2 3 X],1,1);
lmiterm([cnt 4 4 0],-eye(ne1+ncont));
lmiterm([cnt 5 6 ginv],eye(nx),1);
lmiterm([cnt 7 7 0],-eye(nd));
lmiterm([cnt 7 8 ginv],eye(nd),1);
cnt = cnt + 1;

% RRY'*blkdiag(PiRY,Y_0mat)*RRY < 0
lmiterm([cnt 0 0 0],RRY);  % RRY'__RRY outer factor
lmiterm([cnt 1 1 PiRY],1,1); % PiRY multiplier
% Y_0mat
lmiterm([cnt 2 3 Y],1,1);
lmiterm([cnt 4 4 0],-eye(nd1+nmeas));
lmiterm([cnt 5 6 ginv],eye(nx),1);
lmiterm([cnt 7 7 0],-eye(ne));
lmiterm([cnt 7 8 ginv],eye(ne),1);
cnt = cnt + 1;


%     0 < RRXY'*blkdiag(PiXY,XY_0mat)*RRXY
lmiterm([-cnt 0 0 0],RRXY);  % RRXY'__RRXY outer factor
lmiterm([-cnt 1 1 PiXY],1,1);
lmiterm([-cnt 2 2 X],1,1);
lmiterm([-cnt 3 4 ginv],eye(nx),1);
lmiterm([-cnt 5 5 Y],1,1);
cnt = cnt + 1;

if opt.Gammalb>0
    lmiterm([cnt 1 1 ginv],1,1);
    lmiterm([-cnt 1 1 0],1/opt.Gammalb);
    cnt = cnt +1;
end
if isfinite(opt.Gammaub)
    lmiterm([-cnt 1 1 ginv],1,1);
    lmiterm([cnt 1 1 0],1/opt.Gammaub);
    cnt = cnt +1;
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
ndec = decnbr(lmisys);
cobj = zeros(ndec,1);
cobj(1) = -1; 

if ~isequal(opt.Solver,'lmilab')
    error('Specified solver is currently not available.');
end

FeasFlag = 1;
[copt,xopt] = mincx(lmisys,cobj,LMIopt,x0);
if isempty(copt)
    FeasFlag = 0;
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
Gamma = 1/dec2mat(lmisys,xopt,ginv);

Yp = Hp'*dec2mat(lmisys,xopt,Y)*Hp;
partialYp = partialHp'*dec2mat(lmisys,xopt,Y)*Hp + Hp'*dec2mat(lmisys,xopt,Y)*partialHp;
Xp = Gp'*dec2mat(lmisys,xopt,X)*Gp;
partialXp = partialGp'*dec2mat(lmisys,xopt,X)*Gp + Gp'*dec2mat(lmisys,xopt,X)*partialGp;

Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys);

if isequal(opt.Method,'BackOff')
    % Solve Stage 2 LMI: Maximize the Min Eig of X/Y Coupling Constraint   
    opt2 = opt;
    opt2.Method = 'MaxFeas';
    opt2.Gammaub = opt.BackOffFactor*Gamma;
      
    % Construct feasible solution from optimal Stage 1 LMI answer
    % TO DO
%     x0 = [xopt];
    opt2.SolverInit = []; %x0;
    
    % Solve Stage 2 Relaxed LMI
    Info1 = Info;
    Gamma1 = Gamma;
   
    [K,Gamma,Info2] = lpvsyn(sys,nmeas,ncont,Xb,Yb,opt2);
%  lmisys = setmvar(lmisys,gam,opt2.Gammaub);
% [tfeas,xfeas] = feasp(lmisys,[0 1000 0 10 0]);

    Info = struct('MinGamma',Gamma1,'Stage1Info',Info1,'Stage2Info',Info2);
    return
else
%     Construct controller (see thesis of Wu)

     g2 = Gamma^(2);
     g_2 = 1/g2;

     XPinv = Xp\eye(size(Xp,1));
     partialXPinv = -Xp\(partialXp/Xp);

     Q = Yp-g_2*XPinv;

     Om = - d1122 - d1121*((g2*eye(size(d1111,2))\eye(size(d1111,2))) - d1111'* d1111)*d1111'*d1112;

     A_ = a + b2*Om*c2;
     B1_ = b1 + b2*Om*d21;
     C1_ = c1 + d12*Om*c2;
     D11_ = d11 + d12*Om*d21;

     Dh = (eye(size(D11_,1)) - g_2*(D11_*D11_'))\eye(size(D11_,1));
     Dt = (eye(size(D11_,2)) - g_2*(D11_'*D11_))\eye(size(D11_,2));

     F = (-d12'*Dh*d12)\((b2+B1_*D11_'*Dh*d12*g2)'/Xp + d12'*Dh*C1_);
     L = -(Yp\(c2+d21*Dt*D11_'*C1_*g_2)' + B1_*Dt*d21')/(d21*Dt*d21');

     Af = A_ + b2*F;
     Cf = C1_ + d12*F;

     Afx = Xp\Af;
     left = Xp\B1_ + Cf'*D11_;

     H = -(Afx+Afx'+partialXPinv+Cf'*Cf+g_2*left*Dt*left');

     q = Yp - g_2*XPinv;
     qiy = q\Yp;

     m1 = H + F'*(b2'/Xp + d12'*Cf);
     m2 = (q*(-qiy*L*d21 - B1_) + g_2*F'*d12'*D11_)*Dt*left';
     m = m1 + m2;

     %% 02.12.24 EB: a way to reduce number of occurences than using simplify

     ak = Af + qiy*L*c2 - g_2*(q\m);
     bk = -qiy*L;
     ck = F;
     dk = Om;
     K = plftss(ak,bk,ck,dk);

     Info = struct('xopt',xopt,'lmisys',lmisys,...
               'Xopt',Xp,'Yopt',Yp);

end


end








