function [F,Gamma,Info] = lpvsfsyn(P,nu,Xb,opt)
% Parameter-dependent state feedback controller synthesis in LFT
% formualtion
% P is uss and is the generalised plant with weighted in/outputs.
% it must not include the state feedback signal to the controller (or dimensions will not work)
% nu = number of control inputs from controller to plant

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
Gp0 = Gp.Data.nominalvalue;
np = size(Gp,1);
partialGp = zeros(np,nx);
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


Qx = [Gp*Ahat' - partialGp, Gp*C1', zeros(np,nd);...
    Gp, zeros(np,ne-nu+nd);...
    B2', zeros(nu,ne-nu+nd);...
    zeros(ne-nu,nx), eye(ne-nu), zeros(ne-nu,nd);...
    zeros(nx,nx+ne-nu), Bhat; ...
    eye(nx), zeros(nx,ne-nu+nd); ...
    zeros(nd,nx+ne-nu),eye(nd);...
    zeros(nd,nx), D11', zeros(nd,nd)]; 
Qx = simplify(Qx,simplifyopt);


% Initialize LMI and define variables
setlmis([])

if isequal(Method,'MaxFeas')
    % LBC is lower bound on X
    [LBC,ndec] = lmivar(1,[nx 0]);
end

[ginv,ndec] = lmivar(1,[1 1]);

X = lmivar(1,[np 1]);    

% LMI for Dissipation Inequality
cnt = 1;
[QQ,PiQx,ndec,cnt] = fullBlockS(Qx,cnt);
cnt = cnt + 1;

% QQ'*blkdiag(PiQx,X_0mat)*QQ <0
lmiterm([cnt 0 0 0],QQ);            % QQ'__QQ outer factor
lmiterm([cnt 1 1 PiQx],1,1);        % PiQx

% X_0mat
lmiterm([cnt 2 3 X],1,1);
lmiterm([cnt 4 4 0],-eye(ne));
lmiterm([cnt 5 6 ginv],eye(nx),1);
lmiterm([cnt 7 7 0],-eye(nd));
lmiterm([cnt 7 8 ginv],eye(nd),1);
cnt = cnt + 1;

% X_0 > 0
lmiterm([-cnt 1 1 X],1,1);
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
cnt = cnt+1;

if isfinite(opt.Xub)
    xpdupp = opt.Xub;
else
    xpdupp = 1e6;
end
lmiterm([-cnt 1 1 0],xpdupp*eye(nx));
lmiterm([cnt 1 1 X],Gp0',Gp0);

cnt = cnt+1;

% constraint [-I D1'/gam; D1/gam -I] <0
% which implies that gam>= max(svd(D1)), i.e. ginv<=1/max(svd(D1))
gmin = max(0,opt.Gammalb);

% check if D1 is parameter dependent
nparD1 = length(fieldnames(D1.Parameter));
if nparD1 == 0 
    gmin = max(gmin,norm(D1.Data.NominalValue));
else
   D1grid = lft2grid(D1,100);
   normD1 = norm(D1grid);
   gmin = max(gmin,max(normD1.Data));
end  

if gmin>0
    lmiterm([-cnt 1 1 0],1/gmin);
    lmiterm([cnt 1 1 ginv],1,1);
    cnt=cnt+1;
end

% 1/Gammaub <= 1/Gamma
if isfinite(opt.Gammaub)
    lmiterm([-cnt 1 1 ginv],1,1);
    lmiterm([cnt 1 1 0],1/opt.Gammaub);
    cnt = cnt +1;
end


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

% SDP: min gam subject to LMI constraints
lmisys = getlmis;
ndec = decnbr(lmisys);
c = zeros(ndec,1);
c(1) = -1; 


if ~isequal(opt.Solver,'lmilab')
    error('Specified solver is currently not available.');
end

[copt,xopt] = mincx(lmisys,c,LMIopt,x0);
Info = [];

% Handle Infeasible LMI Case
if isempty(xopt)
    F = [];
    Gamma = inf;
    Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys);
    return;
else
    Gamma = 1/dec2mat(lmisys,xopt,ginv);
    if isequal(Method,'BackOff')
        opt2 = opt;
        opt2.Method = 'MaxFeas';
        opt2.Gammaub = opt.BackOffFactor*Gamma;
        opt2.SolverInit = [xopt;0];
        info1 = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,'X',X);
        gam1 = Gamma;
        [F,Gamma,info2] = lpvsfsyn(P,nu,Xb,opt2);
        Info = struct('MinGamma',gam1,'Stage1Info',info1,'Stage2Info',info2);

    else
        % Controller reconstruction
        Xp = Gp'*dec2mat(lmisys,xopt,X)*Gp;
        Zopt= inv(Xp)/Gamma^2;

        % Explicit solution for controller reconstruction
        p12t = C1*Xp + D11*Bhat'/Gamma^2;
        p22 = -eye(ne-nu) + D11*D11'/Gamma^2;
        p23t = D12*D11'/Gamma^2;
        F = -( B2'+ D12*Bhat'/Gamma^2 ...
            -p23t*(p22\p12t) )*Zopt*Gamma^2 - C2;

        % Undo orthog transformation
        F = r2inv*F;

    end
end

    

end