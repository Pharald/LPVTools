% function [F,GAM,INFO] = lpvL2sfsynengine(G,ncont,Fbasis,Fgrad,RateBounds,Opt)
%
% Engine to synthesize a state-feedback controller for grid-based
% LPV systems.
%
% INPUTS
% G: 1-by-1-nmod Array of LTI systems.
% ncont: # of control inputs
% Fbasis: Basis functions used in storage function Nbasis-by-1-by-nmod
% Fgrad: Gradient of the basis functions wrt. the parameters:
%        Nbasis-by-Nparameters-by-nmod
% RateBounds: Upper and lower bounds on the parameter rates: Nparameter-by-2
% Opt: LPVSYNOPTIONS object
%
% OUTPUTS
% F: nu-by-nx-nmod array of state feedback gains
% GAM: Min Hinf cost


% XXX Currently no error checking. The file assumes that all dimensions
% are compatible.
% XXX Add options object to pass through solver options

function [F,gam,info] = lpvL2sfsynengine(G,nu,Fbasis,Fgrad,RateBounds,opt)

% Parse Input
Method = opt.Method;

% Dimensions
G = G(:,:,:);
szG = size(G);
if numel(szG) ==2
    szG = [szG 1];
end
ne = szG(1);   % # of errors
ndu = szG(2);  % # of inputs of G = nd + nu
nmod = szG(3); % # of model in the G array.
nd = ndu-nu;   % # of disturbances
npar = size(RateBounds,1); % # of parameters
nbasis = size(Fbasis,1); % # of basis functions

% Check for non-rate bounded case
ratebndflg = true;
if (nbasis==1 && npar==0)
    ratebndflg = false;
end

% State-space data
[A,B,C,D] = ssdata(G);
nx = size(A,1);
B1 = B(:,1:nd,:);
B2 = B(:,nd+1:end,:);
D1 = D(:,1:nd,:);
D2 = D(:,nd+1:end,:);

% Determine if D2 has full column rank and scale D2 to
%    Q2*D2*R2INV = [0;I].
% Hence if the system is redefined with a unitary transformation
% on ERROR, etilde := Q2 e, and invertible transformation on CONTROL,
% u := R2INV utilde, in the new variables, D2Tilde = [0;I].  Note that
% the unitary transformation on e does not change ||e||, and the invertible
% transformation on u can be included in the overall controller.
r2inv = zeros(nu,nu,nmod);
Ahat = zeros(nx,nx,nmod);
Bhat = zeros(nx,nd,nmod);

C1 = zeros(ne-nu,nx,nmod);
C2 = zeros(nu,nx,nmod);
D11 = zeros(ne-nu,nd,nmod);
D12 = zeros(nu,nd,nmod);
for i=1:nmod
    [q2,r2] = qr(D2(:,:,i));
    rrk = double(rank(r2));
    if (rrk ~= nu)
        error(' D12 DOES NOT HAVE FULL COLUMN RANK over IV')
    end
    q2 = q2(:,[nu+1:end 1:nu]);
    r2inv(:,:,i) = inv(r2(1:nu,:));  
    
    B2(:,:,i) = B2(:,:,i)*r2inv(:,:,i);
    C(:,:,i) = q2'*C(:,:,i);
    D1(:,:,i) = q2'*D1(:,:,i);
    D2(:,:,i) = [zeros(ne-nu,nu); eye(nu)];
    
    C1(:,:,i) = C(1:ne-nu,:,i);
    C2(:,:,i) = C(ne-nu+1:end,:,i);
    D11(:,:,i) = D1(1:ne-nu,:,i);
    D12(:,:,i) = D1(ne-nu+1:end,:,i);
    
    Ahat(:,:,i) = A(:,:,i)-B2(:,:,i)*C2(:,:,i);
    Bhat(:,:,i) = B1(:,:,i)-B2(:,:,i)*D12(:,:,i);
end

% Initialize LMI and define variables
setlmis([])
X = zeros(nbasis,1);
for i1 = 1:nbasis
    X(i1) = lmivar(1,[nx 1]);    
end
%[ginvI,ndec] = lmivar(1,[nd 0]);
[ginv,ndec] = lmivar(1,[1 1]);
ginvI = lmivar(3,ndec*eye(nd));

if isequal(Method,'MaxFeas')
    % LBC is lower bound on X
    [LBC,ndec] = lmivar(1,[nX 0]);
end

% LMI for Dissipation Inequality - Loop through array of models
cnt = 1;
for k1 = 1:nmod
    Ahatk = Ahat(:,:,k1);
    Bhatk = Bhat(:,:,k1);
    B2k = B2(:,:,k1);
    C1k = C1(:,:,k1);
    D11k = D11(:,:,k1);
    
    if ratebndflg
        Fgk1 = Fgrad(:,:,k1);        
    else
        Fgk1 = [];        
    end
    Fbk1 = Fbasis(:,:,k1);
            
    % L2 Bound LMI in block 3-by-3 form with gamma
    for q = 1:2^npar  % number of variables appearing in X's basis fcns        
        for ibasis=1:nbasis
            lmiterm([cnt 1 1 X(ibasis)],[Ahatk; C1k]*Fbk1(ibasis,1),[eye(nx) zeros(nx,ne-nu)],'s');
        end
        lmiterm([cnt 1 1 0],-blkdiag(B2k*B2k',eye(ne-nu)));
        lmiterm([cnt 1 2 ginvI],[Bhatk; D11k],eye(nd));        
        lmiterm([cnt 2 2 0],-eye(nd));        

        rbvec = RateBounds(:,1);
        tmp = dec2bin(q-1,npar);
        idx = find(tmp=='1');
        rbvec(idx) = RateBounds(idx,2);                      
        if ratebndflg
            for k2=1:nbasis
                for k3=1:npar
                    % NOTE PJS 5/15/2011: This term appears to be non-symm
                    % because there is a left factor but no right factor.
                    % However, the left factor is a scalar so this really is
                    % a symmetric term.  This is probably ok but I should just
                    % confirm that LMILab handles this correctly.
                    lmiterm([-cnt 1 1 X(k2)],[eye(nx); zeros(ne-nu,nx)],...
                            0.5*rbvec(k3)*Fgk1(k2,k3)*...
                            [eye(nx) zeros(nx,ne-nu)],'s');
                end
            end
        end
        cnt = cnt+1;
    end
    if ratebndflg || k1 == 1
        % xpdlow*I < X < xpdupp*I
        if opt.Xlb > 0 
            xpdlow = opt.Xlv;
        else
            xpdlow = 1e-6;
        end
        if isequal(Method,'MaxFeas')
            lmiterm([cnt 1 1 LBC],1,1);
        else
            lmiterm([cnt 1 1 0],xpdlow*eye(nx));
        end
        for n1=1:nbasis
            lmiterm([-cnt 1 1 X(n1)],Fbk1(n1),1);
        end
        cnt = cnt+1;
        
        if isfinite(opt.Xub)
            xpdupp = opt.Xub;
        else
            xpdupp = 1e6;
        end
        lmiterm([-cnt 1 1 0],xpdupp*eye(nx));
        for n2=1:nbasis
            lmiterm([cnt 1 1 X(n2)],Fbk1(n2),1);
        end
        cnt = cnt+1;
    end
end

% See Section 4.2 of Wu's Ph. D. Thesis:
% State Feedback posed as R+U*F*V<0 which, by matrix elimination
% lemma, is equivalent to: 
%   Uperp'*R*Uperp<0 and Vperp'*R*Vperp<0
% The Uperp LMI leads to the main LMI constraint. The Vper LMI
% simply leads to a constraint [-I D1'/gam; D1/gam -I] <0
% which implies that gam>= max(svd(D1)), i.e. ginv<=1/max(svd(D1))
gmin = max(0,opt.Gammalb);
for k = 1:nmod
    gmin = max(gmin,norm(D1(:,:,k)));
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
    %LMIopt(1) = 1/(gmax-gmin); % Tol setting in old code
    if isequal(Method,'MaxFeas')
        LMIopt(2) = 50;   % Max # of iters for MaxFeas problem
    elseif RateBndFlag
        %LMIopt(2) = 600;  % Setting in old LPVOFSYN1
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

% SDP: min gamsq subject to LMI constraints
lmisys = getlmis;
c = zeros(ndec,1);
c(end) = -1;

% TODO HP 07/12/2023: so far only lmilab supported. Add other solvers?
if ~isequal(opt.Solver,'lmilab')
    error('Specified solver is currently not available.');
end

[copt,xopt] = mincx(lmisys,c,LMIopt,x0);
gam = 1/xopt(end);
info = [];

% Handle Infeasible LMI Case
if isempty(xopt)
    F = [];
    gam = inf;
    info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys);
    return;
else
    if isequal(Method,'BackOff')
        opt2 = opt;
        opt2.Method = 'MaxFeas';
        opt2.Gammaub = opt.BackOffFactor*gam;
        info1 = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,'X',X,'Y',Y);
        gam1 = gam;
        [F,gam,info2] = lpvsyn(sys,nmeas,ncont,Xb,Yb,opt2);
        info = struct('MinGamma',Gamma1,'Stage1Info',info1,'Stage2Info',info2);

    else
        % Controller reconstruction
        F = zeros(nu,nx,nmod);
        Zopt = zeros(nx,nx,nmod);
        for k1=1:nmod
            Fbk1 = Fbasis(:,:,k1);
            Xk1 = zeros(nx);
            for ibasis=1:nbasis
                Xk1 = Xk1 + Fbk1(ibasis,1)*dec2mat(lmisys,xopt,X(ibasis));
            end
            Zopt(:,:,k1) = inv(Xk1)/gam^2;
            %F(:,:,k1) = -( C2(:,:,k1)+gam^2*B2(:,:,k1)'*Zopt(:,:,k1) );

            % Explicit solution for controller reconstruction
            p12t = C1(:,:,k1)*Xk1 + D11(:,:,k1)*Bhat(:,:,k1)'/gam^2;
            p22 = -eye(ne-nu) + D11(:,:,k1)*D11(:,:,k1)'/gam^2;
            p23t = D12(:,:,k1)*D11(:,:,k1)'/gam^2;
            F(:,:,k1) = -( B2(:,:,k1)'+ D12(:,:,k1)*Bhat(:,:,k1)'/gam^2 ...
                -p23t*(p22\p12t) )*Zopt(:,:,k1)*gam^2 - C2(:,:,k1);

            % Undo orthog transformation
            F(:,:,k1) = r2inv(:,:,k1)*F(:,:,k1);
        end
    end
end
   




