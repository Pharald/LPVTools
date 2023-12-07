function [K,Gamma,Info] = lpvsyn(sys,nmeas,ncont,Xb,Yb,opt)
% LPVSYN  Parameter-dependent controller synthesis for PSS
%
% [K,GAM,INFO] = LPVSYN(P,NMEAS,NCON) computes a parameter-varying
% controller K which minimizes the induced L2 norm of the interconnection 
% defined by lft(P,K). K is a PSS with NMEAS inputs and NCON outputs, 
% defined on same domain as P. GAM is the induced L2 norm of lft(P,K).
% This three argument call assumes that the rate-bounds of the independent
% variables in P are [-inf,inf]. INFO is a structure containing data from
% the Linear Matrix Inequalities that are solved to obtain K.
%
% [K,GAM,INFO] = LPVSYN(P,NMEAS,NCON,Xb,Yb) computes the rate-bounded 
% parameter-varying controller K for a system P. K is the controller which 
% minimizes the induced L2 norm of lft(P,K) when the rate-bounds of the  
% independent variables of P are incorporated into the synthesis. 
% Xb and Yb are BASIS objects, which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for K.
%
% [K,GAM,INFO] = LPVSYN(P,NMEAS,NCON,Xb,Yb,OPT) allows the user to pass in
% a LPVSYNOPTIONS object. 
%
% The default algorithm for LPVSYN will solve the given synthesis problem
% twice. The first iteration attempts to find a solution that minimizes the
% induced L2 norm of lft(P,K). The second iteration will solve the 
% optimization problem again, with the caveat that any solution that is 
% found to lie within 15% of the optimal induced L2 norm of lft(P,K) from 
% the first iteration, is satisfactory. This formulation has been found to 
% yield controllers that are better numerically conditioned. The back-off 
% factor of 15% can be reset to a different value in LPVSYNOPTIONS.
%
% See also: lpvsynOptions, lpvsfsyn, lpvestsyn, lpvncfyn, lpvmixsyn, lpvloopshape.

% TODO PJS 5/17/2011:  Combine rate bounds into a single niv-by-3 cell array?
% RateUB, (<=niv)-by-2, cell array, {PMAT/DOUBLE variableName}
% RateLB, (<=niv)-by-2, cell array, {PMAT/DOUBLE variableName}

% Parse Inputs
nin = nargin;
narginchk(3, 8)
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

if isempty(Yb)
    Yb = basis(1,0);
elseif ~isa(Yb,'basis')
    error('Yb must be a BASIS object')
end
if isempty(Xb)
    Xb = basis(1,0);
elseif ~isa(Xb,'basis')
    error('Xb must be a BASIS object')
end

% TODO - The LMI conditions can be written to handle parameter dependent
% rate bounds but gridded object do not currently support this.

% Single balancing transformation
% TODO PJS 5/20/2011: Do this before or after orthogonalization?
nd = size(sys,2) - ncont;
ne = size(sys,1) - nmeas;
blk = [nd ne; ncont nmeas];

if contains(opt.balance,'off')
    sysb = sys;
else
sysb = lpvbalance(sys,blk);
end

% Orthogonalize general OLIC interconnection
[P,TL,TR,FT] = orthog4syn(sysb,nmeas,ncont);

% Transform system and basis to non-LPVTools format.
[Pdata,RBx,BFx,Px,RBy,BFy,Py] = basis2data(P,Xb,Yb);

% Get state space data of PMAT system P
[a,b,c,d] = ssdata(P);
nX = size(a,1);
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


% Dimensions
Pdata = Pdata(:,:,:);
szP = size(Pdata);
if numel(szP) ==2
    szP = [szP 1];
end
nmod = szP(3); % # of model in the Pdata array.
nparx = size(RBx,1); % # of parameters in the basis functions
nbasisx = size(BFx,1); % # of basis functions
npary = size(RBy,1); % # of parameters in the basis functions
nbasisy = size(BFy,1); % # of basis functions
Xbivn = Xb.IVName;
Ybivn = Yb.IVName;

% Check for non-rate bounded case
RateBndFlag = 1;
if (nbasisx==1 && nparx==0) && (nbasisy==1 && npary==0)
    RateBndFlag = 0; 
else
    % User has input non-constant basis functions, intends Rate-bounded syn
    % Make sure that all system parameters used in basis functios have 
    % finite ratebounds.
    RBmat = [RBx;RBy];
    RBcheck = isfinite(RBmat);
    RBresult = all(RBcheck(:));
    if ~RBresult
        error(['While attempting rate-bounded synthesis for system P. '...
               'BASIS functions (Xb or Yb) depend on a parameter that '...
               'has non-finite rate-bounds in P.'])
    end    
end



% Create LMI variables
setlmis([]);
Xdec = zeros(nbasisx,1);
for k=1:nbasisx
    Xdec(k) = lmivar(1,[nX 1]);
end
Ydec = zeros(nbasisy,1);
for k=1:nbasisy
    Ydec(k) = lmivar(1,[nX 1]);
end

if isequal(Method,'PoleCon')
    % Hold fixed value of Gamma
    Gamma = opt.Gammaub;
    ginv = 1/Gamma;
    
    %     [ginv,ndec] = lmivar(1,[1 1]);
else
    % Variable representing 1/Gamma
    [ginv,ndec] = lmivar(1,[1 1]);
    
    if isequal(Method,'MaxFeas');
        % LBC is lower bound on X in the coupling constraint
        [LBC,ndec] = lmivar(1,[nX 0]);
    end
end

cnt = 1;
for i=1:nmod
    % Grab state-space data for i-th model
    AhatV = Ahat.Data(:,:,i);
    AtilV = Atil.Data(:,:,i);
    BhatV = Bhat.Data(:,:,i);
    CtilV = Ctil.Data(:,:,i);
    b2V = b2.Data(:,:,i);
    c11V = c11.Data(:,:,i);
    d111dotV = d111dot.Data(:,:,i);
    d11dot1V = d11dot1.Data(:,:,i);
    c2V = c2.Data(:,:,i);
    b11V = b11.Data(:,:,i);
    
    % Grab basis function and partial at i-th point
    if RateBndFlag
        Pxi = Px(:,:,i);
    end
    BFxi = BFx(:,:,i);
    if RateBndFlag
        Pyi = Py(:,:,i);
    end
    BFyi = BFy(:,:,i);
    
    % X Ric LMI
    for q = 1:2^nparx  % number of variables appearing in X's basis fcns
        RBV = RBx(:,1);
        tmp = dec2bin(q-1,nparx);
        idx = find(tmp=='1');
        RBV(idx) = RBx(idx,2);
        for k=1:nbasisx
            for j=1:nparx
                % NOTE PJS 5/15/2011: This term appears to be non-symm
                % because there is a left factor but no right factor.
                % However, the left factor is a scalar so this really is
                % a symmetric term.  This is probably ok but I should just
                % confirm that LMILab handles this correctly.
                lmiterm([-cnt 1 1 Xdec(k)],RBV(j)*Pxi(k,j),1);
            end
            lmiterm([cnt 1 1 Xdec(k)],BFxi(k)*AhatV,1,'s');
            lmiterm([cnt 1 2 Xdec(k)],1,BFxi(k)*c11V');
        end
        lmiterm([cnt 1 1 0],-b2V*b2V');
        lmiterm([cnt 2 2 0],-eye(ne1));
        lmiterm([cnt 3 3 0],-eye(nd));
        if isequal(Method,'PoleCon')
            lmiterm([cnt 1 3 0],ginv*BhatV);
            lmiterm([cnt 2 3 0],ginv*d111dotV);
        else
            lmiterm([cnt 1 3 ginv],1,BhatV);
            lmiterm([cnt 2 3 ginv],1,d111dotV);
        end
        cnt = cnt+1;
    end
    
    % Y Ric LMI
    for q = 1:2^npary  % number of variables appearing in Y's basis fcns
        RBV = RBy(:,1);
        tmp = dec2bin(q-1,npary);
        idx = find(tmp=='1');
        RBV(idx) = RBy(idx,2);
        if nd1==0
            % Handle case for LPVMIXSYN with no input distrubance
            % XXX: Also handle the special case above for ne1==0?
            for k=1:nbasisy
                for j=1:npary
                    lmiterm([cnt 1 1 Ydec(k)],RBV(j)*Pyi(k,j),1);
                end
                lmiterm([cnt 1 1 Ydec(k)],1,BFyi(k)*AtilV,'s');
            end
            lmiterm([cnt 1 1 0],-c2V'*c2V);
            lmiterm([cnt 2 2 0],-eye(ne));
            if isequal(Method,'PoleCon')
                lmiterm([cnt 1 2 0],ginv*CtilV');
            else
                lmiterm([cnt 1 2 ginv],1,CtilV');
            end            
        else
            % Standard Case
            for k=1:nbasisy
                for j=1:npary
                    lmiterm([cnt 1 1 Ydec(k)],RBV(j)*Pyi(k,j),1);
                end
                lmiterm([cnt 1 1 Ydec(k)],1,BFyi(k)*AtilV,'s');
                lmiterm([cnt 1 2 Ydec(k)],1,BFyi(k)*b11V);
            end
            lmiterm([cnt 1 1 0],-c2V'*c2V);
            lmiterm([cnt 2 2 0],-eye(nd1));
            lmiterm([cnt 3 3 0],-eye(ne));
            if isequal(Method,'PoleCon')
                lmiterm([cnt 1 3 0],ginv*CtilV');
                lmiterm([cnt 2 3 0],ginv*d11dot1V');
            else
                lmiterm([cnt 1 3 ginv],1,CtilV');
                lmiterm([cnt 2 3 ginv],1,d11dot1V');
            end
        end
        cnt = cnt+1;
    end
    
    % TODO PJS 5/16/2011--Revisit this. Constraints on X/Y do not need
    % to be repeated at each grid point for non-rate bounded problems.
    % There is probably a more general way to avoid repeated constraints.
    if RateBndFlag || (i==1)
        % X/Y Coupling
        if isequal(Method,'PoleCon')
            lmiterm([-cnt 1 2 0],ginv*eye(nX));
            %             lmiterm([-cnt 1 2 ginv],1,eye(nX));
        else
            lmiterm([-cnt 1 2 ginv],1,eye(nX));
            if isequal(Method,'MaxFeas')
                lmiterm([cnt 1 1 LBC],1,1);
            end
        end
        
        % XXX - These LMIs are redundant if Xlb >0 or Ylb>0
        for k=1:nbasisx
            lmiterm([-cnt 1 1 Xdec(k)],BFxi(k),1);
        end
        for k=1:nbasisy
            lmiterm([-cnt 2 2 Ydec(k)],BFyi(k),1);
        end
        cnt = cnt + 1;
        
        % Xlb*I < X < Xub*I
        if opt.Xlb>0
            lmiterm([cnt 1 1 0],opt.Xlb*eye(nX));
            for k=1:nbasisx
                lmiterm([-cnt 1 1 Xdec(k)],BFxi(k),1);
            end
            cnt = cnt+1;
        end
        if isfinite(opt.Xub)
            lmiterm([-cnt 1 1 0],opt.Xub*eye(nX));
            for k=1:nbasisx
                lmiterm([cnt 1 1 Xdec(k)],BFxi(k),1);
            end
            cnt = cnt+1;
        end
        
        % Ylb*I < Y < Yub*I
        if opt.Ylb>0
            lmiterm([cnt 1 1 0],opt.Ylb*eye(nX));
            for k=1:nbasisy
                lmiterm([-cnt 1 1 Ydec(k)],BFyi(k),1);
            end
            cnt = cnt+1;
        end
        if isfinite(opt.Yub)
            lmiterm([-cnt 1 1 0],opt.Yub*eye(nX));
            for k=1:nbasisy
                lmiterm([cnt 1 1 Ydec(k)],BFyi(k),1);
            end
            cnt = cnt+1;
        end
    end
end

if isequal(Method,'PoleCon')
    % Pole Constraint: eta*[X ginv; ginv Y] > -0.5*(PHI+PHI')
    for i=1:nmod
        % Grab state-space data for i-th model
        AhatV = Ahat.Data(:,:,i);
        AtilV = Atil.Data(:,:,i);
        b11V = b11.Data(:,:,i);
        b12V = b12.Data(:,:,i);
        b2V = b2.Data(:,:,i);
        c11V = c11.Data(:,:,i);
        c12V = c12.Data(:,:,i);
        c2V = c2.Data(:,:,i);
        
        % Grab basis function and partial at i-th point
        BFxi = BFx(:,:,i);
        BFyi = BFy(:,:,i);
        
        lmiterm([-cnt 1 2 0],eye(nX)/Gamma);
        for k=1:nbasisx
            lmiterm([-cnt 1 1 Xdec(k)],BFxi(k),1);
            
            lmiterm([cnt 1 1 Xdec(k)],-0.5*BFxi(k)*AhatV,1,'s');
            lmiterm([cnt 2 1 Xdec(k)],0.5*(c11V'*c11V)/Gamma,1);
        end
        lmiterm([cnt 1 1 0],b2V*b2V');
        
        for k=1:nbasisy
            lmiterm([-cnt 2 2 Ydec(k)],BFyi(k),1);
            
            lmiterm([cnt 2 2 Ydec(k)],1,-0.5*AtilV*BFyi(k),'s');
            lmiterm([cnt 2 1 Ydec(k)],1,0.5*(b11V*b11V')/Gamma);
        end
        lmiterm([cnt 2 2 0],c2V'*c2V);
        lmiterm([cnt 1 2 0],-0.5*(b2V*c12V+b12V*c2V)/Gamma);
        
        cnt = cnt+1;
    end
else
    % 1/Gammaub <= 1/Gamma <= 1/Gammalb
    %%gmin = 0.9;	gmax = 500;	gtol = 1;
    %gmin = .1; gmax = 10000; gtol = 1;
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
end

% TODO PJS 5/25/2011: Add lb/ub constraints on LBC? These constraints are
% in the old code

if ~isequal(Method,'PoleCon')
    % Create objective function:
    % Method = 'MaxFeas':   min -LOB <---> max LOB
    % Method = 'MinGamma' or 'BackOff':  min -1/Gamma <---> min Gamma
    cobj = zeros(ndec,1);
    cobj(end) = -1;
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
    Gamma = inf;
    Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys);
    return;
end

% Get optimal X/Y/Gamma variables
if ~isequal(Method,'PoleCon')
    Gamma = 1/dec2mat(lmisys,xopt,ginv);
end

gamsq = Gamma^2;
gamsqi = 1/gamsq;

X = LOCALdec2pmat(nX,lmisys,xopt,Xdec,P.Domain,BFx);
Y = LOCALdec2pmat(nX,lmisys,xopt,Ydec,P.Domain,BFy);

if isequal(Method,'BackOff');
    % Solve Stage 2 LMI: Maximize the Min Eig of X/Y Coupling Constraint
    opt2 = opt;
    opt2.Method = 'MaxFeas';
    opt2.Gammaub = opt.BackOffFactor*Gamma;
    
    % Construct feasible solution from optimal Stage 1 LMI answer
    x0 = [xopt; 0];
    opt2.SolverInit = x0;
    
    % Solve Relaxed LMI
    Info1 = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,'X',X,'Y',Y);
    Gamma1 = Gamma;
    [K,Gamma,Info2] = lpvsyn(sys,nmeas,ncont,Xb,Yb,opt2);
    Info = struct('MinGamma',Gamma1,'Stage1Info',Info1,'Stage2Info',Info2);
else
    % Reconstruct Controller
    % TODO PJS 5/20/2011: Put this in a try/catch?
    
    % Compute d(X^-1)/dt = -inv(X)*[ sum_i drho_i/dt* dX/drho_i ]*inv(X)
    Xdot = zeros(nX,nX);
    for j=1:nparx
        % XXX This needs to be fixed if we do parameter depended ratebounds
        newVar = pgrid([Xbivn{j} 'Dot'],RBx(j,:) );
        tmpMat = LOCALdec2pmat(nX,lmisys,xopt,Xdec,P.Domain,Px(:,j,:));
%         tmpMat = zeros(nX,nX);
%         for k=1:nbasisx
%             tmpMat = tmpMat + Pxpmat(k,j)*dec2mat(lmisys,xopt,Xdec(k));
%         end
        Xdot = Xdot + newVar*tmpMat;
    end
    xinvpart = -X\(Xdot/X);
    
    %     xinvpart = zeros(nX,nX);
    %     for j=1:nparx
    %         % XXX This needs to be fixed if we do parameter depended ratebounds
    %         newVar = idf([Xbivn{j} 'Dot'],double([XivLB{j} XivUB{j}]));
    %         tmpMat = zeros(nX,nX);
    %         for k=1:nbasisx
    %             tmpMat = tmpMat + Xp{k,j}*dec2mat(lmisys,xopt,Xdec(k));
    %         end
    %         xinvpart = xinvpart + newVar*tmpMat;
    %     end
    %     xinvpart = -X\(xinvpart/X);
    
    % Controller Formulas
    omeg = -d1122 - d1121*((gamsq*eye(nd1,nd1) - d1111'*d1111)\d1111')*d1112;
    b2omeg = b2*omeg;
    d12omeg = d12*omeg;
    
    abar = a + b2omeg*c2;
    b1bar = b1 + b2omeg*d21;
    c1bar = c1 + d12omeg*c2;
    d11bar = d11 + d12omeg*d21;
    
    dh = inv( eye(ne,ne)- gamsqi*(d11bar*d11bar') );
    dt = inv( eye(nd,nd) - gamsqi*(d11bar'*d11bar) );
    
    F = (-d12'*dh*d12)\((b2+b1bar*d11bar'*dh*d12*gamsqi)'/X + d12'*dh*c1bar);
    L = -(Y\(c2+d21*dt*d11bar'*c1bar*gamsqi)' + b1bar*dt*d21')/(d21*dt*d21');
    
    Af = abar + b2*F;
    cf = c1bar + d12*F;
    
    Afx = X\Af;
    left = X\b1bar + cf'*d11bar;
    
    H = -(Afx+Afx'+xinvpart+cf'*cf+gamsqi*left*dt*left');
    
    % Controller State Matrices reconstruction
    % XXX 24Jan14 GJB
    % Issue currently with interpolating the controller generated
    % with the reconstruction technique. Use old controller reconstruction
    % algorithms
    copt = 1;
    if copt==1
        % Shortcourse
        q = Y - gamsqi*inv(X);
        qiy = q\Y;
        
        m1 = H + F'*(b2'/X + d12'*cf);
        %m2 = gamsqi*(gamsq*q*(-qiy*L*d21 - b1bar) + F'*d12'*d11bar)*dt*left';
        m2 = (q*(-qiy*L*d21 - b1bar) + gamsqi*F'*d12'*d11bar)*dt*left';
        m = m1 + m2;
        
        AC = Af + qiy*L*c2 - gamsqi*(q\m);
        BC = -qiy*L;
        CC = F;
        DC = omeg;
    elseif copt==0
        % Lawton Lee's PhD Thesis/Gahinet Reconstruction
        % Compute Ydot (Shortcourse notation)
        
        Ydot = zeros(nX,nX);
        for j=1:npary
            % XXX This needs to be fixed if we do parameter depended ratebounds
            newVar = pgrid([Ybivn{j} 'Dot'],RBy(j,:));
            tmpMat = LOCALdec2pmat(nX,lmisys,xopt,Ydec,P.Domain,Py(:,j,:));
            %     tmpMat = zeros(nX,nX);
            %     for k=1:nbasisy
            %         tmpMat = tmpMat + Pypmat(k,j)*dec2mat(lmisys,xopt,Ydec(k));
            %     end
            Ydot = Ydot + newVar*tmpMat;
        end
        
        %
        ahat = abar + gamsqi*b1bar*d11bar'*dh*c1bar;
        b2hat = b2 + gamsqi*b1bar*d11bar'*dh*d12;
        c2hat = c2 + gamsqi*d21*d11bar'*dh*c1bar;
        
        % Switch to X/Y notation in Lee's thesis
        gami = 1/Gamma;
        Xlee = Y*Gamma;
        Ylee = X*Gamma;
        Xleedot = Ydot*Gamma;
        Yleedot = Xdot*Gamma;
        
        % Factorization of X/Y solution
        %         [N,Q] = schur(Xlee-inv(Ylee));
        %         M = -Ylee*N*Q;
        %
        %         [U,Sig,V] = svd( eye(nX) - Xlee*Ylee );
        %         Srt = Sig.^0.5;
        %         N = U*Srt;
        %         M = V*Srt;
        
        Yrt = sqrtm(Ylee);
        Ynrt = inv(Yrt);
        Z = Yrt*Xlee*Yrt-eye(nX);
        Zrt = sqrtm(Z);
        N = Yrt\Zrt;
        M = -Yrt*Zrt;
        
        % Compute derivatives
        Yrtdot = lyap( Yrt, -Yleedot);
        Ynrtdot = lyap( Ynrt, Ylee\(Yleedot)/Ylee );
        
        tmp = -(Yrtdot*Xlee*Yrt + Yrt*Xleedot*Yrt + Yrt*Xlee*Yrtdot);
        Zrtdot = lyap(Zrt,tmp);
        
        Nleedot = Ynrtdot*Zrt + Ynrt*Zrtdot;
        % Mleedot = -(Yrtdot*Zrt+Yrt*Zrtdot);
        % Comment: Residual should be zero
        % Residual=Xleedot*Ylee+Nleedot*M'+Xlee*Yleedot+N*Mleedot';
        
        % Compute state-matrices
        %         AC = -N\( ahat' + Xlee*(ahat+b2hat*F+L*c2hat)*Ylee ...
        %             + gami*c1bar'*dh*(c1bar+d12*F)*Ylee ...
        %             + gami*Xlee*(b1bar+L*d21)*dt*b1bar' + ...
        %             Xleedot*Ylee + Nleedot*M' )/(M');
        %         BC = N\(Xlee*L);
        %         CC = (F*Ylee)/(M');
        %         DC = omeg;
        
        % Compute Controller state-matrices minimizing X, Y inverses
        % Construct F*X and Y*L to eliminate X, Y inverses
        %F = (-d12'*dh*d12)\((b2+b1bar*d11bar'*dh*d12*gamsqi)'/X + d12'*dh*c1bar);
        FX = (-d12'*dh*d12)\((b2+b1bar*d11bar'*dh*d12*gamsqi)'   + d12'*dh*c1bar*X);
        %L = -(Y\(c2+d21*dt*d11bar'*c1bar*gamsqi)' + b1bar*dt*d21')/(d21*dt*d21');
        YL = -(  (c2+d21*dt*d11bar'*c1bar*gamsqi)' + Y*b1bar*dt*d21')/(d21*dt*d21');
        AC = -N\( ahat' ...
            + Y*ahat*X*Gamma^2 + Y*b2hat*FX*Gamma^2 + YL*c2hat*X*Gamma^2 ...
            + c1bar'*dh*(c1bar*X+d12*FX) ...
            + (Y*b1bar+YL*d21)*dt*b1bar' ...
            + Xleedot*X*Gamma)/(M') - N\Nleedot;
        BC = N\(YL*Gamma);
        CC = (FX*Gamma)/(M');
        DC = omeg;
    end
    
    % Undo orthogonalizing transformation
    ortcont = ss(AC,BC,CC,DC);
    K = lft([TR;eye(ncont)]*ortcont*[TL' eye(nmeas)],FT);
    
    % Form Info Structure
    Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,'X',X,'Y',Y);
end

function Y = LOCALdec2pmat(nX,lmisys,xopt,Ydec,Domain,BFy)
% helper function that transforms the decision variables into PMAT objects.
% Y = zeros(nX,nX);
% for k=1:nbasisy
%     Y = Y + BFypmat(k)*dec2mat(lmisys,xopt,Ydec(k));
% end
% Y = lpvsplit(Y,sys.DomainPrivate); XXX Why this?
        
        

nbasisy = size(BFy,1);
nmod = size(BFy,3);

Ymat = zeros(nX,nX,nbasisy);
for k=1:nbasisy
    Ymat(:,:,k) = dec2mat(lmisys,xopt,Ydec(k));
end
    
Yopt = zeros(nX,nX,nmod);
for i1 =1:nmod
    for i2 = 1:nbasisy
        Yopt(:,:,i1) = Yopt(:,:,i1)+BFy(i2,1,i1)*Ymat(:,:,i2);
    end
end
Yopt = reshape(Yopt,[nX;nX;Domain.LIVData]');
if isempty(Yopt)
    Y = pmat([]);
else
    Y = pmat(Yopt,Domain);
end



