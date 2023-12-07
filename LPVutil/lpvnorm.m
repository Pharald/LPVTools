function [gam,X,Info] = lpvnorm(P,Xb,alg)
% LPVNORM  Compute bound on norm for PSS systems.
%
% LPVNORM computes a bound on the norm of a PSS system over the set of
% all permissible trajectories of the independend variables which the PSS 
% depends on. 
%
% [Gamma,X] = lpvnorm(P,'L2') computes an upper bound Gamma on the induced 
% L2 norm of the PSS P. The upper bound Gamma and a constant (parameter
% independent) matrix X are computed to satisfy the induced norm linear
% matrix inequality (LMI) condition. X is returned as a PMAT. The L2 norm
% bound is valid for arbitrarily fast variations of the system parameters.
%
% [Gamma,X] = lpvnorm(P,'LQG') computes an upper bound Gamma on the 
% stochastic LPV bound. The stochastic LPV bound is defined as the expected 
% value of the average instantaneous power of the output of P, assuming its
% inputs are zero mean, white-noise processes with unit intensity.
%
% [Gamma,X] = lpvnorm(P,Xb,ALG) computes a tighter (less conservative) 
% bound on the norm by using a parameter dependent matrix X(p) 
% and bounds on the parameter rates of variation. The basis functions used 
% to construct X(p) are specified with the BASIS object Xb. ALG can be 
% either 'L2' or 'LQG'. A call without the ALG argument is equivalent to 
% [Gamma,X] = lpvnorm(P,Xb,'L2').
%
% See also: norm, lpvwcgain, lpvsyn.


% Enforce constraints LMI constraints on domain of P. Allow Basis
% Functions to be specified on another dense grid. This would allow
% partials to be numerically computed.  We would need to use LPVSPLIT
% to evaluate state matrices of P on domain grid points and LPVINTERP
% to evaluate the Basis Functions (and partials, if specified) on these
% domain grid points.


% Parse Inputs
nin = nargin;
narginchk(1, 3)
if nin==1
    Xb = [];
    alg = 'L2';
elseif nin==2
    if isa(Xb,'basis')
        alg = 'L2';
    elseif isa(Xb,'char')
        alg = Xb;
        Xb = [];
    else
        error(['The second argment in a call to LPVNORM must either be '...
            'a BASIS object or a CHAR'])
    end
end

% Make sure user passed in a legit algorithm choice.
if ~(strcmpi(alg,'L2') || strcmpi(alg,'LQG'))
    error('The algorithm ALG must either be set to ''L2" or ''LQG''')
end

% Always assume one basis function: constant
if isempty(Xb)
    Xb = basis(1,0);
end


% Single balancing transformation
% TODO PJS 5/20/2011: Do this balancing?
if 1%strcmpi(alg,'L2')
    G = lpvbalance(P);
else
    % XXX 11/29/14 Disovered poor numerics in LQG norm computations 
    % post-balancing. Look into this.
    G = P;
end

% Map the input data into LPVNORM engine data form:
[Gdata,RBx,BFdatax,Pdatax] = basis2data(G,Xb);

if strcmpi(alg,'L2')
    % Run the L2 LPVNORM engine:
    [gam,Xopt,Info] = LOCALlpvL2normengine(Gdata,RBx,BFdatax,Pdatax);
elseif strcmpi(alg,'LQG')
    % Run the Stochastic LPVNORM engine:
    [gam,Xopt,Info] = LOCALlpvLQGnormengine(Gdata,RBx,BFdatax,Pdatax);
end

nx = size(Xopt,1);
Xopt = reshape(Xopt,[nx;nx;P.Domain.LIVData]');
if isempty(Xopt)
    X = pmat([]);
else
    X = pmat(Xopt,P.Domain);
end


function [gam,Xopt,Info] = LOCALlpvL2normengine(G,RateBounds,Fbasis,Fgrad)

% Dimensions
G = G(:,:,:);
szG = size(G);
if numel(szG) ==2
    szG = [szG 1];
end
ny = szG(1);  % # of outputs of G
nu = szG(2);  % # of inputs of G
nmod = szG(3); % # of model in the G array.
npar = size(RateBounds,1); % # of parameters in the basis functions
nbasis = size(Fbasis,1); % # of basis functions

% Check for non-rate bounded case
ratebndflg = true;
if (nbasis==1 && npar==0)
    ratebndflg = false;
end

% Grab the state-space data for G
[A,B,C,D] = ssdata(G);
nx = size(A,1);

% Create LMI variables (X,gam)
setlmis([]);
X = zeros(nbasis,1);
for k=1:nbasis
    X(k) = lmivar(1,[nx 1]);
end
[gam,ndec] = lmivar(1,[1 1]);

% Construct LMI Conditions at each grid point
cnt = 1;
for k1 = 1:nmod
    % Grab State-space data at k1-th point
    Ak = A(:,:,k1);
    Bk = B(:,:,k1);
    Ck = C(:,:,k1);
    Dk = D(:,:,k1);
    
    % Grab basis function and partial at k1-th point
    if ratebndflg
        Fgk1 = Fgrad(:,:,k1);
    end
    Fbk1 = Fbasis(:,:,k1);
    
    % L2 Bound LMI in block 3-by-3 form with gamma
    for q = 1:2^npar  % number of variables appearing in X's basis fcns
        rbvec = RateBounds(:,1);
        tmp = dec2bin(q-1,npar);
        idx = find(tmp=='1');
        rbvec(idx) = RateBounds(idx,2);
        
        for k=1:nbasis
            for j=1:npar
                % XXX - Error. The -cnt in this LMI is wrong, should be +cnt.
                % We didn't notice so far because we always had symmetric
                % rate-bounds.
                lmiterm([cnt 1 1 X(k)],rbvec(j)*Fgk1(k,j),1);
            end
            lmiterm([cnt 1 1 X(k)],1,Fbk1(k)*Ak,'s');
            lmiterm([cnt 1 2 X(k)],1,Fbk1(k)*Bk);
        end
        lmiterm([cnt 1 3 0],Ck');
        lmiterm([-cnt 2 2 gam],eye(nu),eye(nu));
        lmiterm([cnt 2 3 0],Dk');
        lmiterm([-cnt 3 3 gam],eye(ny),eye(ny));
        cnt = cnt+1;
    end
    
    if ratebndflg || (k1==1)
        % xpdlow*I < X < xpdupp*I
        xpdlow = 1e-6;
        lmiterm([cnt 1 1 0],xpdlow*eye(nx));
        for k=1:nbasis
            lmiterm([-cnt 1 1 X(k)],Fbk1(k),1);
        end
        cnt = cnt+1;
        
        % TODO PJS 10/22/2011: What value should we choose here?
        xpdupp = 1e6;
        %xpdupp = 1e9;
        lmiterm([-cnt 1 1 0],xpdupp*eye(nx));
        for k=1:nbasis
            lmiterm([cnt 1 1 X(k)],Fbk1(k),1);
        end
        cnt = cnt+1;
    end
end

% Create objective function: min gam
cobj = zeros(ndec,1);
cobj(end) = 1;

% Solve LMI
lmisys = getlmis;
opts = zeros(5,1);
opts(2) = 200; % Max # of iters
opts(5) = 1;   % Toggle display
[copt,xopt] = mincx(lmisys,cobj,opts);
gam = copt;

% Store solver info
% TODO PJS 5/16/2011: Additional info to be returned?
Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys);

% Get optimal X
if isempty(copt)
    gam = inf;
    Xopt = pmat;
else
    Xmat = zeros(nx,nx,nbasis);
    for k=1:nbasis
        Xmat(:,:,k) = dec2mat(lmisys,xopt,X(k));
    end
    
    Xopt = zeros(nx,nx,nmod);
    for i1 =1:nmod
        for i2 = 1:nbasis
            Xopt(:,:,i1) = Xopt(:,:,i1)+Fbasis(i2,1,i1)*Xmat(:,:,i2);
        end
    end
end

function [gam,Xopt,Info] = LOCALlpvLQGnormengine(G,RateBounds,Fbasis,Fgrad)
% XXX - 11/29/14 We are implementing the form described by Theorem 1.4.1 
% in Fen Wu's thesis (pdf page 26). We can also implement the form 
% described by Theorem 1.4.2 but there is no guidance for which is a better
% form to use.

% Dimensions
G = G(:,:,:);
szG = size(G);
if numel(szG) ==2
    szG = [szG 1];
end
ny = szG(1);  % # of outputs of G
nu = szG(2);  % # of inputs of G
nmod = szG(3); % # of model in the G array.
npar = size(RateBounds,1); % # of parameters in the basis functions
nbasis = size(Fbasis,1); % # of basis functions

% Check for non-rate bounded case
ratebndflg = true;
if (nbasis==1 && npar==0)
    ratebndflg = false;
end

% Grab the state-space data for G
[A,B,C,~] = ssdata(G);
nx = size(A,1);

% Create LMI variables (X,gam)
setlmis([]);
X = zeros(nbasis,1);
for k=1:nbasis
    X(k) = lmivar(1,[nx 1]);
end
[Z,ndec,Zdec] = lmivar(1,[nu 1]);

% Construct LMI Conditions at each grid point
cnt = 1;
for k1 = 1:nmod
    % Grab State-space data at k1-th point
    Ak = A(:,:,k1);
    Bk = B(:,:,k1);
    Ck = C(:,:,k1);
    
    % Grab basis function and partial at k1-th point
    if ratebndflg
        Fgk1 = Fgrad(:,:,k1);
    end
    Fbk1 = Fbasis(:,:,k1);
    
    % L2 Bound LMI in block 3-by-3 form with gamma
    for q = 1:2^npar  % number of variables appearing in X's basis fcns
        rbvec = RateBounds(:,1);
        tmp = dec2bin(q-1,npar);
        idx = find(tmp=='1');
        rbvec(idx) = RateBounds(idx,2);
        
        for k=1:nbasis
            for j=1:npar
                lmiterm([-cnt 1 1 X(k)],rbvec(j)*Fgk1(k,j),1);
            end
            lmiterm([cnt 1 1 X(k)],Fbk1(k)*Ak,1,'s');
            lmiterm([cnt 1 2 X(k)],1,Fbk1(k)*Ck');
        end
        lmiterm([-cnt 2 2 0],eye(ny));
        cnt = cnt+1;
    end
    
    if ratebndflg || (k1==1)
        % [X B;B' Z] > xpdlow*I_(nx+nu)
%         xpdlow = 1e-6;
%         lmiterm([cnt 1 1 0],xpdlow*eye(nx));
%         lmiterm([cnt 2 2 0],xpdlow*eye(nu));
        for k=1:nbasis
            lmiterm([-cnt 1 1 X(k)],Fbk1(k),1);
        end
        lmiterm([-cnt 1 2 0],Bk);
        lmiterm([-cnt 2 2 Z],1,1);
        cnt = cnt+1;
        
        % Assume  we don't need to upper bound this in H2 framework.
        %         % TODO PJS 10/22/2011: What value should we choose here?
        %         xpdupp = 1e6;
        %         %xpdupp = 1e9;
        %         lmiterm([-cnt 1 1 0],xpdupp*eye(nx));
        %         for k=1:nbasis
        %             lmiterm([cnt 1 1 X(k)],Fbk1(k),1);
        %         end
        %         cnt = cnt+1;
        
    end
end

% Create objective function: min gam
cobj = zeros(ndec,1);
cobj(diag(Zdec)) = 1;

% Solve LMI
lmisys = getlmis;
opts = zeros(5,1);
opts(2) = 200; % Max # of iters
opts(5) = 1;   % Toggle display
[copt,xopt] = mincx(lmisys,cobj,opts);
gam = sqrt(copt);



% Get optimal X
if isempty(copt)
    gam = inf;
    Xopt = pmat;
    Zopt = [];
else
    Xmat = zeros(nx,nx,nbasis);
    for k=1:nbasis
        Xmat(:,:,k) = dec2mat(lmisys,xopt,X(k));
    end
    
    Xopt = zeros(nx,nx,nmod);
    for i1 =1:nmod
        for i2 = 1:nbasis
            Xopt(:,:,i1) = Xopt(:,:,i1)+Fbasis(i2,1,i1)*Xmat(:,:,i2);
        end
    end
    Zopt = dec2mat(lmisys,xopt,Z);
end

% Store solver info
% TODO PJS 5/16/2011: Additional info to be returned?
Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,'Zopt',Zopt);



