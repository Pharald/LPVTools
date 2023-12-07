% function [F,GAM,INFO] = lpvLQGsfsynengine(G,ncont,Fbasis,Fgrad,RateBounds)
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
%
% OUTPUTS
% F: nu-by-nx-nmod array of state feedback gains
% GAM: Min Hinf cost


% XXX Currently no error checking. The file assumes that all dimensions
% are compatible.
% XXX Add options object to pass through solver options

function [F,gam,Info] = lpvLQGsfsynengine(G,nu,Fbasis,Fgrad,RateBounds)


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
if any(D1(:))
    error('Feedthrough matrix (D11) from disturbance to error must be zero')
end
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

C1 = zeros(ne-nu,nx,nmod);
C2 = zeros(nu,nx,nmod);
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
    D2(:,:,i) = [zeros(ne-nu,nu); eye(nu)];
    
    C1(:,:,i) = C(1:ne-nu,:,i);
    C2(:,:,i) = C(ne-nu+1:end,:,i);
    
    Ahat(:,:,i) = A(:,:,i)-B2(:,:,i)*C2(:,:,i);
end

% Initialize LMI and define variables
setlmis([])
X = zeros(nbasis,1);
for i1 = 1:nbasis
    X(i1) = lmivar(1,[nx 1]);
end
[Z,ndec,Zdec] = lmivar(1,[nd 1]);

% LMI for Dissipation Inequality - Loop through array of models
cnt = 1;
for k1 = 1:nmod
    Ahatk = Ahat(:,:,k1);
    B1k = B1(:,:,k1);
    B2k = B2(:,:,k1);
    C1k = C1(:,:,k1);
    
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
        xpdlow = 1e-6;
        lmiterm([cnt 1 1 0],xpdlow*eye(nx));
        lmiterm([cnt 2 2 0],xpdlow*eye(nd));
        for k=1:nbasis
            lmiterm([-cnt 1 1 X(k)],Fbk1(k),1);
        end
        lmiterm([-cnt 1 2 0],B1k);
        lmiterm([-cnt 2 2 Z],1,1);
        cnt = cnt+1;
    end
end

% Create objective function: min gam
cobj = zeros(ndec,1);
cobj(diag(Zdec)) = 1;

% Solve LMI
lmisys = getlmis;
opt = [0 0 0 0 1];
[copt,xopt] = mincx(lmisys,cobj,opt);
gam = sqrt(copt);
info = [];
if ~isempty(xopt)
    % Controller reconstruction
    F = zeros(nu,nx,nmod);
    Popt = zeros(nx,nx,nmod);
    for k1=1:nmod
        Fbk1 = Fbasis(:,:,k1);
        Xk1 = zeros(nx);
        for ibasis=1:nbasis
            Xk1 = Xk1 + Fbk1(ibasis,1)*dec2mat(lmisys,xopt,X(ibasis));
        end
        Popt(:,:,k1) = inv(Xk1)/gam^2;
        
        % Explicit solution for controller reconstruction
        p12t = C1(:,:,k1)*Xk1;
        p22 = -eye(ne-nu);
        F(:,:,k1) = -B2(:,:,k1)'*Popt(:,:,k1)*gam^2 - C2(:,:,k1);
        
        % Undo orthog transformation
        F(:,:,k1) = r2inv(:,:,k1)*F(:,:,k1);
    end
    
    
    Zopt = dec2mat(lmisys,xopt,Z);
    % Build X optimal
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
    
else
    F = [];
    Popt = [];
end
Info = struct('xopt',xopt,'copt',copt,'lmisys',lmisys,'Zopt',Zopt,'Xopt',Xopt);


