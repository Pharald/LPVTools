% function [GAM, lamopt, Xopt] = iqcgain(G,Fbasis,Fgrad,RateBounds,Ublksize,IQCS)
%
% Compute an upper bound GAM on worst-case gain of Fu(G,Delta) where
% the input/output behavior of the uncertainty is described by IQCs.
%
% G: Array of LTI systems.
% Fbasis: Basis functions used in storage function Nbasis-by-1-by-nmod
% Fgrad: Gradient of the basis functions wrt. the parameters:
%        Nbasis-by-Nparameters-by-nmod
% RateBounds: Upper and lower bounds on the parameter rates: Nparameter-by-2
% Ublksize: Nblk-by-2 matrix of positive integers. Each row is [R,C] where 
%        R and C are the number of rows and columns of the uncertainty.
% IQCS:  Nblk-by-1 structured array specifying the IQCs in factorized
%        form as (psi,M).  The fields of IQCS are:
%        IQCS(i).psi: Ni-by-1 cell array of LTI systems
%        IQCS(i).M: Ni-by-1 cell array of real, symmetric matrices.
%
%
% XXX Currently no error checking. The file assumes that all dimensions
% are compatible.

function [gam,lamopt,Xopt] = iqcgridengine(G,Fbasis,Fgrad,RateBounds,Ublksize,iqcs)


% Dimensions
G = G(:,:,:);
szG = size(G);
if numel(szG) ==2
    szG = [szG 1];
end
nve = szG(1);  % # of outputs of G = nv + ne
nzd = szG(2);  % # of inputs of G = nz + nd
nmod = szG(3); % # of model in the G array.
if isempty(Ublksize)
    Nblk = 0;
    nz=0;
    nv=0;
else
    Nblk = size(Ublksize,1);
    szU = sum(Ublksize,1);
    nz = szU(1);   % # of outputs of Delta
    nv = szU(2);   % # of inputs to Delta
end 
nd = nzd-nz;   % # of disturbances
ne = nve-nv;   % # of errors
npar = size(RateBounds,1); % # of parameters
nbasis = size(Fbasis,1); % # of basis functions


% Check for non-rate bounded case
ratebndflg = true;
if (nbasis==1 && npar==0) 
    ratebndflg = false;
end

% Construct extended system Pext including dynamics of plant G and
% all auxiliary systems Psi. The I/O relationship for Pext is:
%     [y1; ...; yNiqc; e] = Pext*[z1; ...; zNblk; d]
% where zi are the outputs of uncertainty i and yj are the outputs
% of the auxiliary system psi for the j^th IQC.
vidx = cell(Nblk,1);
zidx = cell(Nblk,1);

zptr = 0;
vptr = 0;
AllPsi = [];
AllM = [];
yidx = [];
yptr = 0;
GIidx = [];
Niqcs = 0;
for i1=1:Nblk
    % Grab indices associated with i^th uncertaint I/O channels

    zidx{i1} = [zptr+1:zptr+Ublksize(i1,1)]';
    zptr = zptr + Ublksize(i1,1);
    vidx{i1} = [vptr+1:vptr+Ublksize(i1,2)]';
    vptr = vptr + Ublksize(i1,2);
    
    % Grab IQCs for the i1^th uncertainty and stack together
    psi = iqcs(i1).psi(:);
    AllPsi = blkdiag(AllPsi, vertcat(psi{:}));    
    
    % Store indices of v/z associated with I/O of i1^th uncertainty
    GIidx = [GIidx; vidx{i1}; nve+zidx{i1} ];
    
    % Store output indices associated with all iqcs
    Ni = numel(psi);
    Niqcs = Niqcs+Ni;    
    AllM = [AllM; iqcs(i1).M(:)];
    for i2=1:Ni
        szyi = size(psi{i2},1);
        yidx = [yidx; {yptr + (1:szyi)}];        
        yptr = yptr+szyi;
    end    
end
ny = yptr; % # of y dimensions

% Append indices for error output channels
GIidx = [GIidx; [nve-ne+1:nve]'];
AllPsi = blkdiag(AllPsi,ss(eye(ne)));
AllPsi = repsys(AllPsi,[1 1 nmod]);
GI = [G;repmat(eye(nzd),[1 1 nmod])];
GI = GI(GIidx,:,:);
Pext = AllPsi*GI;

% Partition state matrices of Pext based on inputs [z;d] and outputs [y;e]
[A,B,C,D] = ssdata(Pext);
nx = size(A,1);
Ce = C(end-ne+1:end,:,:);
De = D(end-ne+1:end,:,:);

% Initialize LMI and define variables
setlmis([])
X = zeros(nbasis,1);
for i1 = 1:nbasis
    X(i1) = lmivar(1,[nx 1]);    
end
for i2 = 1:Niqcs
    lam(i2) = lmivar(1,[1 1]);
end
[gsq,ndec] = lmivar(1,[1 1]);
gsqI = lmivar(3,blkdiag(zeros(nx+nz),ndec*eye(nd)));

% LMI for Dissipation Inequality - Loop through array of models
cnt = 1;
for k1 = 1:nmod
    Ak = A(:,:,k1);
    Bk = B(:,:,k1);
    Cek = Ce(:,:,k1);
    Dek = De(:,:,k1);
    
    if ratebndflg
        Fgk1 = Fgrad(:,:,k1);        
    else
        Fgk1 = [];        
    end
    Fbk1 = Fbasis(:,:,k1);
            
    % L2 Bound LMI in block 3-by-3 form with gamma
    for q = 1:2^npar  % number of variables appearing in X's basis fcns
        
        for ibasis=1:nbasis
            lmiterm([cnt 1 1 X(ibasis)],[eye(nx); zeros(nz+nd,nx)],[Ak Bk]*Fbk1(ibasis,1),'s');
        end
        lmiterm([-cnt 1 1 gsqI],1,1);
        lmiterm([cnt 1 1 0],[Cek Dek]'*[Cek Dek]);
        for i1=1:Niqcs
            % Grab i1^th IQC
            Mi = AllM{i1};
            Cyi = C(yidx{i1},:,k1);
            Dyi = D(yidx{i1},:,k1);
            Mtil = [Cyi Dyi]'*Mi*[Cyi Dyi];
            lmiterm([cnt 1 1 lam(i1)],Mtil,1,'s');
        end        
        
        if ratebndflg
            rbvec = RateBounds(:,1);
            tmp = dec2bin(q-1,npar);
            idx = find(tmp=='1');
            rbvec(idx) = RateBounds(idx,2);              
            for k2=1:nbasis
                for k3=1:npar
                    lmiterm([-cnt 1 1 X(k2)],[eye(nx); zeros(nzd,nx)],...
                            0.5*rbvec(k3)*Fgk1(k2,k3)*...
                            [eye(nx) zeros(nx,nzd)],'s');
                end
            end
        end
        cnt = cnt+1;
    end
    if ratebndflg || k1 == 1
        % xpdlow*I < X < xpdupp*I
%         xpdlow = 1e-6;
%         lmiterm([cnt 1 1 0],xpdlow*eye(nx));
        for n1=1:nbasis
            lmiterm([-cnt 1 1 X(n1)],Fbk1(n1),1);
        end
        cnt = cnt+1;
        
        % TODO PJS 10/22/2011: What value should we choose here?
%         xpdupp = 1e6;
%         %xpdupp = 1e9;
%         lmiterm([-cnt 1 1 0],xpdupp*eye(nx));
%         for n2=1:nbasis
%             lmiterm([cnt 1 1 X(n2)],Fbk1(n2),1);
%         end
%         cnt = cnt+1;
    end
end

% LMI to enforce scalings to be non-negative
for i1 = 1:Niqcs
    lmiterm([-(cnt+i1-1) 1 1 lam(i1)],1,1);
    cnt = cnt+1;
end

% SDP: min gamsq subject to LMI constraints
lmisys = getlmis;
c = zeros(ndec,1);
c(end) = 1;

opt = [0 0 0 0 1];
[copt,xopt] = mincx(lmisys,c,opt);
gam = sqrt(copt);
if ~isempty(xopt)
    for i1 = 1:nbasis
        Xopt(:,:,i1) =  dec2mat(lmisys,xopt,X(i1));
    end
    for i2 = 1:Niqcs
        lamopt(i2) = dec2mat(lmisys,xopt,lam(i2));
    end
else
    Xopt = [];
    lamopt = [];
end



