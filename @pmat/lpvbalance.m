function [dMd,dr,dci] = lpvbalance(M,blk)
% LPVBALANCE   Diagonal scaling for PMAT objects
%
% [B,D] = lpvbalance(A) computes a single diagonal similarity transformation
% to improve the conditioning of the N-by-N PMAT A at all points in the
% domain.  The transformation D is computed using a generalized version of
% Osborne's iteration. D is returned as an N-by-N diagonal, DOUBLE matrix.
% The scaled PMAT B is D*A*inv(D).  This applies the transformation D at
% each point in the domain of A.
%
% [B,DR,DC] = lpvbalance(A,BLK) computes structured, diagonal transformations
% for the N-by-M PMAT A. The scaled PMAT B is DR*A*DC where DR is an
% N-by-N matrix and and DC is an M-by-M matrix. BLK is a K-by-2 matrix
% that specifies the block partitioning dimensions of DR and DC. If
% BLK = [c1 r1; ... ; ck rk] then DR and DC are partitioned as:
%    DR = blkdiag( d1*I_r1, ...,  dk*I_rk )
%    DC = blkdiag( (1/d1)*I_c1, ...,  (1/dk)*I_ck )
% where the notation I_r1= eye(r1), represents the r1-by-r1 identity matrix.
% The block partitioning must be consistent with the dimensions of A, i.e.
% the sum across the rows of BLK should equal [M N].
%
% See also: balance.

% Parse Inputs
nin = nargin;
narginchk(1, 2)
if nin==1
    blk = ones( size(M,1) , 2 );
end

% Get PMAT and block data
Data = M.DataPrivate;
Dom  = M.DomainPrivate;
szM  = [size(Data) 1];
Npts = prod(szM(3:end));
nblk = size(blk,1);
if any( sum(blk) ~= [szM(2) szM(1)] )
    error('Dimensions of M are incompatible with the dimensions of blk.')
end
if ~all( blk(:)>0 & floor(blk(:))==ceil(blk(:)) )
    error('Block dimensions must be non-negative integers.');
end
nfull = sum( blk(:,1)>1 | blk(:,2)>1 );

% Get block indices used in Osborne's iteration
cptr = 1;
rptr = 1;
cnt = 0;
fullrows = cell(nfull,1);
fullcols = cell(nfull,1);
afullidx = zeros(nfull,1);
sr = zeros(nblk, szM(1) );
sc = zeros(szM(2), nblk);
for i1 = 1:nblk
    cdim = blk(i1,1);
    rdim = blk(i1,2);
    idxr = rptr:rptr+rdim-1;
    idxc = cptr:cptr+cdim-1;
    rptr = rptr + rdim;
    cptr = cptr + cdim;
    
    if rdim>1 || cdim>1
        % Non-scalar (full) block
        cnt = cnt+1;
        fullrows{cnt} = idxr;
        fullcols{cnt} = idxc;
        afullidx(cnt)=i1;
    end
    
    % Mask used to compute block norms
    sr(i1,idxr) = ones(1,rdim);
    sc(idxc,i1) = ones(cdim,1);
end

% Compute norms of blocks of M at each point in the domain
A = zeros(nblk,nblk,Npts);
for k=1:Npts
    Mk = Data(:,:,k);
    Ak = sr*real(conj(Mk).*Mk)*sc;
    for i1 = 1:nfull,
        ridx = fullrows{i1};
        afullr = afullidx(i1);
        for i2 = 1:nfull
            cidx = fullcols{i2};
            afullc = afullidx(i2);
            Ak(afullr,afullc) = norm( Mk(ridx,cidx) )^2;
        end
    end
    A(:,:,k) = Ak;
end

% Compute peak norms of blocks across domain points
A = max(A,[],3);
A = A-diag(diag(A));

% Perform conventional Osborne's method on A
maxiter = 30;     % Max # of iterations
reltol = 1e-4;    % Relative stopping tolerance
d = ones(1,nblk); % Scaling matrix stored as a vector
cost = sum(A(:));
oldcost = max([2*cost 10*eps]);
itcnt = 0;
Astart = A;
while (itcnt < maxiter) && reltol*oldcost<(oldcost-cost)
    sa = max(sum(A),10*eps);
    sat = max(sum(A,2),10*eps)';
    d = d.*sqrt(sqrt(sa./sat));
    d = min(max(d,1e-8),1e8);
    A = Astart.*(d'*(1 ./d));
    itcnt = itcnt+1;
    oldcost = max([cost 10*eps]);
    cost = sum(A(:));
    
    % Stop if ill-conditioned
    dcond = max(d)/max(min(d),10*eps);
    if dcond>1e10
        itcnt = maxiter;
    end
end

% Create full scaling matrix for output
d = sqrt(d/d(end,end));
drvec = diag(sr'*diag(d)*sr);
dcivec = diag(sc*diag(1./d)*sc')';
dr = diag(drvec);
dci = diag(dcivec);

% Apply scalings to all input matrices
dMd = zeros(szM);
for i0 = 1:Npts
    dMd(:,:,i0) = Data(:,:,i0).*(drvec*dcivec);
end

% Pack data back into a PMAT
dMd = pmat(dMd,Dom);
