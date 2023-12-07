function [dMd,dr,dci] = lpvosbal(M,blk)
% M is Nr-by-Nc-Npts
% blk is nblk-by-2 of block dimensions (no repeated scalar blocks allowed)

% Get sizes
szM = [size(M) 1];
Nr = szM(1);
Nc = szM(2);
Npts = prod(szM(3:end));
nblk = size(blk,1);
if any( sum(blk) ~= [szM(2) szM(1)] )
    error('Dimensions of M are incompatible with the dimensions of blk.')
end
nfull = sum( blk(:,1)>1 | blk(:,2)>1 );

% Get block indices
cptr = 1;
rptr = 1;
cnt = 0;
fullrows = cell(nfull,1);
fullcols = cell(nfull,1);
afullidx = zeros(nfull,1);
sr = zeros(nblk,Nr);
sc = zeros(Nc,nblk);
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

% Compute norms of blocks of M
A = zeros(nblk,nblk,Npts);
for k=1:Npts
    Mk = M(:,:,k);
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
A = max(A,[],3);
A = A-diag(diag(A));

% Perform conventional Osborne's method
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
d = sqrt(d/d(1,1));
drvec = diag(sr'*diag(d)*sr);
dcivec = diag(sc*diag(1./d)*sc')';
dr = diag(drvec);
dci = diag(dcivec);

% Apply scalings to all input matrices
dMd = zeros(szM);
for i0 = 1:Npts
    dMd(:,:,i0) = M(:,:,i0).*(drvec*dcivec);
end
