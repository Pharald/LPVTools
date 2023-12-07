function [feas,S,xopt,perfparm] = iqcsolve4(G,blk,IQCinfoS,lmisys,omega,cidx,PerfFlag)
% Solve IQC LMI using cutting plane method.
%
%  It is assumed that the uncertainties are ordered with udyns at the end.
%  The IQC LMI is
%      [G(j wk); Im]'*S(j wk)*[G(j wk);Im] <= -tol*I    k=1,2,...,Nw
%  The IQC LMI constraints are enforced on a coarse grid of frequencies
%  specified by cidx. cidx is a subset of {1,...,Nw}. If the LMI problem
%  is feasible then S(j wk) is formed at all indices (k=1,..,Nw).
%
% XXX Clean up documentation for performance block case
%
% INPUTS
%  G is an n-by-m-Nw array.  G(:,:,i) is the frequency response data
%      at the i^th frequency point.
%  blk is an Nblk-by-2 array of block dimensions in the old mu notation
%      [-r 0] for repeated reals, (r>0)
%      [c 0] for repeated complex (c>0)
%      [r c] for either full complex (r,c>0)
%  P: Nudyn-by-1 cell array.  allPi{i} is an ni-by-1 cell array.
%       allPi{i}{j} is an (n+m)-by-(n+m)-by-Nw array specifying the j^th
%       IQC multiplier for the i^th block.
%    XXXX dimens are wrong for P{i}{j}.  These are not (n+)-by-(n+m)
%    but they are related to the I/O dimensions of the i^th udyn.
%  Pblk: Nudyn-by-2 array of block dimensions
%  omega: Nw-by-1 array of frequencies
%  cidx: Nc-by-1 array of indices specifying the frequencies at which
%       the LMI constraints are enforced.
%
% OUTPUTS
%  feas: feas=1 if the LMI problem is feasible. Otherwise feas=0.
%  S: (n+m)-by-(n+m)-by-Nw array of IQC multipliers. If feas=1 then S will
%    satisfy the IQC LMI constraints at the indices specified by cidx but
%    may not satisfy the constraints at the other indices.
%  xopt: The decision variable vector returned by MINCX/FEASP.
%  gain:

% XXX We set up and solve [G;I]'*PI*[G;I]<0 as one block LMI.  This
% means LMILab is doing the multiplications by I.  It would probably
% be faster to implement this as:
%    G'*PI11*G + G'*PI12 + PI12'*G + PI22<0
% In this form the mublock terms will look similar to mulmiub. Also,
% we won't have to do the re-shuffling.

% XXX This implementation assumes the IQC multiplier can be expressed
% as a linear combo of basis multipliers:  PI = sum_k x_k PI_k,  x_k>=0.
% This is sufficiently generic but it will lead to inefficient encoding
% in some cases. Specifically, some IQCs can be compactly parameterized
% affinely with matrices of free variables.  This implementation would
% require these types of IQC multipliers to be expanded in the linear
% combo form.

% AKP, Nov 2012.  Added last "perf" flag argument.  If this is TRUE, then the
% last UDYN is a performance IQC.  It is either 1-by-1 cell, or 1-by-2 cell.
% If 1-by-1, then that multiplier is NOT scaled, and the problem is only a
% FEASP.   If 1-by-2, then the first entry is not scaled, and the second entry is
% scaled by positive decision-variable, and that decision variable is also the
% objective to be minimized.  These are assumed to be (nd+ne)-by-(nd+ne), and
% are reshuffled as appropriate to fit correctly into the LMI solved here.
% When this is called, P has the performance IQC multipliers in it, and Pblk
% has the I/O dimension of the performance channels (these are both at the END
% of those arrays).

%
% MUSYNID

% Check for performance flag
nin = nargin;
if nin==6
    PerfFlag = false;
end

% Tolerances
tol = 0;            % IQC LMI tolerance
xtol = 1e-5;        % Tolerance for enforcing pos def of IQC multipliers

% Dimensions
[n,m,Nw] = size(G);
Nblk = size(blk,1);
cidx = unique(cidx);
Nc = length(cidx);
Nudyn = numel(IQCinfoS);
Pblk = reshape( [IQCinfoS.ExtBlkDim], [2,Nudyn])';

% Handle Performance IQCs
existobj = false;
if PerfFlag
    nd = Pblk(end,1);
    ne = Pblk(end,2);
    IQCPerf1 = IQCinfoS(end).PsiPi{1};
    
    % Rearrange the i/o channels and correct ths sign of performance IQC
    flipioch = [nd+1:nd+ne 1:nd];
    ConstantMultiplier = -IQCPerf1(flipioch,flipioch,:);  % 3-d double
    
    % Update Nudyn to remove the performance block
    Nudyn = Nudyn-1;
    
    if numel(IQCinfoS(end).PsiPi)==2
        IQCPerf2 = IQCinfoS(end).PsiPi{2};
        existobj = true;
        
        % Rearrange the i/o channels and correct ths sign of performance IQC
        ObjectiveMultiplier = -IQCPerf2(flipioch,flipioch,:);
        IDcenterPerf = IQCinfoS(end).IDcenter(1);
    end    
end

% mudim is the dimension of IQC multiplier for mu blocks
mudim = 0;
for k=1:Nblk
    if blk(k,2) == 0 % repeated real & repeated complex
        blksz = abs(blk(k,1))*[1 1];
    else
        blksz = blk(k,1:2);
    end
    mudim = mudim+blksz(1)+blksz(2);
end

% Shuffle [G;I] so that input/outputs for blocks are grouped together
% [In this ordering S(jw) is block diagonal conformable with Delta]
% Note: There may be extras signals associated with each IQC
GI = [G; repmat(eye(m),[1 1 Nw]) ];

Gptr = 0;
Iptr = 0;
Sptr = 0;
idx = [];
allblk = [blk;Pblk];
for k=1:size(allblk,1); % XXX Pblk includes perf block %(Nblk+Nudyn)
    if allblk(k,2) == 0 % repeated real & repeated complex
        blksz = abs(allblk(k,1))*[1 1];
    else
        blksz = allblk(k,1:2);
    end
    idx = [idx (Gptr+1):(Gptr+blksz(2))];
    Gptr = Gptr+blksz(2);
    idx = [idx n+( (Iptr+1):(Iptr+blksz(1)) ) ];
    Iptr = Iptr+blksz(1);
end
H = GI(idx,:,:);
muidx = idx(1:mudim);
udynidx = idx(mudim+1:end);

% LMI vars for mu blocks (rep. real, rep. complex, and full complex)
setlmis(lmisys);
vars = decnbr(lmisys);
if Nblk>0
    Sall = cell(Nc,1);
    Xall = cell(Nc,1);
    for i1=1:Nc
        R = [];
        I = [];
        X = [];
        for i2=1:Nblk
            if blk(i2,2) == 0 % repeated real & repeated complex
                blksz = abs(blk(i2,1));
                Xr  = symdec(blksz,vars);
                vars = vars+blksz*(blksz+1)/2;
                Xi = skewdec(blksz,vars);
                vars = vars+blksz*(blksz-1)/2;
                
                if blk(i2,1)<0  % repeated real
                    Yr = skewdec(blksz,vars);
                    vars = vars+blksz*(blksz-1)/2;
                    Yi  = symdec(blksz,vars);
                    vars = vars+blksz*(blksz+1)/2;
                else % repeated complex
                    Yr = zeros(blksz);
                    Yi = zeros(blksz);
                end
                R = blkdiag(R,[Xr Yr; Yr' -Xr]);
                I = blkdiag(I,[Xi Yi; -Yi' -Xi]);
                X = blkdiag(X,[Xr Xi; -Xi Xr]);
            else   % full block
                blksz = blk(i2,:);
                Xr = (vars+1)*eye(blksz(1));
                Xc = (vars+1)*eye(blksz(2));
                
                tmp = blkdiag(Xc,-Xr);
                R = blkdiag(R,tmp);
                I = blkdiag(I,zeros(size(tmp)));
                X = blkdiag(X,vars+1);
                
                vars = vars+1;
            end
        end
        Sall{i1} = lmivar(3,[R I; -I R]);
        Xall{i1} = lmivar(3,X);
    end
end

% Create LMI
% XXX LMI is linear in decision vars. Normalize one variable to be = 1?
% XXX Have to be a bit careful with the scaling because some vars
%  could be = 0 and trying to normalize them to be = 1 could create
%  some odd cases.
cnt = lminbr(lmisys)+1;
for i1=1:Nc
    H1 = H(1:mudim,:,cidx(i1));
    H2 = H(mudim+1:end,:,cidx(i1));
    H1E = [real(H1) imag(H1); -imag(H1) real(H1)];
    
    % H'*S*H <-tol*I (for mu blocks)
    if Nblk>0
        lmiterm([cnt,1,1, Sall{i1}],H1E',H1E);
    end
    lmiterm([-cnt,1,1, 0],-tol);
    
    % LMI terms for udyns
    ptr = 0;
    for i2=1:Nudyn
        % Grab IQCs associated with block i2
        IDcenter = IQCinfoS(i2).IDcenter;
        PsiPi  = IQCinfoS(i2).PsiPi;
        PsiFlag  = IQCinfoS(i2).PsiFlag;
        ni = numel(IDcenter);
        
        % Pull rows of [G;I] associated with this block.
        rdim = Pblk(i2,1)+Pblk(i2,2);
        H2i = H2(ptr+(1:rdim),:);
        ptr = ptr+rdim;
        
        for i3=1:ni
            if PsiFlag(i3)
                %fac = freqresp( PsiPi{i3}, omega(cidx(i1)) )*H2i;
                fac = PsiPi{i3}(:,:,cidx(i1))*H2i;
                R = real(fac);
                I = imag(fac);
                F1 = [R I];
                F2 = [-I R];
                
                lmiterm([cnt 1 1 IDcenter(i3)],F1',F1);
                lmiterm([cnt 1 1 IDcenter(i3)],F2',F2);
            else
                tmp = H2i'*PsiPi{i3}(:,:,cidx(i1))*H2i;
                tmpE = [real(tmp) imag(tmp); -imag(tmp) real(tmp)];
                lmiterm([cnt,1,1, IDcenter(i3)],tmpE/2,1,'s');
            end
        end
    end
    
    % LMI Terms for performance
    if PerfFlag  
        % Pull rows of [G;I] associated with perf channels
        rdim = nd+ne;
        H2i = H2(ptr+(1:rdim),:);

        % IO channels and Sign of ConstantMultiplier have already been adjusted
        tmp = H2i'*ConstantMultiplier(:,:,cidx(i1))*H2i;
        tmpE = [real(tmp) imag(tmp); -imag(tmp) real(tmp)];
        lmiterm([cnt,1,1, 0],(tmpE+tmpE')/2);        
        if existobj
            % Objective Term
            tmp = H2i'*ObjectiveMultiplier(:,:,cidx(i1))*H2i;
            tmpE = [real(tmp) imag(tmp); -imag(tmp) real(tmp)];
            lmiterm([cnt,1,1, IDcenterPerf],tmpE/2,1,'s');
        end
    end    
    cnt = cnt+1;
    
    % xtol*I < X
    if Nblk>0
        lmiterm([-cnt,1,1,Xall{i1}],1,1);
        lmiterm([cnt,1,1,0],xtol*eye(size(X)));
        cnt = cnt+1;
    end
end

% Solve LMI
lmisys = getlmis;
opts = [0 0 0 0 1];
%R = 1e4; opts = [0 0 R 0 1];
if existobj
    obj = zeros(vars,1);
    % XXXX the performance gain-variable, which multiplies the perfoamnce
    % IQC, is not in the last place, we think.  we think this does the
    % right indexing, but you need to check it.  This is the correction
    % that got the GAIN with UREALs working again.
    zzz = decinfo(lmisys,IQCinfoS(end).IDcenter(1));
    obj(zzz) = 1;
    opts(2) =300;
    [copt,xopt] = mincx(lmisys,obj,opts);
    if ~isempty(xopt)
        feas=1;
    else
        feas=0;
        S=[];
        xopt=[];
        perfparm = [];
        return;
    end
    perfparm = copt;
    
    %     % XXX PJS Play around with slightly increasing the gain used in the
    %     % IQC multiplier.  This seems to speed up convergence of the iter-
    %     % ation by making it more likely for the dense grid check to pass
    %     xopt(end) = 1.02*xopt(end);
else
    [tmin,xopt] = feasp(lmisys,opts);
    if tmin<=0
        feas = 1;
    else
        feas = 0;
    end
    perfparm = [];
end

% Form mu blocks of S
S = zeros([(n+m) (n+m) Nw]);
if Nblk>0
    % Use xopt at frequencies specified by cidx
    for i1=1:Nc
        Smu = dec2mat(lmisys,xopt,Sall{i1});
        S(muidx,muidx,cidx(i1)) = Smu(1:mudim,1:mudim)+...
            1j*Smu(1:mudim,mudim+1:end);
    end
    
    % Use linear interpolation at other frequencies
    % S is held constant for indices outside of [cidx(1),cidx(end)]
    %
    % XXX Also try to compute S at other frequencies using optimal
    % mussv scalings.  This should be easy to set up and run.
    didx = setdiff(1:Nw,cidx);
    for i1=1:length(didx)
        ptr = didx(i1);
        prev = max( cidx(cidx<ptr) );
        if isempty(prev)
            prev = ptr;
        end
        next = min( cidx(cidx>ptr) );
        if isempty(next)
            next = ptr;
        end
        
        % XXX linear interpolation in freq may not be the best choice
        % Possibly scl should be computed based on the log of freq?
        scl = ( omega(ptr)-omega(prev) ) / ( omega(next)-omega(prev) );
        S(muidx,muidx,ptr) = (1-scl)*S(muidx,muidx,prev) + ...
            scl*S(muidx,muidx,next);
    end
end

% Form udyn and perf blocks of S
for i1=1:Nw
    ptr = 0;
    for i2=1:Nudyn        
        IDcenter = IQCinfoS(i2).IDcenter;
        PsiFlag = IQCinfoS(i2).PsiFlag;
        PsiPi = IQCinfoS(i2).PsiPi;
        ni = numel(IDcenter);
        blksz = sum(Pblk(i2,:));
        Stmp = zeros(blksz);
        
        for i3=1:ni
            if PsiFlag(i3)
                %Psitmp = freqresp( PsiPi{i3}, omega(i1) );
                Psitmp = PsiPi{i3}(:,:,i1);
                Mtmp = dec2mat(lmisys,xopt,IDcenter(i3));
                Stmp = Stmp+Psitmp'*Mtmp*Psitmp;
            else
                xtmp = dec2mat(lmisys,xopt,IDcenter(i3));
                Stmp = Stmp + xtmp*PsiPi{i3}(:,:,i1);
            end
        end
        
        blkidx = udynidx( ptr+(1:sum(blksz)) );
        S(blkidx,blkidx,i1)=Stmp;
        ptr = ptr+sum(blksz);
    end
    if PerfFlag
        blksz = ne + nd;
        blkidx = udynidx( ptr+(1:sum(blksz)) );
        S(blkidx,blkidx,i1)=S(blkidx,blkidx,i1)+ConstantMultiplier(:,:,i1);
        if existobj
            xtmp = dec2mat(lmisys,xopt,IDcenterPerf);
            S(blkidx,blkidx,i1)=S(blkidx,blkidx,i1)+xtmp*ObjectiveMultiplier(:,:,i1);
        end
    end    
end
