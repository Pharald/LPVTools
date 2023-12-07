function [perfparm,Sopt,xopt] = iqcperf4(A,omega,perfInfo)

% Check for gain flag
nin = nargin;
if nin==2
    perfInfo = 'gain';
end
Nomega = numel(omega);

% Re-order uncertainties so that udyns appear last
[G, allPiC, mublk, Nudyn, udynsize, udyncop] = udyn2end(A);
[ne,nd] = size(A);

% Augment allPiC with the correct performance IQC.  Here, the performance
% IQC is as the user would enter it, i.e. to solve
%   min alpha
%   s.t. [d;e]*(PIperf1+ alpha*Piperf2)*[d;e] >=0
% Later (iqcsolve) it gets rearranged with the correct sign for the LMI
gainCalc = false;
if isequal(perfInfo,'gain') || isequal(perfInfo,{'gain'})
    perfInfo = {blkdiag(zeros(nd),-eye(ne)) blkdiag(eye(nd),zeros(ne))};
    gainCalc = true;
elseif numel(perfInfo)>2 || ~isa(perfInfo,'cell')
    % XXX Also allow struct array (function handle, etc) form?
    error('Invalid PerfInfo');
end
allPiC{Nudyn+1} = IQCcell(perfInfo,true);
udynsize = [udynsize;nd ne];
udyncop = [udyncop;1];
Nudyn = Nudyn + 1;

% Initialize IQCInfoS structured array
IQCinfoS = struct('IDcenter',[],'PsiPi',{cell(0,1)},...
        'PsiFlag',{true(0,1)},'OrigBlkDim',[],'ExtBlkDim',[]);
IQCinfoS = repmat(IQCinfoS,[Nudyn 1]);

% Evaluate IQCs for udyns
setlmis([]);
lmisys = getlmis;
% PsiPi= cell(Nudyn,1);
% IDcenter= cell(Nudyn,1);
% PsiFlag = false(Nudyn,1);
% Pblk = zeros(Nudyn,2);
allsigdim = zeros(Nudyn,1);
for k=1:Nudyn
    % Grab User Data associated with k^th block (struct array)
    UserDatakS = allPiC{k};
    
    nfhs = numel(UserDatakS);
    sigdim = zeros(nfhs,1);
    origPidim = udyncop(k)*sum(udynsize(k,:));
    nIQC = zeros(nfhs,1);
    for j=1:nfhs        
        % Each fh can return multiple IQCs
        fh = UserDatakS(j).IQCfunction;
        param = UserDatakS(j).IQCparams;
        PsiFlag = UserDatakS(j).PsiFlag;
        
        % Convert fh description to a list of IQCs
        [lmisys,IDcenterTMP,PsiPiTMP] = fh(lmisys,udyncop(k),omega,param);
        IQCinfoS(k).IDcenter = [IQCinfoS(k).IDcenter; IDcenterTMP(:)];
        IQCinfoS(k).PsiPi = [IQCinfoS(k).PsiPi; PsiPiTMP(:)];
        nIQC(j) = numel(PsiPiTMP);
        IQCinfoS(k).PsiFlag = [IQCinfoS(k).PsiFlag; repmat(PsiFlag,[nIQC(j) 1])];        
        
        % Handle IQCs with extra signal dimensions
        % XXX Assumes that all values of PsiTMP have the same col dim
        tmp = cellfun(@size,PsiPiTMP,'UniformOutput',false);
        tmp = cell2mat( tmp(:) );
        sigdim(j) = tmp(1,2) - origPidim;
    end
    
    % Reshape the IQCs to include any extra signals assuming
    % all extra signals are included at the end
    sumsigdim = cumsum([0; sigdim(:)]);
    allsigdim(k) = sumsigdim(end);
    finalPidim = origPidim + sumsigdim(end);
    Itmp = eye(finalPidim);
    ptr = 0;
    for j = 1:nfhs
        idx = origPidim + 1 + (sumsigdim(j):sumsigdim(j+1)-1);
        Mscl = Itmp([1:origPidim idx],:);
        PsiFlag = UserDatakS(j).PsiFlag;
        nMscl = size(Mscl,2);        
        for i=1:nIQC(j)
            if PsiFlag
                Psitmp = zeros( size(IQCinfoS(k).PsiPi{i+ptr},1), nMscl, Nomega);
                for n=1:Nomega
                    Psitmp(:,:,n) = IQCinfoS(k).PsiPi{i+ptr}(:,:,n)*Mscl;
                end
                IQCinfoS(k).PsiPi{i+ptr} = Psitmp;
            else
                Pitmp = zeros( nMscl , nMscl, Nomega);
                for n=1:Nomega
                    Pitmp(:,:,n) = Mscl'*IQCinfoS(k).PsiPi{i+ptr}(:,:,n)*Mscl;
                end
                IQCinfoS(k).PsiPi{i+ptr} = Pitmp;
            end
        end
        
%         for i=1:nIQC(j)
%             for n=1:size(IQCinfoS(k).PsiPi{i+ptr},3)
%                 if PsiFlag
%                     IQCinfoS(k).PsiPi{i+ptr}(:,:,n) = IQCinfoS(k).PsiPi{i+ptr}(:,:,n)*Mscl;
%                 else
%                     IQCinfoS(k).PsiPi{i+ptr}(:,:,n) = Mscl'*IQCinfoS(k).PsiPi{i+ptr}(:,:,n)*Mscl;
%                 end
%             end
%         end
        ptr = ptr+nIQC(j);
    end
    
    % Store updated block dimensions accounting for any extra signals
    c = udynsize(k,2)*udyncop(k);
    IQCinfoS(k).OrigBlkDim = udyncop(k)*udynsize(k,:);
    IQCinfoS(k).ExtBlkDim = [finalPidim-c c];
end

% Create matrix to expand G to account for extra IQC signals
tmpsz = sum(udyncop.*udynsize(:,1));
Gscl = zeros( tmpsz, tmpsz+sum(allsigdim) );
rptr = 0;
cptr = 0;
for k=1:Nudyn
    dim = udyncop(k)*udynsize(k,1);
    Id = eye(dim);
    ridx = rptr + (1:dim);
    cidx = cptr + (1:dim);
    Gscl(ridx,cidx) = Id;
    rptr = rptr+dim;
    cptr = cptr+dim+allsigdim(k);
end
Gscl = blkdiag( eye(sum(abs(mublk(:,1)))) , Gscl);
G = G*Gscl;

% Get frequency and response data of G
% XXX Need to be careful about frequency units here. Check with other code.
Gg = freqresp(G,omega);
szGg = size(Gg);
GI = [Gg; repmat(eye(szGg(2)),[1 1 Nomega]) ];

% Iteration: Solve IQC LMI Feasibility problem on coarse frequency grid,
% check results on dense grid, and add frequencies as needed.
if Nomega<=2
    cidx = 1:Nomega;
else
    % Endpoints are needed to ensure proper interpolation.
    cidx = unique([1 Nomega]);
end
%cidx = 1:Nomega;  % Uncomment to use all frequency points in iqcsolve

go = 1;
cnt = 0;
feas = 0;
Sopt = [];
xopt = [];
ev = zeros(Nomega,1);
PerfFlag = true;
while go==1
    % Solve feasibility problem on coarse grid
    cnt = cnt+1;
    [cfeas,S,xopt,perfparm] = iqcsolve4(Gg, mublk, IQCinfoS, lmisys, omega, cidx, PerfFlag);
    if cfeas==0
        % Infeasible on coarse grid: Stop
        go = 0;
        feas = 0;
        newidx = -1;
    else
        % Check IQC LMI on dense grid
        % XXX Perhaps we should be wrapping in the mutlitpliers for true
        % uncertainty blocks and then computing the gain at each freq?
        % (instead of checking the stab condition at each freq)
        for i1=1:Nomega
            IQClmi = GI(:,:,i1)'*S(:,:,i1)*GI(:,:,i1);
            ev(i1)=max(real(eig(IQClmi)));
        end
        fidx = find(ev>=0);  % XXX Add tolerance?
        
        if isempty(fidx)
            % S is feasible on dense grid: Stop
            feas = 1;
            Sopt = S;
            go = 0;
            newidx = -1;
        else
            % Add new frequencies
            % XXX For now, simply add freq corresponding to worst
            %    constraint violation.
            [maxev,newidx] = max(ev);
            if any(cidx==newidx)
                % XXX This would mean the LMI was feasible with constraint
                % newidx in the list cidx but the subsequent IQC LMI check
                % found that newidx is the constraint with the worst
                % violation.
                % This should never happend but if it does then some
                % numerical issue has occured.
                feas = 0;
                go = 0;
            else
                cidx = [cidx newidx];
            end
        end
    end
    
    % Diagnostic Information
    %     fprintf(' cnt = %d \t cfeas = %d Nomega = %d Ncidx=%d\t newidx=%d\n',...
    %         cnt,cfeas,Nomega,length(cidx),newidx);
    
    %figure(1); plot(1:Nomega,ev); title(num2str(cnt)); keyboard
end

if gainCalc
    perfparm = sqrt(perfparm);
end

if feas==0;
    perfparm = [];
    Sopt = [];
    xopt = [];
end