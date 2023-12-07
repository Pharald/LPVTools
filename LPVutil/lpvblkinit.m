function [nU,ABMat,CDMat,IVData,dIVData,niv,LIVData] = lpvblkinit(sys,plist)
% LPVBLKINIT   Initialize the Simulink LPVTools blocks.

if isa(sys,'pss')
    % Mask Initialization for PSS
    szS = size(sys);
    nU = szS(2);
    nad = numel(szS)-2;
    if nad>0
        error('LPVBLK does not permit systems with array dimensions.');
    end
    IVData = sys.Domain.IVData;
    niv = length(IVData);
    dIVData = sys.Domain.DIVData;
    LIVData = sys.Domain.LIVData;
    IVName = sys.Domain.IVName;
    [a,b,c,d] = ssdata(sys);
    ABpmat = [a b];
    CDpmat = [c d];
    [tfidx,loc] = ismember(plist,IVName);
    if all(tfidx) && numel(plist)==numel(IVName)
        IVData = IVData(loc);
        dIVData = dIVData(loc);
        LIVData = LIVData(loc);
        ABMat = permute(ABpmat.Data,[1 2 2+loc']);
        CDMat = permute(CDpmat.Data,[1 2 2+loc']);
    else
        error('ParameterVectorInputList must match parameters in system')
    end
elseif isa(sys,'plftss')    
    % Mask Initialization for PLFTSS
    szS = size(sys);
    nU = szS(2);
    niv = numel(plist);
    
    % dIVData is unused. Set as NaN to flag the PLFTSS case.
    dIVData = nan; 

    % Get parameter ranges
    IVName = fieldnames(sys.Parameter);
    [tfidx,loc] = ismember(plist,IVName);
    if all(tfidx) && numel(plist)==numel(IVName)
        IVData = cell(numel(loc),1);
        for i=1:numel(plist)
           IVData{i} = sys.Parameter.(IVName{loc(i)}).Range;
        end
    else
        error('Input parameter list must match parameters in system')
    end
        
     % Get SS Data
    [a,b,c,d]=ssdata(sys);
    ABpmat = [a b];
    CDpmat = [c d];
    
    % Pull apart LFT
    [Mab,Delab,Blkab]=lftdata(ABpmat,[],'Parameters');
    [Mcd,Delcd,Blkcd]=lftdata(CDpmat,[],'Parameters');
    if isuncertain(Mab) || isuncertain(Mcd)
        error('Input system may not include uncertainties');
    end    
    
    % Process blocks for Mab to easily compute 
    %   Fu(M,Del(rho)) = M22 + M21*Del*inv(I-M11*Del)*M12
    % Note that Delta(rho) will be square and block diag. 
    % Store partitioned 2-by-2 matrices of M in ABmat.
    % Store indices to construct Del from rho in LIVData
    idxab = [];
    for i=1:numel(Blkab)
        [~,idx] = ismember(Blkab(i).Name,plist);
        if idx==0
            % XXX Can this happen? There is already a check above that 
            % PLIST matches the list of TVREALs
            error(['Variable ' Blkab(i).Name ' does not appear in input parameter list.']);
        end
        idxab = [idxab, repmat(idx,[1 Blkab(i).Occurrences])];
    end    
    r = size(Delab,1);
    Mab11 = Mab(1:r,1:r);
    Mab12 = Mab(1:r,r+1:end);
    Mab21 = Mab(r+1:end,1:r);
    Mab22 = Mab(r+1:end,r+1:end);    
    ABMat = {Mab11, Mab12, Mab21, Mab22};
    LIVData{1} = idxab;
    
    % Process blocks for Mcd to easily compute Fu(Mcd,Delcd)
    idxcd = [];
    for i=1:numel(Blkcd)
        [~,idx] = ismember(Blkcd(i).Name,plist);
        if idx==0
            error(['Variable ' Blkcd(i).Name ' does not appear in input parameter list.']);
        end
        idxcd = [idxcd, repmat(idx,[1 Blkcd(i).Occurrences])];
    end    
    r = size(Delcd,1);
    Mcd11 = Mcd(1:r,1:r);
    Mcd12 = Mcd(1:r,r+1:end);
    Mcd21 = Mcd(r+1:end,1:r);
    Mcd22 = Mcd(r+1:end,r+1:end);    
    CDMat = {Mcd11, Mcd12, Mcd21, Mcd22};
    LIVData{2} = idxcd;
else
    error('Input system must be a PSS or PLFTSS object.');
end