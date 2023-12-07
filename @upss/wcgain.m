function [WCG,WCU,Info] = wcgain(sys,varargin)
% WCGAIN   Pointwise calculates worst-case gain of a UPSS.
% 
% [WCG,WCU,INFO] = wcgain(SYS,...) computes the worst-case gain of SYS
% at each point in the domain of SYS. See DynamicSystem/wcgain for details.
%
% see also: wcgain, robuststab, robustperf, loopsens.

% Define common variables
Dom = sys.Domain;

% Look for options object
idxOpt = find(cellfun(@(x) isa(x,'rctoptions.wcgain'),varargin));
if numel(idxOpt)==1
    opt = varargin{idxOpt};
else
    opt = wcgainOptions();
end
MaxOverFrequency = strcmp(opt.MaxOverFrequency,'on');
MaxOverArray = strcmp(opt.MaxOverArray,'on');


% Reorder as [AD IV]
niv  = Dom.NumIV;
nad = numel(size(sys))-2;
if nad>0
    Data = permute(sys.Data,[(niv+1:niv+nad) (1:niv)]);
else
    Data = sys.Data;
end
LIVData = [Dom.LIVData]';
Npts = prod(LIVData);
idx = cellstr(repmat(':',nad+2,1))';
if nad==0
    idxad  = {':'};
else
    idxad = cellstr(repmat(':',nad,1))';
end

% Perform wcgain analysis at each point in the domain
if MaxOverArray
    for i = 1:Npts
        if i==1
            [WCG,WCU,Info] = wcgain(Data(idx{:},i),varargin{:});
        else
            [WCG(i),WCU(i),Info(i)] = wcgain(Data(idx{:},i),varargin{:});
        end
    end
    WCG  = pstruct(reshape(WCG,[1,1,LIVData]),Dom);
    WCU  = pstruct(reshape(WCU,[1,1,LIVData]),Dom);
    Info = pstruct(reshape(Info,[1,1,LIVData]),Dom);
    
else
    for i = 1:Npts
        if i==1
            [WCG,WCU,Info] = wcgain(Data(idx{:},i),varargin{:}); % 
        else
            [WCG(idxad{:},i),WCU(idxad{:},i),Info(idxad{:},i)] = wcgain(Data(idx{:},i),varargin{:});
        end
    end
    szsys = size(sys);
    szad = szsys(3:end);
    if nad==0
        szad = [1 1];
    end
    
    % Reshape from [row,col,AD,prod(IV)] to [row,col,AD,IV] 
    WCG = reshape(WCG, [szad LIVData]);
    WCU = reshape(WCU, [szad LIVData]);
    Info = reshape(Info,[szad LIVData]);
    
    % Reorder from [AD, IV] to [AD1,AD2,IV,AD3toN]
    numAD = numel(szad(3:end));
    WCG    = permute(WCG,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);
    WCU = permute(WCU,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);
    Info      = permute(Info,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);
    
    % Form PSTRUCT to hold results
    WCG     = pstruct(WCG,Dom);
    WCU     = pstruct(WCU,Dom);
    Info    = pstruct(Info,Dom);
end

