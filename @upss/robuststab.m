function [STABMARG,DESTABUNC,REPORT,INFO] = robuststab(sys,opt)
%ROBUSTSTAB   Pointwise robust performance margins of a UPSS.
%
% [STABMARG,DESTABUNC,REPORT,INFO] = ROBUSTSTAB(SYS) computes the robust 
% stability margins for the uncertain parameter-varying system SYS at each
% point in the domain of SYS. The outputs STABMARG, DESTABUNC, REPORT, 
% and INFO are PSTRUCTs describing the results at each point on the domain 
% of SYS. STABMARG contains the bounds on the stability margin. 
% DESTABUNC contains the value of uncertainty which leads to instability. 
% REPORT is a string of text describing the robust stability analysis. 
% INFO contains additional results from the analysis. 
% See DynamicSystem/robustperf for details.
%
% [STABMARG,DESTABUNC,REPORT,INFO] = ROBUSTSTAB(SYS,OPTS) allows the user 
% to pass in a ROBUSTSTABOPTIONS object OPTS.
%
% See also: robuststab, robuststabOptions, robustperf, wcgain, loopsens.

% Define common variables
Dom = sys.Domain;

% Define default options object.
if nargin == 1
    opt = robustperfOptions();
end


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

for i = 1:Npts
    if i==1        
        [STABMARG,DESTABUNC,rpt,INFO] = robuststab(Data(idx{:},i),opt);
        % rpt is a char not a struct. We want to organize the reports in 
        % rpt into a structure that can be formed into a PSTRUCT. 
        % Make a struct called REPORT to place the reports in rpt into, by 
        % creating an empty struct called REPORT and placing rpt into a 
        % field called "REPORT.Report".        
        RPTtemplate = rmfield(STABMARG,fieldnames(STABMARG));
        rptad = prod(size(RPTtemplate));
        REPORT = RPTtemplate;
        for k = 1:rptad
            REPORT(k).Report = rpt(:,:,k);
        end
    else
        [STABMARG(idxad{:},i) DESTABUNC(idxad{:},i) rpt INFO(idxad{:},i)] = robuststab(Data(idx{:},i),opt);        
        Rtemp = RPTtemplate;
        for k = 1:rptad
            Rtemp(k).Report = rpt(:,:,k);
        end
        REPORT(idxad{:},i) = Rtemp;
    end
end

szsys = size(sys);
szad = szsys(3:end);
if nad==0
    szad = [1 1];
end


% Reshape to [row,col,AD,IV] from [row,col,AD,prod(IV)] 
STABMARG = reshape(STABMARG, [szad LIVData]);
DESTABUNC = reshape(DESTABUNC, [szad LIVData]);
REPORT = reshape(REPORT,[szad LIVData]);
INFO = reshape(INFO,[szad LIVData]);

% Reorder from [row, col, AD, IV] to [row,col,IV,AD]
numAD = numel(szad(3:end));
STABMARG  = permute(STABMARG,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);
DESTABUNC = permute(DESTABUNC,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);
REPORT    = permute(REPORT,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);
INFO      = permute(INFO,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);

% Form PSTRUCT to hold results
STABMARG     = pstruct(STABMARG,Dom);
DESTABUNC    = pstruct(DESTABUNC,Dom);
REPORT       = pstruct(REPORT,Dom);
INFO         = pstruct(INFO,Dom);
    
end
