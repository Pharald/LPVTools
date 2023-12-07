function [PERFMARG,PERFMARGUNC,REPORT,INFO] = robustperf(sys,opt)
%ROBUSTPERF   Pointwise robust performance margins of a UPSS.
%
% [PERFMARG,PERFMARGUNC,REPORT,INFO] = ROBUSTPERF(SYS) computes the robust 
% performance margin of the uncertain parameter-varying system SYS at each
% point in the domain of SYS. The outputs PERFMARG, PERFMARGUNC, REPORT, 
% and INFO are PSTRUCTs describing the results at each point on the domain 
% of SYS. PERFMARG contains the bounds on the performance margin. 
% PERFMARGUNC contains the value of the worst-case uncertainty. REPORT is a
% string of text describing the robust performance analysis. INFO contains 
% additional results from the analysis. See DynamicSystem/robustperf for 
% details.
%
% [PERFMARG,PERFMARGUNC,REPORT,INFO] = ROBUSTPERF(SYS,OPTS) allows the user 
% to pass in a ROBUSTPERFOPTIONS object OPTS.
%
% See also: robustperf, robustperfOptions, robuststab, wcgain, loopsens.


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

% Version 1 - Forcing REPORT to be a struct by creating a dummy field:
for i = 1:Npts
    if i==1        
        [PERFMARG,PERFMARGUNC,rpt,INFO] = robustperf(Data(idx{:},i),opt);
        % rpt is a char not a struct. We want to organize the reports in 
        % rpt into a structure that can be formed into a PSTRUCT. 
        % Make a struct called REPORT to place the reports in rpt into, by 
        % creating an empty struct called REPORT and placing rpt into a 
        % field called "REPORT.Report".        
        RPTtemplate = rmfield(PERFMARG,fieldnames(PERFMARG));
        rptad = prod(size(RPTtemplate));
        REPORT = RPTtemplate;
        for k = 1:rptad
            REPORT(k).Report = rpt(:,:,k);
        end
    else
        [PERFMARG(idxad{:},i) PERFMARGUNC(idxad{:},i) rpt INFO(idxad{:},i)] = robustperf(Data(idx{:},i),opt);        
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
PERFMARG = reshape(PERFMARG, [szad LIVData]);
PERFMARGUNC = reshape(PERFMARGUNC, [szad LIVData]);
REPORT = reshape(REPORT,[szad LIVData]);
INFO = reshape(INFO,[szad LIVData]);

% Reorder from [row, col, AD, IV] to [row,col,IV,AD]
numAD = numel(szad(3:end));
PERFMARG    = permute(PERFMARG,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);
PERFMARGUNC = permute(PERFMARGUNC,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);
REPORT      = permute(REPORT,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);
INFO        = permute(INFO,[1 2 (3+numAD):(2+niv+numAD) (3:2+numAD)]);

% Form PSTRUCT to hold results
PERFMARG     = pstruct(PERFMARG,Dom);
PERFMARGUNC  = pstruct(PERFMARGUNC,Dom);
REPORT       = pstruct(REPORT,Dom);
INFO         = pstruct(INFO,Dom);
end


% % Version 2, returning a REPORT that is a char array:
%
% for i = 1:Npts
%     if i==1
%         % TODO - Original from visit to Minneapolis
%         %             [WCG(idxad{:}),WCU(idxad{:}),Info(idxad{:})] = wcgain(Data(idx{:},i),varargin{:});
%         % TODO - This works for simple example in Q2 report.
%         % [PERFMARG(idxad{:}),PERFMARGUNC(idxad{:}),REPORT(idxad{:}),INFO(idxad{:})] = robustperf(Data(idx{:},i),opt)
%         [PERFMARG,PERFMARGUNC,REPORT,INFO] = robustperf(Data(idx{:},i),opt);
%         % REPORT is a char not a struct. Hence it needs to be the same
%         % size for each iteration. Contents of the report will change
%         % for each point that is analyzed. Solution is to make pad the
%         % report to be larger than is needed and then enforce all the
%         % reports to be of that overlarge size.
%         spz = '     ';
%         REPORT = [REPORT [spz;spz;spz;spz;spz]];
%         szRT = size(REPORT);
%     else
%         %             [PERFMARG(idxad{:},i),PERFMARGUNC(idxad{:},i),REPORT(idxad{:},i),INFO(idxad{:},i)] = robustperf(Data(idx{:},i),opt)
%         [PERFMARG(idxad{:},i) PERFMARGUNC(idxad{:},i) rpt INFO(idxad{:},i)] = robustperf(Data(idx{:},i),opt);
%         % resize rpt to be same size as REPORT
%         szrpt = size(rpt);
%         diffrpt = szRT(2)-szrpt(2);
%         buff = repmat(' ',5,diffrpt);
%         augrpt = [rpt buff];
%         REPORT(idx{:},i) = augrpt;
%         
%     end
% end
% szsys = size(sys);
% szad = szsys(3:end);
% if nad==0
%     szad = [1 1];
% end
% 
% PERFMARG  = pstruct(reshape(PERFMARG, [szad,LIVData]),Dom);
% PERFMARGUNC  = pstruct(reshape(PERFMARGUNC, [szad,LIVData]),Dom);
% REPORT = reshape(REPORT,[szRT,LIVData]);
% INFO = pstruct(reshape(INFO,[szad,LIVData]),Dom);
% end
