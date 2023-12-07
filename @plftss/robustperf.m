function [PERFMARG,PERFMARGUNC,REPORT,INFO] = robustperf(varargin)
%ROBUSTPERF   Pointwise robust performance margins of a PLFTSS.
%
% [PERFMARG,PERFMARGUNC,REPORT,INFO] = ROBUSTPERF(SYS,DOMAIN) computes the  
% robust performance margin of the uncertain parameter-varying system SYS 
% at each point in the domain specified by the RGRID object DOMAIN. The 
% outputs PERFMARG, PERFMARGUNC, REPORT, and INFO are PSTRUCTs describing 
% the results at each point on the domain of SYS. PERFMARG contains the 
% bounds on the performance margin. PERFMARGUNC contains the value of the 
% worst-case uncertainty. REPORT is a string of text describing the robust 
% performance analysis. INFO contains additional results from the analysis. 
% See DynamicSystem/robustperf for details.
%
% [PERFMARG,PERFMARGUNC,REPORT,INFO] = ROBUSTPERF(SYS,OPTS,DOMAIN) allows 
% the user to pass in a ROBUSTPERFOPTIONS object OPTS.
%
% See also: robustperf, robustperfOptions, robuststab, wcgain, loopsens.


DOM = varargin{end};
varargin(end) = [];

if ~isa(DOM,'rgrid')
    error(['Last input argument must be an RGRID object that specifies'...
           ' the parameter domain grid.'])
end

for i = 1:numel(varargin)
    if isa(varargin{i},'plftss')
        varargin{i} = lft2grid(varargin{i},DOM);
    end
end

[PERFMARG,PERFMARGUNC,REPORT,INFO] = robustperf(varargin{:});