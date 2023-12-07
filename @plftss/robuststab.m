function [STABMARG,DESTABUNC,REPORT,INFO] = robuststab(varargin)
%ROBUSTSTAB   Pointwise robust performance margins of a PLFTSS.
%
% [STABMARG,DESTABUNC,REPORT,INFO] = ROBUSTSTAB(SYS,DOMAIN) computes the robust 
% stability margins for the uncertain parameter-varying system SYS at each
% point in the domain specified by the RGRID object DOMAIN. The outputs 
% STABMARG, DESTABUNC, REPORT, and INFO are PSTRUCTs describing the results 
% at each point on the domain of SYS. STABMARG contains the bounds on the 
% stability margin. DESTABUNC contains the value of uncertainty which leads 
% to instability. REPORT is a string of text describing the robust stability 
% analysis. INFO contains additional results from the analysis. 
% See DynamicSystem/robustperf for details.
%
% [STABMARG,DESTABUNC,REPORT,INFO] = ROBUSTSTAB(SYS,OPTS,DOMAIN) allows 
% the user to pass in a ROBUSTSTABOPTIONS object OPTS.
%
% See also: robuststab, robuststabOptions, robustperf, wcgain, loopsens.


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

[STABMARG,DESTABUNC,REPORT,INFO] = robuststab(varargin{:});