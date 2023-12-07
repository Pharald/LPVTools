function [WCG,WCU,Info] = wcgain(varargin)
% WCGAIN   Pointwise calculates worst-case gain of a PLFTSS.
% 
% [WCG,WCU,INFO] = wcgain(SYS,...,DOMAIN) computes the worst-case gain of 
% SYS at each point in the domain specified by the RGRID object DOMAIN. 
% See DynamicSystem/wcgain for details.
%
% see also: wcgain, robuststab, robustperf, loopsens.


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

[WCG,WCU,Info] = wcgain(varargin{:});