function [GAP,NUGAP] = gapmetric(varargin)
% gapmetric Pointwise gap gap and the Vinnicombe gap metric for PLFTSS objects
%
%[GAP,NUGAP] = gapmetric(P1,P2,DOMAIN) calculates the gap metric between 
% the systems P1 and P2 at each point in the domain specified by the rgrid 
% object DOMAIN. See lti/gapmetric for details.
% 
% See also: norm, loopmargin, wcmargin.

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

[GAP,NUGAP] = gapmetric(varargin{:});