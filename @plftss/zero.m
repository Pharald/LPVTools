function [Z,GAIN] = zero(varargin)
% ZERO  Pointwise computes the transmission zeros of a PLFTSS object
%  
% Z = ZERO(SYS,DOMAIN) returns the transmission zeros of SYS at each point in
% the domain specified by the RGRID object DOMAIN. Z is returned as a 
% column vector PMAT containing the transmission zeros at each point in 
% the domain of SYS.
%
% [Z,GAIN] = ZERO(SYS,DOMAIN) also returns the transfer function gain for
% SISO models SYS.
%
% See also: zero, damp, pole.

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
if nargout>1
    [Z,GAIN]=zero(varargin{:});
else
    [Z]=zero(varargin{:});
end