function [Wn,Z,P] = damp(varargin)
% DAMP  Pointwise compute natural frequency and damping of PLFTMAT objects.
% 
% [Wn,Z] = DAMP(M,DOMAIN) returns vector PMATs Wn and Z containing the natural 
% frequencies and damping factors of the PLFTMAT M at each point in the 
% domain specified by the RGRID object DOMAIN.
%
% [Wn,Z,P] = DAMP(M,DOMAIN) also returns the poles P of M evaluated
% at each point in the domain specified by the RGRID object DOMAIN.
%
% See also: damp, pole, zero.


DOM = varargin{end};
varargin(end) = [];

if ~isa(DOM,'rgrid')
    error(['Last input argument must be an RGRID object that specifies'...
           ' the parameter domain grid.'])
end

for i = 1:numel(varargin)
    if isa(varargin{i},'plftmat')
        varargin{i} = lft2grid(varargin{i},DOM);
    end
end

[Wn,Z,P]=damp(varargin{:});