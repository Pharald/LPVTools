function [Wn,Z,P] = damp(varargin)
% DAMP  Pointwise compute natural frequency and damping of PLFTSS objects.
% 
% [Wn,Z] = DAMP(SYS,DOMAIN) returns vector PMATs Wn and Z containing the 
% natural  frequencies and damping factors of SYS at each point in the 
% domain specified by the RGRID object DOMAIN. For discrete-time models, 
% the equivalent s-plane natural frequency and damping ratio of an 
% eigenvalue lambda are:
%     Wn = abs(log(lambda))/Ts ,   Z = -cos(angle(log(lambda))) 
% Wn and Z are empty vectors if the sample time Ts is undefined.
%
% [Wn,Z,P] = DAMP(SYS,DOMAIN) also returns the poles P of SYS evaluated
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
    if isa(varargin{i},'plftss')
        varargin{i} = lft2grid(varargin{i},DOM);
    end
end

[Wn,Z,P]=damp(varargin{:});