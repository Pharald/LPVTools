function [Gm,Pm,Wcg,Wcp] = margin(varargin)
% MARGIN  Pointwise gain and phase margins, and crossover frequencies for a PLFTSS
%  
% [Gm,Pm,Wcg,Wcp] = MARGIN(SYS,DOMAIN) computes the gain margin Gm, phase 
% margin Pm, and the associated frequencies Wcg and Wcp, for the SISO
% open-loop PLFTSS SYS (continuous or discrete).  The margins are computed
% at each point in the domain specified by the RGRID object DOMAIN.
% Gm, Pm, Wcg, Wcp are returned as PMATs that provide the margins and 
% crossover frequencies at each point in the domain specified in DOMAIN. 
% See DynamicSystem/MARGIN for the definitions of the gain and phase margin.
%
% See also: margin.

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

[Gm,Pm,Wcg,Wcp] = margin(varargin{:});