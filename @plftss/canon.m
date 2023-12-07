function [csys,T] = canon(varargin)
% CANON  Pointwise canonical state-space realizations for PLFTSS objects 
%
% CSYS = CANON(SYS,TYPE,DOMAIN) computes a canonical state-space 
% realization CSYS at each point in the domain specified by the RGRID
% object DOMAIN. The string TYPE selects the type of realization and can 
% be 'modal' (default) or 'companion'.
%
% [CSYS,T] = CANON(SYS,TYPE,DOMAIN) also returns the transformation 
% PMAT T that relates the canonical state vector z to the original 
% state vector x, by z = Tx at each point in the domain of SYS.
%
% CSYS = CANON(SYS,'modal',CONDT,DOMAIN) specifies an upper bound CONDT on 
% the condition number of the block-diagonalizing transformation T. 
%
% See also: canon, ss2ss.

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

[csys,T] = canon(varargin{:});