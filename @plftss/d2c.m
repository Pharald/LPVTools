function sysc = d2c(sysd,varargin)
% D2C   Point-wise discrete-time to continuous-time conversion for PLFTSS
%
% SYSC = d2c(SYSD,METHOD,DOMAIN) performs the discrete- to continuous-time
% convertion at each point in the domain specified by the RGRID object 
% DOMAIN. The string METHOD selects the discretization method among 
% 'zoh' (default), 'tustin', and 'matched'.
%
% D2C(SYSC,OPTIONS,DOMAIN) gives access to additional discretization 
% options. Use D2COPTIONS to create and configure the option set OPTIONS.
%
% See also: d2c, d2cOptions, c2d.

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

[sysd,G] = d2c(varargin{:})
