function K = dcgain(varargin)
% DCGAIN  Pointwise DC gain of a PLFTSS object
%
% K = DCGAIN(SYS,DOMAIN) computes the steady-state gain of SYS at each 
% point in the domain specified by the RGRID object DOMAIN. If SYS has NY 
% outputs and NU inputs then K is returned as an NY-by-NU PMAT containing 
% the steady-state gains of SYS at each point in the domain.
%
% See also: dcgain, norm, evalfr, freqresp.

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

K = dcgain(varargin{:});