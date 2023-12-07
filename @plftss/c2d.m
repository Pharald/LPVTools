function [sysd,G] = c2d(varargin)
% C2D   Point-wise continuous-time to discrete-time conversion for PLFTSS
%
% SYSD = c2d(SYSC,TS,METHOD,DOMAIN) performs the continuous to discrete
% convertion at each point in the domain specified by the RGRID object 
% DOMAIN. SYSD is a discrete-time PSS with sampling time TS that 
% approximates the continuous-time PSS SYSC. The string METHOD selects the 
% discretization method among 'zoh' (default), 'foh', 'impulse', 'tustin', 
% and 'matched'.
%
% C2D(SYSC,TS,OPTIONS,DOMAIN) gives access to additional discretization 
% options. Use C2DOPTIONS to create and configure the option set OPTIONS.
%
% For state-space models without delays, [SYSD,G] = C2D(SYSC,Ts,METHOD,DOMAIN), 
% also returns the matrix G mapping the states xc(t) of SYSC to the states 
% xd[k] of SYSD: xd[k] = G * [xc(k*Ts) ; u[k]]  
%
% See also: c2d, c2dOptions, d2c.

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

[sysd,G] = c2d(varargin{:})
