function [sysd,G] = c2d(sysc,Ts,varargin)
% C2D   Continuous-time to discrete-time conversion for PSS objects
%
% SYSD = C2D(SYSC,TS,METHOD) computes a discrete-time PSS SYSD with 
% sampling time TS that approximates the continuous-time PSS SYSC.  The
% conversion is performed at each point in the domain of SYSC.  The string 
% METHOD selects the discretization method among 'zoh' (default), 'foh',
% 'impulse', 'tustin', and 'matched'.
%
% C2D(SYSC,TS,OPTIONS) gives access to additional discretization options. 
% Use C2DOPTIONS to create and configure the option set OPTIONS.
%
% For state-space models without delays, [SYSD,G] = C2D(SYSC,Ts,METHOD), 
% also returns the matrix G mapping the states xc(t) of SYSC to the states 
% xd[k] of SYSD: xd[k] = G * [xc(k*Ts) ; u[k]]  
%
% See also:  c2d, c2dOptions, d2c.

% Check inputs/outputs
nin = nargin;
narginchk(2, 3);
nout = nargout;
nargoutchk(0, 2);

% Conversion
D = sysc.DomainPrivate;
if nout==2
    [sysd,G] = c2d(sysc.DataPrivate,Ts,varargin{:});
    G = pmat(G,D);
    sysd = pss(sysd,D);
else
    sysd = c2d(sysc.DataPrivate,Ts,varargin{:});
    sysd = pss(sysd,D);
end
