function sysc = d2c(sysd,varargin)
% D2C   Discrete-time to continuous-time conversion for PSS objects
%
% SYSC = D2C(SYSD,METHOD) computes a continuous-time PSS SYSC that
% approximates the discrete-time PSS SYSD.  The conversion is performed at
% each point in the domain of SYSD.  The string METHOD selects the
% discretization method among 'zoh' (default), 'tustin', and 'matched'.
%
% D2C(SYSC,OPTIONS) gives access to additional discretization options.
% Use D2COPTIONS to create and configure the option set OPTIONS.
%
% See also:  d2c, d2cOptions, c2d.

% TODO PJS 5/18/2011: In some cases d2c on an Data can return models
% of varying state dimension (e.g. a discrete system with created with
% DRSS can have negative poles and d2c will increase the model order to
% handle these poles).  PSS objects are not allowed to have varying state
% dimension and this function will error out at the PSS construction.
% We should probably add a more intelligent error message if this occurs.

% Check inputs/outputs
nin = nargin;
narginchk(1, 2);
nout = nargout;
nargoutchk(0, 1);

% Conversion
% XXX PJS TODO 1/24/2014: Add functionality for second output argument?
D = sysd.DomainPrivate;
sysc = d2c(sysd.DataPrivate,varargin{:});
sysc = pss(sysc,D);

