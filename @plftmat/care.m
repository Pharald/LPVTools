function varargout = care(varargin)
% CARE  Pointwise solve continuous-time algebraic Riccati equations for PLFTMATs
%
% [X,L,G,REPORT] = CARE(A,B,Q,R,S,E,DOMAIN) computes the unique 
% stabilizing solution X of the continuous-time algebraic Riccati equation 
%        A'XE + E'XA - (E'XB + S) inv(R)  (B'XE + S') + Q = 0 .
% at each point in the domain specified by the RGRID object DOMAIN.
% X is returned as a PMAT that specifies the solution at each point in the 
% combined domain. Similarly, G is a PMAT that specifies the gain matrix, 
% L is a PMAT that contains the eigenvalues,  and REPORT is a PSTRUCT that 
% contains the diagnosis reports. See numerics/care for additional calling 
% options.
%
% See also: care, dare.


DOM = varargin{end};
varargin(end) = [];

if ~isa(DOM,'rgrid')
    error(['Last input argument must be an RGRID object that specifies'...
           ' the parameter domain grid.'])
end

for i = 1:numel(varargin)
    if isa(varargin{i},'plftss') || isa(varargin{i},'plftmat')
        varargin{i} = lft2grid(varargin{i},DOM);
    end
end

varargout = cell(nargout,1);
[varargout{:}]=care(varargin{:});