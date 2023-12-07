function varargout = dare(varargin)
% DARE  Pointwise solve discrete-time algebraic Riccati equations for PLFTMATs
%
% [X,L,G,REPORT] = DARE(A,B,Q,R,S,E,DOMAIN) computes the unique 
% stabilizing solution X of the discrete-time algebraic Riccati equation 
%        E'XE = A'XA - (A'XB + S)*inv(B'XB + R)*(A'XB + S)' + Q 
% at each point in the domain specified by the RGRID object DOMAIN.
% X is returned as a PMAT that specifies the solution at each point in 
% DOMAIN. Similarly, G is a PMAT that specifies the gain matrix, 
% L is a PMAT that contains the eigenvalues,  and REPORT is a PSTRUCT that 
% contains the diagnosis reports. See numerics/dare for additional calling 
% options.
%
% See also: dare, care.


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
[varargout{:}]=dare(varargin{:});