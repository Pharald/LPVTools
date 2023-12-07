function varargout = qr(varargin)
% QR    Orthogonal-triangular decomposition for PLFTMAT objects.
%
% [Q,R]=QR(A,DOMAIN) computes the QR decomposition of A at each point in 
% the domain specified by the RGRID object DOMAIN. See additional calling 
% options in the QR help.
%
% See also: qr.

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
[varargout{:}]=qr(varargin{:});