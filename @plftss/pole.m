function P = pole(varargin)
% POLE  Pointwise computes the poles of a PLFTSS object
%
% P = POLE(SYS,DOMAIN) returns the poles P of the system SYS at each 
% point in the domain specified by the RGRID object DOMAIN.  P is 
% returned as a column vector PMAT containing the eigenvalues of the 
% state matrix A at each point in the domain of SYS.
%
% See also: pole, damp, zero.

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

P=pole(varargin{:});