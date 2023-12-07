function varargout = norm(varargin)
% NORM   Dynamic system norm of a PLFTSS object.
%
% NORM(S,DOMAIN) is the H2 norm of S at each point in the domain specified in the 
% RGRID object DOMAIN. NORM(S,2,DOMAIN) also returns the H2 norm.
%
% NORM(S,inf,DOMAIN) is the Hinfty norm of S at each point in the domain of S. 
% NORM(S,inf,TOL,DOMAIN) specifies the relative tolerance for the computed norm.
%
% [NINF,FPEAK] = NORM(S,inf,DOMAIN) returns PMATs NINF and FPEAK.  NINF if 
% the Hinfty norm and FPEAK the frequency at which the system achives the  
% peak for each point in the domain of S.
%
% See also: norm, sigma.


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
[varargout{:}]=norm(varargin{:});