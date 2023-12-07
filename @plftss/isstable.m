function s = isstable(varargin)
% ISSTABLE  Pointwise stability test for a PLFTSS object.
%
% ISSTABLE(SYS) tests the stability of SYS at each point in the domain 
% specified by the RGRID object DOMAIN. Returns a PMAT on DOMAIN, 
% containing True at those points where SYS is stable, and False otherwise.
%  
% See also: isstable, pole.

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

s = isstable(varargin{:});