function g = frd(varargin)
% FRD  Point-wise frequency response for PLFTSS
%
% G = frd(SYS,FREQS,DOMAIN) converts the PLFTSS SYS into a PFRD frequency 
%   response data model G, defined at each point in the domain specified by 
%   the RGRID object DOMAIN. FREQS is a Nf-by-1 DOUBLE of frequency values. 
%   G is returned as a PFRD containing the frequency response of SYS at 
%   each point in FREQS.
%
% See also: frd, pfrd.


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

g = frd(varargin{:});