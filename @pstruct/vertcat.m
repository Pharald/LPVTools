function out = vertcat(varargin)
% VERTCAT  Vertical concatenation of PSTRUCT objects.
%
% S = VERTCAT(S1,S2,...) performs a concatenation operation of
% S = [S1; S2; ...]  at each point in the combined domains of S1, S2, ...
%
% See also: vertcat, horzcat.

% Check # of input arguments
narginchk(1, inf)

if nargin==1
    out = varargin{1};
else
    varargin{1} = pstruct(varargin{1});
    out = binop(varargin{1},varargin{2},'vertcat');
    if nargin>2
        out = vertcat(out,varargin{3:end});
    end    
end


