function out = horzcat(varargin)
% HORZCAT   Horizontal concatenation of PFRD objects
%
% S = HORZCAT(S1,S2,...) performs a concatenation operation of
% S = [S1 , S2, ...] at each point in the combined domains of S1, S2, ...
%
% See also: horzcat, vertcat.

% Check # of input arguments
narginchk(1, inf)

if nargin==1
    out = varargin{1};
else
    varargin{1} = pfrd(varargin{1});
    out = binop(varargin{1},varargin{2},'horzcat');
    if nargin>2
        out = horzcat(out,varargin{3:end});
    end    
end

