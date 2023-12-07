function out = vertcat(varargin)
% VERTCAT  Vertical concatenation of PFRD objects.
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
    varargin{1} = pfrd(varargin{1});
    out = binop(varargin{1},varargin{2},'vertcat');
    if nargin>2
        out = vertcat(out,varargin{3:end});
    end    
end



% ----- OLD CODE

% function out = vertcat(varargin)
% 
% if nargin==1
%    out = varargin{1};
% elseif nargin==2
%    if isa(varargin{1},'pss') && isa(varargin{2},'pss')
%       A = varargin{1};
%       sza = size(A);
% 	   B = varargin{2};
%       szb = size(B);
%       if sza(1)==0 && (sza(2)==szb(2) || sza(2)==0)
%          out = B;
%       elseif szb(1)==0 && (sza(2)==szb(2) || szb(2)==0)
%          out = A;
%       elseif sza(2)==szb(2)
%          try
%             out = binop(A,B,'vertcat');
%          catch
%             error(lasterr)
%          end
%       else
%          error('Incompatible Column Dimensions')
%       end      
%    else
%       out = vertcat(pss(varargin{1}),pss(varargin{2}));
%    end
% else
%    out = vertcat(vertcat(pss(varargin{1}),pss(varargin{2})),varargin{3:end});
% end
% 


