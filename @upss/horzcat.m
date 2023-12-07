function out = horzcat(varargin)
% HORZCAT   Horizontal concatenation of UPSS objects
%
% SYS = HORZCAT(SYS1,SYS2,...) performs a concatenation operation of
% SYS = [SYS1 , SYS2, ...] at each point in the combined domains
% of SYS1, SYS2, ...
%
% See also: horzcat, vertcat.

% Check # of input arguments
narginchk(1, inf)

if nargin==1
    out = varargin{1};
else
    varargin{1} = upss(varargin{1});
    out = binop(varargin{1},varargin{2},'horzcat');
    if nargin>2
        out = horzcat(out,varargin{3:end});
    end    
end


% ----- OLD CODE

% function out = horzcat(varargin)
% 
% if nargin==1
%    out = varargin{1};
% elseif nargin==2
%    if isa(varargin{1},'pss') && isa(varargin{2},'pss')
%       A = varargin{1};
%       sa = size(A);
% 	   B = varargin{2};
%       sb = size(B);
%       if sa(1)==sb(1)
%          if sa(2)==0
%             out = B;
%          elseif sb(2)==0
%             out = A;
%          else
%             try
%                out = binop(A,B,'horzcat');
%             catch
%                error(lasterr)
%             end
%          end
%       else
%          error('Invalid Row dimensions');
%       end
%    else
%       out = horzcat(pss(varargin{1}),pss(varargin{2}));
%    end
% else
%    out = horzcat(horzcat(pss(varargin{1}),pss(varargin{2})),varargin{3:end});
% end
% 
% 
% 
