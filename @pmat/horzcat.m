function out = horzcat(varargin)
% HORZCAT   Horizontal concatenation of PMAT objects
%
% MAT = HORZCAT(MAT1,MAT2,...) performs a concatenation operation of
% MAT = [MAT1 , MAT2, ...] at each point in the combined domains
% of MAT1, MAT2, ...
%
% See also: horzcat, vertcat.

% TODO PJS 4/1/2011: CAT is currently not implemented. Does it make
% sense to implement this when PMATs are viewed as 2-dim matrix with IVs?

% Check # of input arguments
narginchk(1, inf)

if nargin==1
    out = varargin{1};
else
    varargin{1} = pmat(varargin{1});
    out = binop(varargin{1},varargin{2},'horzcat');
    if nargin>2
        out = horzcat(out,varargin{3:end});
    end    
end



% ----- OLD CODE

% if nargin==1
%     out = varargin{1};
% else
%     if isa(varargin{1},'pmat') && isa(varargin{2},'pmat')
%         sza = size(varargin{1});
%         szb = size(varargin{2});
%         if sza(2)==0
%             out = varargin{2};
%             if sza(1)>0 && sza(1)~=szb(1)
%                 warning('Concatenation involves an empty array with incorrect number of rows.');
%             end
%         elseif szb(2)==0
%             out = varargin{1};
%             if szb(1)>0 && szb(1)~=sza(1)
%                 warning('Concatenation involves an empty array with incorrect number of rows.');
%             end
%         elseif sza(1)==szb(1)
%             out = binop(varargin{1},varargin{2},'horzcat');
%         else
%             error('pmat Row dimensions are not equal')
%         end
%         %out = genredu(out);
%     else
%         out = horzcat(pmat(varargin{1}),pmat(varargin{2}));
%     end
%     
%     if nargin>2
%         out = horzcat(out,varargin{3:end});
%     end    
% end
