function out = vertcat(varargin)
% VERTCAT  Vertical concatenation of PMAT objects.
%
% MAT = VERTCAT(MAT1,MAT2,...) performs a concatenation operation of
% MAT = [MAT1; MAT2; ...]  at each point in the combined domains
% of MAT1, MAT2, ...
%
% See also: vertcat, horzcat.

% Check # of input arguments
narginchk(1, inf)

if nargin==1
    out = varargin{1};
else
    varargin{1} = pmat(varargin{1});
    out = binop(varargin{1},varargin{2},'vertcat');
    if nargin>2
        out = vertcat(out,varargin{3:end});
    end    
end



% ----- OLD CODE

% if nargin==1
%     out = varargin{1};
% else
%     if isa(varargin{1},'pmat') && isa(varargin{2},'pmat')
%         sza = size(varargin{1});
%         szb = size(varargin{2});
%         if sza(1)==0
%             out = varargin{2};
%             if sza(2)>0 && sza(2)~=szb(2)
%                 warning('Concatenation involves an empty array with incorrect number of columns.');
%             end
%         elseif szb(1)==0
%             out = varargin{1};
%             if szb(2)>0 && szb(2)~=sza(2)
%                 warning('Concatenation involves an empty array with incorrect number of columns.');
%             end
%         elseif sza(2)==szb(2)
%             out = binop(varargin{1},varargin{2},'vertcat');
%         else
%             error('pmat Column dimensions are not equal')
%         end
%         %out = genredu(out);
%     else
%         out = vertcat(pmat(varargin{1}),pmat(varargin{2}));
%     end
%     
%     if nargin>2
%         out = vertcat(out,varargin{3:end});
%     end
% end

