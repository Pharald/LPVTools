function out = lpvsplit(m,varargin)
% LPVSPLIT  Extract PMAT data based on independent variable range
%
% B = LPVSPLIT(A,NAME1,RANGE1,NAME2,RANGE2,....) extracts data from the PMAT
% A on the domain specified by the NAME/RANGE pairs. Each RANGE is a 1-by-2 
% row vector [min, max], that specifies the values of the independent 
% variable to be extracted along the domain direction NAME. The data at all
% points in A.Parameter."NAME".GridData which lies inside RANGE is extracted. 
% If RANGE is a scalar then the data is extracted where the variable NAME 
% is exactly equal to RANGE. If an independent variable of A is not listed 
% in the inputs then all values along this domain direction are retained in B.
%
% B = LPVSPLIT(A,NAME1,INDEX1,NAME2,INDEX2,....,'index') extracts data from
% the PMAT A on the domain specified by the NAME/INDEX pairs. Each INDEX
% is a vector of integers or a logical array that specifies the indices
% of A.Parameter."NAME".GridData to be extracted.
%
% B = LPVSPLIT(A,NAME1,VALUES1,NAME2,VALUES2,....,'value') extracts data 
% from the PMAT A on the domain specified by the NAME/VALUES pairs. VALUES 
% specifies the A.Parameter."NAME".GridData to be extracted.
%
% B = LPVSPLIT(A,NAME,RANGE) is an alternative syntax. NAME is an N-by-1
% cell array of characters and RANGE is an N-by-1 cell array of ranges.
% B = LPVSPLIT(A,NAME,INDEX,'index') and B = LPVSPLIT(A,NAME,VALUES,'value')
% also apply for INDEX or VALUES as a cell array. 
%
% B = LPVSPLIT(A,DOMAIN) is an another alternative syntax. DOMAIN is
% an RGRID object. This extracts data from the PMAT A based on the
% independent variables and data ranges in DOMAIN.
%
% See also: lpvinterp.

% Check wether we have index or value call:
nin = nargin;
cflag = 1;
if nin>=4 && ischar(varargin{end})
    if strcmpi(varargin{end}(1),'i')
        cflag = 2;
    elseif strcmpi(varargin{end}(1),'v')
        cflag = 3;
    else
        error('Undefined option')
    end
    varargin(end) = [];
    nin = nin-1;
end

% Parse inputs
if nin==2 && isa(varargin{1},'rgrid')
    % B = LPVSPLIT(A,DOMAIN)
    D = varargin{1};
    cnames = D.IVName;
    NumIV = D.NumIV;
    newIvRange = [cellfun(@min,D.IVData), cellfun(@max,D.IVData)];
    newIvRange = mat2cell(newIvRange,ones(NumIV,1),2);
elseif nin==3 && isa(varargin{1},'cell') && isa(varargin{2},'cell')
    % B = LPVSPLIT(A,NAME,RANGE)
    cnames = varargin{1};
    newIvRange = varargin{2};
else
    % B = LPVSPLIT(A,NAME1,RANGE1,NAME2,RANGE2,....)
    cnames = varargin(1:2:end-1);
    newIvRange = varargin(2:2:end);
end


% Get PMAT data
niv = m.DomainPrivate.NumIV;
ivn = m.DomainPrivate.IVName;
ivd = m.DomainPrivate.IVData;
ivrb = m.DomainPrivate.IVRateBounds;

% Find locations of NAMEs in the PMAT list of IVs
ctmp = char([ivn;cnames(:)]);
civn = ctmp(1:niv,:);
ccnames = ctmp(niv+1:end,:);
[~,idxorig,idxnew] = intersect(civn,ccnames,'rows');

% Extract Data
if niv==0
    out = m;
else
    INDEX = cellstr(repmat(':',2+niv,1))';
    if cflag==1
        % Split functionality by range
        for i=1:length(idxnew)
            minmax = newIvRange{idxnew(i)};
            if isscalar(minmax)
                minmax = [minmax minmax];
            elseif numel(minmax)>2
                error('RANGE must be a scalar or 1-by-2 vector');
            end
            idx = find(ivd{idxorig(i)}>=minmax(1) &  ivd{idxorig(i)}<=minmax(2));
            if isempty(idx)
                out = pmat;
                return
            end
            ivd{idxorig(i)} = ivd{idxorig(i)}(idx);
            INDEX{2+idxorig(i)} = idx;
        end
    elseif cflag==2
        % Split functionality by index
        for i=1:length(idxnew)
            idx = newIvRange{idxnew(i)};
            try
                ivd{idxorig(i)} = ivd{idxorig(i)}(idx);
            catch
                error(['Index is not compatible with Domain of IV ' ...
                    ivn{idxorig(i)} ]);
            end
            INDEX{2+idxorig(i)} = idx;
        end
    else
        % Split functionality by value.
        for i=1:length(idxnew)
            [tf,idx] = ismember( newIvRange{idxnew(i)} , ivd{idxorig(i)} );
            idx = idx(tf);
            if ~all(tf)
               warning(['Input value is not in the domain grid and'...
                       ' will be ignored. Use LPVINTERP to obtain '...
                       'interpolated model at values off the domain grid.']) 
            end
            if isempty(idx)
                out = pmat;
                return
            end
            ivd{idxorig(i)} = ivd{idxorig(i)}(idx);
            INDEX{2+idxorig(i)} = idx;
        end
    end
    
    out = pmat(m.DataPrivate(INDEX{:}),rgrid(ivn,ivd,ivrb));
end




