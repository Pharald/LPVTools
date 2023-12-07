function C = binop(A,B,op,varargin)
% BINOP Binary operations for a PLFT
% 
% BINOP is a utility function that handles binary operations (e.g. times,
% ldivide, lft, etc.) on PLFTs. The user should not call BINOP directly. 
% Instead each binary operation has a front end function that the user 
% calls (e.g. *,+,lft,horzcat). The frontend function then call calls BINOP 
% to do the required calculations/operations

%% PLFT/BINOP

% NOTE PJS 4/7/2011:
% FEVAL is a bit slower than directly performing the operation, eg.
%   if M is a matrix then M*M is about 4x faster than
%   feval('mtimes',M,M) and M*M is about 2x faster than
%   fh(M,M) where fh=@mtimes.
% The main BINOPs are implemented directly, e.g. plus uses '+'.

% Grab all rate bound information
if isa(A,'plft')    
    ARateBounds = A.RateBounds;
    AData = A.Data;
else
    ARateBounds = cell(0,2);
    AData = A;
end

if isa(B,'plft')  
    BRateBounds = B.RateBounds;
    BData = B.Data;
else
    BRateBounds = cell(0,2);
    BData = B;
end

% Perform binary operation on umat/uss/ufrd data
switch op
    case 'subsasgn'
        % Single step subsasgn (L is 1-by-1) on Data: 
        %    C = binop(A,B,'subsasgn',L);
        L = varargin{1};
        CData = subsasgn(AData,L,BData);        
    case 'horzcat'
        CData = [AData BData];
    case 'vertcat'
        CData = [AData; BData];
    case 'plus';
        CData = AData+BData;
    case 'minus';
        CData = AData-BData;
    case 'mtimes'
        CData = AData*BData;
%     case 'times'
%         CData = AData.*BData;
    case 'mpower'
        CData = AData^BData;
%     case 'power'
%         CData = AData.^BData;
    case 'mldivide'
        CData = AData\BData;
    case 'mrdivide'
        CData = AData/BData;
%     case 'ldivide'
%         CData = AData.\BData;
%     case 'rdivide'
%         CData = AData./BData;
    case 'blkdiag'
        CData = blkdiag(AData,BData);
    case 'append'
        CData = append(AData,BData);
    case 'stack'
        CData = stack(varargin{1},AData,BData);
    case 'starp'
        % TODO PJS: Is this obsolete if we implement lft binop?
        CData = starp(AData,BData,varargin{:});
    case 'lft'
        CData = lft(AData,BData,varargin{:});
    case 'feedback'
        CData = feedback(AData,BData,varargin{:});
    case 'series'
        CData = series(AData,BData,varargin{:});
    case 'parallel'
        CData = parallel(AData,BData,varargin{:});
    case 'eq'
        CData = eq(AData,BData);
    case 'ne'
        CData = ne(AData,BData);
    case 'lt'
        CData = lt(AData,BData);
    case 'gt'
        CData = gt(AData,BData);
    case 'le'
        CData = le(AData,BData);
    case 'ge'
        CData = ge(AData,BData);
    case 'and'
        CData = and(AData,BData);
    case 'or'
        CData = or(AData,BData);
    case 'xor'
        CData = xor(AData,BData);
    otherwise
        % TODO PJS Put this in a try-catch?
        funchan = eval(['@' op]);
        CData = funchan(AData, BData , varargin{:} );
end

% Find unique TVREAL list 
CRateBounds = [ARateBounds; BRateBounds];
[~,idx1,idx2] = unique(CRateBounds(:,1));
RBUnique = CRateBounds(idx1,:);

% Check for conflicting rate bound info
RBInterval1 = cell2mat( CRateBounds(:,2) );
RBInterval2 = cell2mat( RBUnique(idx2,2) );
[i,j] = find( RBInterval1~=RBInterval2 );
if ~isempty(i)
    Name = CRateBounds{i(1),1};
    error(['The block with name "' Name '"  has two or more' ...
        ' conflicting rate bound definitions.']);
end

% Remove RB info if corresponding tvreal does not appear in the result
if isuncertain(CData)    
    if isa(CData,'umat') || isa(CData,'uss') || ...
            isa(CData,'ufrd')
        U = fieldnames(CData.Uncertainty);
    else
        U = CData.Name;
    end        
else
    U = [];
end
[~,idx]=setdiff(RBUnique(:,1),U);
RBUnique(idx,:)=[];

% Package data in the proper object
switch class(CData)
    case {'double','umat'}
        C = plftmat(umat(CData),RBUnique);
    case {'ss','tf','zpk','uss'}
        C = plftss(uss(CData),RBUnique);
    case {'frd','ufrd'}
        C = plftfrd(ufrd(CData),RBUnique);
    otherwise
        % XXX PJS: Update error to provide more info.
        error('Unsupported binary operation.');
end

