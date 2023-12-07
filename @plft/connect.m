function out = connect(varargin)
%CONNECT  Block-diagram interconnections of PLFT
%
% CONNECT overloads the standard MATLAB function to work with PLFT objects.
% See DynamicSystem/connect for details.
%
% See also: connect, append, series, parallel, sysic, feedback, lft.

RateBounds = cell(0,2);
for i=1:numel(varargin)
    vi = varargin{i};
    if isa(vi,'plft')
        RateBounds = [RateBounds;vi.RateBounds];
        varargin{i} = vi.Data;
    end
end

% Apply CONNECT
out = connect(varargin{:});


% Find unique TVREAL list 
[~,idx1,idx2] = unique(RateBounds(:,1));
RBUnique = RateBounds(idx1,:);

% Check for conflicting rate bound info
RBInterval1 = cell2mat( RateBounds(:,2) );
RBInterval2 = cell2mat( RBUnique(idx2,2) );
[i,j] = find( RBInterval1~=RBInterval2 );
if ~isempty(i)
    Name = RateBounds{i(1),1};
    error(['The block with name "' Name '"  has two or more' ...
        ' conflicting rate bound definitions.']);
end

% Remove RB info if corresponding tvreal does not appear in the result
if isuncertain(out)    
    if isa(out,'umat') || isa(out,'uss') || ...
            isa(out,'ufrd')
        U = fieldnames(out.Uncertainty);
    else
        % XXX review this: what is U.Name?
        U = out.Name;
    end        
else
    U = [];
end
[~,idx]=setdiff(RBUnique(:,1),U);
RBUnique(idx,:)=[];

% Package data in the proper object
switch class(out)
    case {'double','umat'}
        out = plftmat(umat(out),RBUnique);
    case {'ss','tf','zpk','uss'}
        out = plftss(uss(out),RBUnique);
    case {'frd','ufrd'}
        out = plftfrd(ufrd(out),RBUnique);
    otherwise
        % XXX PJS: Update error to provide more info.
        error('Unsupported binary operation.');
end
        

