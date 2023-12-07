function B = usubs(A,varargin)
%  usample   Generates random samples of a plft.

A = plft(A);
ARateBounds = A.RateBounds;
AData = A.Data;

for i =1:length(varargin)
    if isa(varargin{i},'plft')                
        ARateBounds = [ARateBounds;varargin{i}.RateBounds];
        varargin{i} = varargin{i}.Data;
    end
end
BData = usubs(AData,varargin{:});    


% Find unique TVREAL list 
[~,idx1,idx2] = unique(ARateBounds(:,1));
RBUnique = ARateBounds(idx1,:);

% Check for conflicting rate bound info
RBInterval1 = cell2mat( ARateBounds(:,2) );
RBInterval2 = cell2mat( RBUnique(idx2,2) );
[i,j] = find( RBInterval1~=RBInterval2 );
if ~isempty(i)
    Name = ARateBounds{i(1),1};
    error(['The block with name "' Name '"  has two or more' ...
        ' conflicting rate bound definitions.']);
end



% Remove RB info if corresponding tvreal does not appear in the result
if isuncertain(BData)      
    if isa(BData,'umat') || isa(BData,'uss') || ...
            isa(BData,'ufrd')
        U = fieldnames(BData.Uncertainty);
    else
        U = BData.Name;
    end
else
    U = [];
end
[~,idx]=setdiff(RBUnique(:,1),U);
RBUnique(idx,:)=[];

% Package data in the proper object
switch class(BData)
    case {'double','ss','frd'}
        B = BData;
    case 'umat'
        B = plftmat(BData,RBUnique);
    case 'uss'
        B = plftss(BData,RBUnique);
    case 'ufrd'
        B = plftfrd(BData,RBUnique);
    otherwise
        % XXX PJS: Update error to provide more info.
        error('Unsupported operation.');
end

