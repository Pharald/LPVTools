function B = reshape(A,varargin)
% RESHAPE   Reshape PSTRUCT arrays
%
% B = RESHAPE(A,M,N) reshapes A into an M-by-N PSTRUCT. The elements of A are 
% taken columnwise from A.  A must have M*N elements. B = RESHAPE(A,[M N]) 
% is an alternative syntax
%
% B = RESHAPE(A,M,N,P,..,Q) reshapes A into the M-by-N-by-P-by-...-by-Q
% array B. M*N*P*...*Q must equal PROD(SIZE(A)).
%
% B = RESHAPE(A, ..., [], ... ) leaves one desired dimension unspecified, 
% and replaces the corresponding entry with []. This dimension is computed 
% automatically so that the product of the desired dimensions matches PROD(SIZE(A)).
% Only one occurance of [] can be used.
%
% See also: squeeze, shiftdim, colon.


% Reorder as [row col AD IV]
niv = A.Domain.NumIV;
nad = numel(size(A))-2;
liv = A.Domain.LIVData;
Adata = permute(A.Data,[1 2 (3+niv:2+niv+nad) (3:2+niv)]);

% Convert single vector reshape dimension to cell array list 
if nargin==1
   error('Not enough input arguments');
elseif nargin==2 
   if isnumeric(varargin{1})
      rsdim = varargin{1}(:);
      rsdim = num2cell( rsdim' );
   else
      error('Size vector must be a row vector with integer elements.');
   end
else
   rsdim = varargin;
end

% Append list of IV dimensions to the reshape dimension list
alldim = [ rsdim   num2cell(liv)'];

% Do the reshape
Adata = reshape( Adata, alldim{:});

% Reorder as [row col IV AD]
if length(rsdim)==1
   Adata = permute(Adata,[1 niv+2 2:niv+1]);   
else
   nad = length(rsdim)-2;
   Adata = permute(Adata,[1 2 (3+nad:2+nad+niv)  (3:2+nad)]);
end
            
% Pack up result
B = pstruct(Adata,A.Domain);

