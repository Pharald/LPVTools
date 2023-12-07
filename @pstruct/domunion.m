function varargout = domunion(varargin)
% DOMUNION   Define PSTRUCTs on a common domain
%
% Let A and B be PSTRUCTs.  If A depends on independent variables (X,Y)
% and B depends on independent variables (X,Z) then
% [Aext,Bext]=domunion(A,B) returns PSTRUCTs Aext and Bext that have
% a common domain with independent variables (X,Y,Z). Aext evaluated at
% point (x,y,z) is given by A evaluated at (x,y). Bext evaluated at
% point (x,y,z) is given by B evaluated at (x,z).
%
% Given PSTRUCTs A1,...,AN, the syntax
%   [A1ext,...,ANext] = domunion(A1,...,AN)
% constructs A1ext,...,ANext that are defined on a common domain.


% If last argument is 0 or 1 then:
% 0:= Work as documented
% 1:= Only find union of IV domains and ignore the array dimensions

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(2, inf)
nargoutchk(0, nin+1)

if isa(varargin{end},'pstruct')
   flg = 0;
else
   flg = varargin{end};
   nin = nin-1;
end

if ~flg
   % Get domains for each PSTRUCT
   dcell = cell(nin,1);
   for i=1:nin
      varargin{i} = pstruct( varargin{i} );
      dcell{i} = varargin{i}.DomainPrivate;
   end
   
   % Construct single domain containing union of IVs in input domains
   idxcell = cell(nin,1);
   [Udom,idxcell{:}] = domunion( dcell{:} );
   szdom = size(Udom);
   
   % Expand each input PSTRUCT to be constant along new IV dimensions
   varargout = cell(nin,1);
   for i=1:nin
      A = varargin{i};
      Aidx = idxcell{i};
      
      Adata = permute(A.DataPrivate,[1 2 2+Aidx]);
      repval = ones(1,length(szdom));
      idx = Aidx > A.DomainPrivate.NumIV;
      repval(idx) = szdom(idx);
      Adata = repmat(Adata,[1 1 repval]);
      Adata = adscalarexp(Adata,Udom);
      Aext = pstruct(Adata,Udom);
      varargout{i} = Aext;
   end
else
   % Get public domains for each PSTRUCT
   dcell = cell(nin,1);
   for i=1:nin
      varargin{i} = pstruct( varargin{i} );
      dcell{i} = varargin{i}.Domain;
   end
   
   % Construct single domain containing union of IVs in input domains
   idxcell = cell(nin,1);
   [Udom,idxcell{:}] = domunion( dcell{:} );
   szdom = size(Udom);
   
   % Expand each input PSTRUCT to be constant along new IV dimensions
   varargout = cell(nin,1);
   for i=1:nin
      A = varargin{i};
      Aidx = idxcell{i};
      niv = A.Domain.NumIV;
      nad = numel(size(A))-2;
            
      % Reorder as [row col AD IV]
      Adata = permute(A.Data,[1 2 (3+niv:2+niv+nad) (3:2+niv)]);
      
      % Reorder IVs to align with common domain
      Adata = permute(Adata,[1:2+nad 2+nad+Aidx]);
      
      % Fan out singleton IVs to proper dimension
      repval = ones(1,length(szdom));
      idx = Aidx > niv;
      repval(idx) = szdom(idx);
      Adata = repmat(Adata,[ones(1,2+nad) repval]);
      
      % Reorder as [row col IV AD]
      nuv = Udom.NumIV;
      Adata = permute(Adata,[1 2 (3+nad:2+nad+nuv)  (3:2+nad)]);
            
      % Pack up as PSTRUCT
      Adata = adscalarexp(Adata,Udom);
      Aext = pstruct(Adata,Udom);
      varargout{i} = Aext;
   end
end


