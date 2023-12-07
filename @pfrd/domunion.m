function varargout = domunion(varargin)
% DOMUNION   Define PFRDs on a common domain
%
% Let A and B be PFRDs.  If A depends on independent variables (X,Y)
% and B depends on independent variables (X,Z) then
% [Aext,Bext]=domunion(A,B) returns PFRDs Aext and Bext that have
% a common domain with independent variables (X,Y,Z). Aext evaluated at
% point (x,y,z) is given by A evaluated at (x,y). Bext evaluated at
% point (x,y,z) is given by B evaluated at (x,z).
%
% Given PFRDs A1,...,AN, the syntax
%   [A1ext,...,ANext] = domunion(A1,...,AN)
% constructs A1ext,...,ANext that are defined on a common domain.


% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(2, inf)
nargoutchk(0, nin+1)

if islogical(varargin{end})
    flg = varargin{end};
    nin = nin-1;
else
    flg = false;
end

if ~flg
   
   % Get domains for each PFRD
   dcell = cell(nin,1);
   for i=1:nin
      varargin{i} = pfrd( varargin{i} );
      dcell{i} = varargin{i}.DomainPrivate;
   end
   
   % Construct single domain containing union of IVs in input domains
   idxcell = cell(nin,1);
   [Udom,idxcell{:}] = domunion( dcell{:} );
   szdom = size(Udom);
   
   % Expand each input PFRD to be constant along new IV dimensions
   varargout = cell(nin,1);
   for i=1:nin
      A = varargin{i};
      Aidx = idxcell{i};
      if numel(Aidx)==0
         Adata = A.DataPrivate;
      elseif numel(Aidx)==1
         Adata = permute(A.DataPrivate,[Aidx 2]);
      else
         Adata = permute(A.DataPrivate,[Aidx]);
      end
      repval = ones(1,length(szdom));
      idx = Aidx > A.DomainPrivate.NumIV;
      repval(idx) = szdom(idx);
      Adata = repsys(Adata,[1 1 repval]);
      Adata = adscalarexp(Adata,Udom);
      Aext = pfrd(Adata,Udom);
      varargout{i} = Aext;
   end
else
   % Get public domains for each PFRD
   dcell = cell(nin,1);
   for i=1:nin
      varargin{i} = pfrd( varargin{i} );
      dcell{i} = varargin{i}.Domain;
   end
   
   % Construct single domain containing union of IVs in input domains
   idxcell = cell(nin,1);
   [Udom,idxcell{:}] = domunion( dcell{:} );
   szdom = size(Udom);
   
   % Expand each inM1put PFRD to be constant along new IV dimensions
   varargout = cell(nin,1);
   for i=1:nin
      A = varargin{i};
      Aidx = idxcell{i};
      niv = A.Domain.NumIV;
      nad = numel(size(A))-2;
      
      % Reorder as [row col AD IV]
      if nad+niv>=2
         Adata = permute(A.Data,[(1+niv:niv+nad) (1:niv)]);
         % Reorder IVs to align with common domain
         Adata = permute(Adata,[1:nad nad+Aidx]);
      else
         Adata = A.Data;
      end
      
      % Fan out singleton IVs to proper dimension
      repval = ones(1,length(szdom));
      idx = Aidx > niv;
      repval(idx) = szdom(idx);
      Adata = repmat(Adata,[ones(1,nad) repval]);
      
      % Reorder as [row col IV AD]
      nuv = Udom.NumIV;
      if nad+nuv>=2
         Adata = permute(Adata,[(1+nad:nad+nuv)  (1:nad)]);
      end
      
      % Pack up as PFRD
      Adata = adscalarexp(Adata,Udom);
      Aext = pfrd(Adata,Udom);
      varargout{i} = Aext;
   end
end
