function R = ndLinInterp(M,nIV,N,IVData,dIVData,pValue)
% M: nR-by-nC-by-[IVDims] DOUBLE array
% nIV: scalar, double, number of IVs
% N: nIV-by-1, length of each IVGrid
% IVData: nIV-by-1 CELL, IVData{k} has k'th IVGrid, N(k)-by-1
% dIVData: nIV-by-1 CELL, dIVData{k} has diff(k'th IVGrid), (N(k)-1)-by-1
% pValue: nIV-by-1 DOUBLE, parameter value at which interpolation should
% take place.

szM = size(M);
nVertices = 2^nIV;
factor = ones(nVertices,1);
idxcell = cell(1,nIV);
for k=1:nIV
   [i,alpha] = LOCALfindslotalpha(N(k),IVData{k},pValue(k),dIVData{k});
   if ~isnan(i)
      if N(k)==1
         idxcell{k} = [1 1];
      else
         if i==N(k)
            i = N(k)-1;
            alpha = 1;
         end
         idxcell{k} = [i i+1];
      end
      vecOnes = ones(2^(k-1),1);
      
      % factor = factor.*kron(ones(2^(nIV-k),1),[(1-alpha)*vecOnes;alpha*vecOnes]);
      ma = 2^(nIV-k);
      mb = 2^k;
      yy = (1:mb).';
      ib = yy(:,ones(1, ma));
      oma = 1-alpha;
      B = [oma(vecOnes); alpha(vecOnes)];
      K = B(ib,1);
      factor = factor.*K;
   else
       error('Parameter value is outside of the domain.')
   end
   
       
end
R = reshape(reshape(M(:,:,idxcell{:}),[szM(1)*szM(2) nVertices])*factor,[szM(1) szM(2)]);


function [i,alpha] = LOCALfindslotalpha(N,vec,val,dvec)
% N integer
% vec 1-by-N (or N-by-1), sorted
% val, scalar, vec(1) <= val <= vec(N)
% dvec = diff(vec)

i = max(find(val>=vec));   % don't follow advice - it is slower.
if ~isempty(i)
   if i<N
      alpha = (val - vec(i))/dvec(i);
   elseif val==vec(N)
      alpha = 0;
   else
      i = NaN;
      alpha = NaN;
   end
else
   i = NaN;
   alpha = NaN;
end

   
