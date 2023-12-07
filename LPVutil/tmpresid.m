function Pout = tmpresid(P,poleloc)
%% TODO allow 3rd "cutoff" sv tolerance, then no user interaction
%% TODO better name for function
[fesp,edim] = chkeigsp(P,poleloc);
if ~isempty(fesp)
   szF = size(fesp);
   % "normally", there are more samples than states, but not always, for
   % instance if there are just a few samples.  This impacts how we use the
   % economy SVD
   if szF(1)<=szF(2)
      % [U,S] = svd(fesp,'econ');
      [U,S] = svd([imag(fesp) real(fesp)],'econ');
   else
      [U,S] = svd(fesp);
   end
   dS = diag(S);
   disp(['Singular values of Fast (' int2str(edim) ') Eigenspace across parameter'])
   disp([(1:length(dS))' dS(:)])
   % TODO loop for valid entry
   dim = input('Enter dimension of fast eigenspace: ');
   if dim>0
      % Changed from P to P.DataPrivate to address error 6 May 13
      tmp = ss2ss(P.DataPrivate,U');
      sztmp = size(tmp);     
      ztmp = tmp;
      for j=1:prod(sztmp(3:end))
          ztmp(:,:,j) = modred(tmp(:,:,j),1:dim);
      end
      Pout = pss(ztmp,P.DomainPrivate);
   else
      Pout = P;
   end
else
   % TODO revisit these message/output choices
   Pout = P;
   disp('Nothing to residulize');
end

