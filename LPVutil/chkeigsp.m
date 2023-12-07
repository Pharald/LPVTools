function [fesp,edim] = chkeigsp(P,evlcu)
%
% Find a common-dimension eigenspace for all eigenvalues of a PSS system
% whose real part is more negative than a given value.
% TODO what portion belongs in @PSS

szP = size(P.Data);
fesp = [];
for i=1:prod(szP(3:end))
   [evc,evl] = eig(P.Data(:,:,i).a);
   bigneg = find(diag(real(evl))<-abs(evlcu));
               length(bigneg)  % JT
   if i==1
      edim = length(bigneg);
   else
      if edim~=length(bigneg)
         fesp = zeros(0,0);
         edim = 0;
         return
      end
   end
   fesp = [fesp evc(:,bigneg)];
end
