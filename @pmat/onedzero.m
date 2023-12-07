function zz = onedzero(mat,imeth)
% finds zeros of scalar 1-d function

% TODO PJS 4/4/2011: Revisit

if nargin==1
   imeth = 'linear';
end

% Restrict function to no array dimensions
if hasArray(mat)
   error('ONEDZERO is restricted to PMATs without array dimensions')
end

ivd = mat.DomainPrivate.IVData;
niv = mat.DomainPrivate.NumIV;
szm = [privatesize(mat) 1];
if szm(1)==1 && szm(2)==1 && niv==1 
   optim = optimset('display','off');
   zz = findzero(mat.DataPrivate(:),ivd{1});
   for i=1:length(zz)
      zz(i) = fminsearch(@absevalpmat,zz(i),optim,mat.DataPrivate(:),ivd{1},imeth);
   end
else
   error('matrix should be scalar, and 1 IV')
end

function zz = findzero(val,ivd)
val = val(:);
ivd = ivd(:);
loc = find(val(1:end-1).*val(2:end)<=0);
zz = zeros(1,length(loc));
for k=1:length(loc)
   kk = loc(k);
   if val(kk)==0
      %zz = [zz val(kk)];
      zz(k) = val(kk);
   else
      r = ivd(kk) + val(kk)*(ivd(kk)-ivd(kk+1))/(val(kk+1)-val(kk));
      %zz = [zz r];
      zz(k) = r;
   end
end

function val = absevalpmat(x,vals,ivd,imeth)

if x<min(ivd) || x>max(ivd)
  val = 5+10*max(abs(vals));
else
  val = abs(interp1(ivd,vals,x,imeth));
end
