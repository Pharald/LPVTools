function [Ke,Gamma,Info] = ...
   lpvof2mi(sys,nmeas,ncont,BOF,...
   Xb,Yb,derUpperBound,derLowerBound)
% 
% derUpperBound, (<=niv)-by-2, cell array, {PMAT/DOUBLE variableName}
% derLowerBound, (<=niv)-by-2, cell array, {PMAT/DOUBLE variableName}

nin = nargin;
if nin<3
   error('correct syntax requires (Sys,Nmeas,NControl)');
elseif nin==3
   BOF = 1.01;
elseif nin==4
   Xb(1,1).BasisFunction = pmat(1);
   Xb(1,1).Partial = cell(0,2);
   Yb = Xb;
   derUpperBound = cell(0,2);
   derLowerBound = cell(0,2);
elseif nin~=8
   error('Must specify Basis and Parameter Rate bounds')
end

[P,TL,TR,FT,ne,nd,~,~] = orthog4syn(sys,nmeas,ncont);

% Get state space data
[a,b,c,d] = ssdata(P);
nX = size(a,1);
ne1 = ne-ncont;
ne2 = ncont;
nd1 = nd-nmeas;
nd2 = nmeas;
d11 = d(1:ne,1:nd);
d11dot1 = d(1:ne,1:nd1);
d11dot2 = d(1:ne,nd1+1:nd);
d111dot = d(1:ne1,1:nd);
d112dot = d(ne1+1:ne,1:nd);
d1111 = d(1:ne1,1:nd1);
d1112 = d(1:ne1,nd1+1:nd);
d1121 = d(ne1+1:ne,1:nd1);
d1122 = d(ne1+1:ne,nd1+1:nd);

b11 = b(:,1:nd1);
b12 = b(:,nd1+1:nd);
b1 = [b11 b12];
b2 = b(:,nd+1:end);
c11 = c(1:ne1,:);
c12 = c(ne1+1:ne,:);
c1 = [c11;c12];
c2 = c(ne+1:end,:);
Ahat = a - b2*c12;
Bhat = b1 - b2*d112dot;
Atil = a - b12*c2;
Ctil = c1 - d11dot2*c2;

Pivn = P.DomainPrivate.IVName;
szXb = size(Xb);
nXb = szXb(1);
szYb = size(Yb);
nYb = szYb(1);

% Make sure basis function don't depend on variables that are not in P
Xbivn = cell(0,1);
for i=1:nXb
   Xbivn = union(Xbivn,Xb(i).BasisFunction.DomainPrivate.IVName);
end
nivX = length(Xbivn);

Ybivn = cell(0,1);
for i=1:nYb
   Ybivn = union(Ybivn,Yb(i).BasisFunction.DomainPrivate.IVName);
end
nivY = length(Ybivn);


[tfX,iX] = ismember(Xbivn,Pivn);  % Xbivn=P(iX)
if ~all(tfX)
   error('XBasis has dependence that P does not have.');
end
[tfY,iY] = ismember(Ybivn,Pivn);  % Ybivn=P(iY)
if ~all(tfY)
   error('YBasis has dependence that P does not have.');
end

% Arrange partials and rate bounds as listed in Xbivn and Ybivn
%   Xb(i).Partial is an something-by-2 with the partial in the first
%   column and the partial is taken wrt the variable in the second col.
Xp = cell(nXb,nivX);
XivLB = cell(1,nivX);
XivUB = cell(1,nivX);
for j = 1:nivX
   % Partials
   for k = 1:nXb
      tf = strcmp( Xbivn{j} , Xb(k).Partial(:,2));
      ntf = numel(find(tf));
      if ntf==1
         Xp{k,j} = Xb(k).Partial{tf,1};
      elseif ntf==0
         Xp{k,j} = pmat(0);
      else
         error(['Invalid partial: j = ' int2str(j) ', k = ' int2str(k)]);
      end
   end
   
   % Upper/Lower bounds on rate of Xiv
   tf = strcmp( Xbivn{j} , derLowerBound(:,2));
   ntf = numel(find(tf));
   if ntf==1
      XivLB{1,j} = derLowerBound{tf,1};
   else
      error(['Invalid LowerRate bound for Xbasis: j = ' int2str(j)]);
   end
   tf = strcmp( Xbivn{j} , derUpperBound(:,2));
   ntf = numel(find(tf));
   if ntf==1
      XivUB{1,j} = derUpperBound{tf,1};
   else
      error(['Invalid UpperRate bound for Xbasis: j = ' int2str(j)]);
   end
end

Yp = cell(nYb,nivY);
YivLB = cell(1,nivY);
YivUB = cell(1,nivY);
for j = 1:nivY
   % Partials
   for k = 1:nYb
      tf = strcmp( Ybivn{j} , Yb(k).Partial(:,2));
      ntf = numel(find(tf));
      if ntf==1
         Yp{k,j} = Yb(k).Partial{tf,1};
      elseif ntf==0
         Yp{k,j} = pmat(0);
      else
         error(['Invalid partial: j = ' int2str(j) ', k = ' int2str(k)]);
      end
   end
   
   % Upper/Lower bounds on rate of Yiv
   tf = strcmp( Ybivn{j} , derLowerBound(:,2));
   ntf = numel(find(tf));
   if ntf==1
      YivLB{1,j} = derLowerBound{tf,1};
   else
      error(['Invalid LowerRate bound for Ybasis: j = ' int2str(j)]);
   end
   tf = strcmp( Ybivn{j} , derUpperBound(:,2));
   ntf = numel(find(tf));
   if ntf==1
      YivUB{1,j} = derUpperBound{tf,1};
   else
      error(['Invalid UpperRate bound for Ybasis: j = ' int2str(j)]);
   end
end

% Create LMI variables
setlmis([]);
Xdec = zeros(nXb,1);
for k=1:nXb
   Xdec(k) = lmivar(1,[nX 1]);
end
Ydec = zeros(nYb,1);
for k=1:nYb
   Ydec(k) = lmivar(1,[nX 1]);
end
% Variable representing 1/Gamma
[ginv,ndec] = lmivar(1,[1 1]);

D = P.DomainPrivate;
ngpts = prod(D.LIVData);
cnt = 1;
for i=1:ngpts
   % Get parameter value
   pvaluec = num2cell(D(i));  % single index into RGRID gives value
   
   % Evaluate state matrices at parameter
   AhatV = double(lpvinterp(Ahat,D.IVName,pvaluec));
   AtilV = double(lpvinterp(Atil,D.IVName,pvaluec));
   BhatV = double(lpvinterp(Bhat,D.IVName,pvaluec));
   CtilV = double(lpvinterp(Ctil,D.IVName,pvaluec));
   b2V = double(lpvinterp(b2,D.IVName,pvaluec));
   c11V = double(lpvinterp(c11,D.IVName,pvaluec));
   d111dotV = double(lpvinterp(d111dot,D.IVName,pvaluec));
   d11dot1V = double(lpvinterp(d11dot1,D.IVName,pvaluec));
   c2V = double(lpvinterp(c2,D.IVName,pvaluec));
   b11V = double(lpvinterp(b11,D.IVName,pvaluec));
   
   % Evaluate basis functions
   XbV = zeros(nXb,1);
   for k=1:nXb
      XbV(k)=double(lpvinterp(Xb(k).BasisFunction,D.IVName,pvaluec));
   end
   YbV = zeros(nYb,1);
   for k=1:nYb
      YbV(k)=double(lpvinterp(Yb(k).BasisFunction,D.IVName,pvaluec));
   end
   
   % Evaluate partials and rate bounds of basis functions
   XpV = zeros(nXb,nivX);
   XivLBV = zeros(1,nivX);
   XivUBV = zeros(1,nivX);
   for j = 1:nivX
      XivLBV(1,j)=double(lpvinterp(XivLB{1,j},D.IVName,pvaluec));
      XivUBV(1,j)=double(lpvinterp(XivUB{1,j},D.IVName,pvaluec));
      for k = 1:nXb
         XpV(k,j)=double(lpvinterp(Xp{k,j},D.IVName,pvaluec));
      end
   end
   XivBV = [XivLBV; XivUBV];
   
   YpV = zeros(nXb,nivY);
   YivLBV = zeros(1,nivY);
   YivUBV = zeros(1,nivY);
   for j = 1:nivY
      YivLBV(1,j)=double(lpvinterp(YivLB{1,j},D.IVName,pvaluec));
      YivUBV(1,j)=double(lpvinterp(YivUB{1,j},D.IVName,pvaluec));
      for k = 1:nYb
         YpV(k,j)=double(lpvinterp(Yp{k,j},D.IVName,pvaluec));
      end
   end
   YivBV = [YivLBV; YivUBV];     % 2-by-nivY
   
   % X Ric LMI
   for q = 1:2^nivX  % number of variables appearing in X's basis fcns
      XPattern = ones(1,nivX);  % everything Lower
      tmp = dec2bin(q-1,nivX);
      XPattern(tmp=='1') = 2;  % grab Upper
      ind = sub2ind([2 nivX],XPattern,1:nivX);
      RBV = XivBV(ind);
      for k=1:nXb
         for j=1:nivX
            lmiterm([-cnt 1 1 Xdec(k)],RBV(j)*XpV(k,j),1);
         end
         lmiterm([cnt 1 1 Xdec(k)],XbV(k)*AhatV,1,'s');
         lmiterm([cnt 1 2 Xdec(k)],1,XbV(k)*c11V');
      end
      lmiterm([cnt 1 1 0],-b2V*b2V');
      lmiterm([cnt 1 3 ginv],1,BhatV);
      lmiterm([cnt 2 2 0],-eye(ne1));
      lmiterm([cnt 2 3 ginv],1,d111dotV);
      lmiterm([cnt 3 3 0],-eye(nd));
      cnt = cnt+1;
   end
   
   % Y Ric LMI
   for q = 1:2^nivY  % number of variables appearing in Y's basis fcns
      YPattern = ones(1,nivY);  % everything Lower
      tmp = dec2bin(q-1,nivY);
      YPattern(tmp=='1') = 2;  % grab Upper
      ind = sub2ind([2 nivY],YPattern,1:nivY);
      RBV = YivBV(ind);
      for k=1:nYb
         for j=1:nivY
            lmiterm([cnt 1 1 Ydec(k)],RBV(j)*YpV(k,j),1);
         end
         lmiterm([cnt 1 1 Ydec(k)],1,YbV(k)*AtilV,'s');
         lmiterm([cnt 1 2 Ydec(k)],1,YbV(k)*b11V);
      end
      lmiterm([cnt 1 1 0],-c2V'*c2V);
      lmiterm([cnt 1 3 ginv],1,CtilV');
      lmiterm([cnt 2 2 0],-eye(nd1));
      lmiterm([cnt 2 3 ginv],1,d11dot1V');
      lmiterm([cnt 3 3 0],-eye(ne));
      cnt = cnt+1;
   end
   
   % X/Y Coupling
   lmiterm([-cnt 1 2 ginv],1,eye(nX));
   for k=1:nXb
      lmiterm([-cnt 1 1 Xdec(k)],XbV(k),1);
   end
   for k=1:nYb
      lmiterm([-cnt 2 2 Ydec(k)],YbV(k),1);
   end
   cnt = cnt + 1;
   
end

% Create objective function and get lmis
cobj = zeros(ndec,1);
cobj(end) = -1;
lmisys = getlmis;

% Solve LMI
opts = zeros(5,1);
[copt,xopt] = mincx(lmisys,cobj,opts);


%
X = zeros(nX,nX);
for k=1:nXb
   X = X + Xb(k).BasisFunction*dec2mat(lmisys,xopt,Xdec(k));
end
Y = zeros(nX,nX);
for k=1:nYb
   Y = Y + Yb(k).BasisFunction*dec2mat(lmisys,xopt,Ydec(k));
end

xinvpart = zeros(nX,nX);
for j=1:nivX
   newVar = idf([Xbivn{j} 'Dot'],[-100 100]);
   tmpMat = zeros(nX,nX);
   for k=1:nXb
      tmpMat = tmpMat + Xp{k,j}*dec2mat(lmisys,xopt,Xdec(k));
   end
   xinvpart = xinvpart + newVar*tmpMat;
end
xinvpart = -inv(X)*xinvpart*inv(X);

% Controller Reconstruction
d22 = zeros(nmeas,ncont);
d12 = [zeros(ne1,ne2);eye(ne2,ne2)];
d21 = [zeros(nd2,nd1) eye(nd2,nd2)];

gam = -1/copt;
gams = gam^2;
gamsi = 1/gams;

omeg = -d1122 - d1121*((gams*eye(nd1,nd1) - d1111'*d1111)\d1111')*d1112;

q = Y - gamsi*inv(X);

b2omeg = b2*omeg;
d12omeg = d12*omeg;
abar = a + b2omeg*c2;
b1bar = b1 + b2omeg*d21;
c1bar = c1 + d12omeg*c2;
d11bar = d11 + d12omeg*d21;

dh = inv(eye(ne,ne)- gamsi*d11bar*d11bar');
dt = inv(eye(nd,nd) - gamsi*d11bar'*d11bar);

F = (-d12'*dh*d12)\((b2+b1bar*d11bar'*dh*d12*gamsi)'/X + d12'*dh*c1bar);
L = -(Y\(c2+d21*dt*d11bar'*c1bar*gamsi)' + b1bar*dt*d21')/(d21*dt*d21');

Af = abar + b2*F;
cf = c1bar + d12*F;

Afx = X\Af;
left = X\b1bar + cf'*d11bar;

H = -(Afx+Afx'+xinvpart+cf'*cf+gamsi*left*dt*left');

qiy = q\Y;

m1 = H + F'*(b2'/X + d12'*cf);
m2 = gamsi*(gams*q*(-qiy*L*d21 - b1bar) + F'*d12'*d11bar)*dt*left';
m2 = (q*(-qiy*L*d21 - b1bar) + gamsi*F'*d12'*d11bar)*dt*left';
m = m1 + m2;

AC = Af + qiy*L*c2 - gamsi*(q\m);
BC = -qiy*L;
CC = F;
DC = omeg;

ortcont = ss(AC,BC,CC,DC);
K = lft([TR;eye(ncont)]*ortcont*[TL' eye(nmeas)],FT);

% Solve Relaxed LMI
%
Gamma = gam;
[Ke,Xe,Ye,copte,xopte] = ...
   lpvof2miPartB(sys,nmeas,ncont,gam*BOF,...
   Xb,Yb,derUpperBound,derLowerBound);

Info = struct('Xe',Xe,'Ye',Ye,'xopte',xopte,'copte',copte,...
   'K',K,'X',X,'Y',Y,'xopt',xopt,'copt',copt);

