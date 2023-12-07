% Orthogonalizes a general OLIC interconnection for Ratebounded LPV design
% function [osys,Tleft,Tright,ftterm,unil,unir] = orthog4syn(P,nmeas,nctrl)
function [osys,Tleft,Tright,ftterm,ne,nd,unil,unir] = orthog4syn(P,nmeas,nctrl)

if nargin ~= 3
   disp('usage: [osys,invl,invr,ftterm,unil,unir] = orthog4syn(sys,nmeas,nctrl)');
   return
end

szP = size(P);
nd = szP(2) - nctrl;
ne = szP(1) - nmeas;

[a,b,c,d,Ts] = ssdata(P);
if Ts~=0
   error('Discrete-Time is not addressed yet.  Likely the same code is OK.');
end
b1 = b(:,1:nd);
b2 = b(:,nd+1:end);
c1 = c(1:ne,:);
c2 = c(ne+1:end,:);
d11 = d(1:ne,1:nd);
d12 = d(1:ne,nd+1:end);
d21 = d(ne+1:end,1:nd);
d22 = d(ne+1:end,nd+1:end);
% Determine if D12 has full column rank and scale D12 to
%    Q12*D12*R12INV = [0;I].
% Hence if the system is redefined with a unitary transformation
% on ERROR, etilde := Q12 e, and intertible transformation on CONTROL,
% u := R12INV utilde, in the new variables, D12Tilde = [0;I].  Note that
% the unitary transformation on e does not change ||e||, and the invertible
% transformation on u can be included in the overall controller.
[q12,r12] = qr(d12);
rrk = double(lpvmin(rank(r12)));
if (rrk ~= nctrl)
   error(' D12 DOES NOT HAVE FULL COLUMN RANK over IV')
end
q12 = q12(:,[nctrl+1:ne 1:nctrl]);
Tright = inv(r12(1:nctrl,:));  % R12INV
unil = q12';
% Determine if D21 has full col rank and scale D21 to
%    R21INV'*D21*Q21 = [0 I].
[q21,r21] = qr(d21');
crk = double(lpvmin(rank(r21)));
if (crk ~= nmeas)
   error(' D21 DOES NOT HAVE FULL ROW RANK')
end
unir = q21(:,[nmeas+1:nd 1:nmeas]);
Tleft = inv(r21(1:nmeas,:));   % R21INV

ftterm = -Tleft'*d22*Tright;

c1 = unil*c1;
c2 = Tleft'*c2;
b1 = b1*unir;
b2 = b2*Tright;
d11 = unil*d11*unir;
d12 = [zeros(ne-nctrl,nctrl);eye(nctrl)];
d21 = [zeros(nmeas,nd-nmeas) eye(nmeas)];
d22 = zeros(nmeas,nctrl);
c = [c1;c2];
b = [b1 b2];
d = [d11 d12;d21 d22];
osys = ss(a,b,c,d);

