function [Pout,Ratio,MinAbove,MaxBelow,T] = tmpresid2(P,Nremove)
%% TODO  Help, user input

A = P.DataPrivate.a;

[Ratio,MinAbove,MaxBelow,T] = splitlpv( A );
nr = length(Ratio);
i2 = 1:nr-1;
semilogy(i2,MinAbove(i2),'+',i2,MaxBelow(i2),'*',i2,Ratio(i2),'r')
legend('MinAbove','MaxBelow','Ratio')
sztmp = size(P.DataPrivate);

tmp = ss2ss(P.DataPrivate,inv(T(:,:,Nremove)));
ztmp = tmp;
for j=1:prod(sztmp(3:end))
    ztmp(:,:,j) = modred(tmp(:,:,j),1:Nremove);
end
Pout = pss(ztmp,P.DomainPrivate);
