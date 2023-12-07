function [Ratio,MinAbove,MaxBelow,T] = splitlpv( A )

% A is n x n x npts
n = size(A,1);
ns = size(A);
npts = prod(ns(3:end));
V = zeros(n,n*npts);
VL = zeros(n,npts);
T = zeros(n,n,npts);

for i=1:npts
     [evc,evl] = eig(A(:,:,i));
     [VL(:,i),idx] = sort(diag(real(evl)));
     V(:,(i-1)*n+1:(i*n)) = evc(:,idx);
end

MaxBelow = zeros(1,n-1);
MinAbove = zeros(1,n-1);
Ratio = zeros(1,n-1);
B = 0:n:(n-1)*npts;

for i=1:n-1
    fast = 1:i;
    evcidx = repmat(fast,[1 npts]) + kron(B,ones(1,i));
    if numel(evcidx)<n
        [T(:,:,i),S] = svd(V(:,evcidx));
    else
        [T(:,:,i),S] = svd(V(:,evcidx),'econ');
    end
    Ratio(i) = S(i+1,i+1)/S(i,i);
    MinAbove(i) = -max(VL(i,:));
    MaxBelow(i) = -min(VL(i+1,:));
end
