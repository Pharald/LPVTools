function [D,DN] = LOCALlhcdesign(n,p,Range)
% Latin Hypercube Design
% [0 1]^n, p points total
%
% % Example
% % Inputs:
% n = 4;  % 4-dimensional space
% p = 50; % 50 samples
% Range = [-1 1;0 5;-2 0;4 5];
% D = LOCALlhcdesign(n,p,Range);
% D  % 4-by-50 sample matrix in Range

V = rand(p,n);
ILength = 1/p;
IBegin = (0:p-1)'*ILength;

V = repmat(IBegin,[1 n]) + ILength*V;
for i=1:n
   [~,idx] = sort(rand(1,p));
   V(:,i) = V(idx,i);
end
DN = V';
% Range n-by-2
Left = Range(:,1);
Width = Range(:,2) - Range(:,1);
D = repmat(Left,[1 p])+DN.*repmat(Width,[1 p]);
