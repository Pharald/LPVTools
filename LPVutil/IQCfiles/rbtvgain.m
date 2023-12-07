function [IQCcell] = rbtvgain(RB,cpoleVec)
% delta in in [-1 1]; |deltadot| <= RB
% cpoleVec: vector of real, negative poles associated with the separate
% first-order terms

nIQC = numel(cpoleVec);
IQCcell1 = cell(1,nIQC);
IQCcell2 = cell(1,nIQC);
IQCcell0 = [1 0;0 -1];
for i=1:nIQC
    beta = -cpoleVec(i);
    % Compute phi for h(t) = e^(-beta*t)
    phi = RB/beta^2*(1 - exp(-2*beta/RB));
    H = tf(1,[1 beta]);
    rho = 0.25;
    % Eq 28 in MR
    IQCcell1{i} = [(1+rho)*(H'*H + phi^2/rho) 0; 0 -H'*H];
    % Eq 29 in MR
    H = H - H';
    IQCcell2{i} = [2*phi H;H' 0];
end

IQCcell = [IQCcell0 IQCcell1 IQCcell2];
    
    
    



