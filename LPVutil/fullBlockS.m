function [Gpart,Pi,ndec,cnt,info] = fullBlockS(Gp,ndec,cnt)

% for X(p) < 0 that satisfies X(p) = G'(p) X0 G(p) and G(p) = Fu(G0,Delta(p))
% and X0 is symmetrical
%
% an equivalent set of LMIs is
%
% Gpart'*blkdiag(Pi,X0)*Gpart < 0
% [I;Delta(p)]' * Pi *[I; Delta(p)] >= 0
% 
% Where Pi is a multiplier

% cnt is the current LMI number

nargoutchk(4,5); % not necessary to output info

% partition Gp and determine size of Delta(p)
[Gpart,ndelta] = partition_mat(Gp);

info.ndelta = ndelta;

% Define multiplier
% ---- Can add options and more types of multipliers

% symmetric multiplier such that parameter dependent LMI is fullfilled 
sS = skewdec(ndelta,ndec);
[~,ndec,sS] = lmivar(3,sS);         % skew symmetric
sR = diag(ndec+1:ndelta+ndec);
[R,ndec,~] = lmivar(3,sR);                 % diagonal
Pi = lmivar(3,[-sR, sS; sS', sR]);

info.Pidim = [ndelta*2, ndelta*2]; % dimensions of multiplier

% lmi condition
lmiterm([cnt 1 1 R],1,1); 
% only one condition so cnt doesn't change
% -------------------------------------------------



end

function [Gpart,ndelta] = partition_mat(G)

% separate G into lft(G0,deltar); 

if isa(G,'uss') || isa(G,'umat')
G = simplify(G,'full');
    [G0,delta] = lftdata(G);

elseif isa(G,'plftss') || isa(G,'plftmat')
    G = simplify(G.Data,'full');
[G0,delta,~,~] = lftdata(G);
end

% get dimensions
n_in = size(G,2);
n_out = size(G,1);
ndelta = size(delta,1);

% partition G0
G011 = G0(1:ndelta,1:ndelta);
G021 = G0(ndelta+1:end,1:ndelta);
G022 = G0(ndelta+1:end,ndelta+1:end);
G012 = G0(1:ndelta,ndelta+1:end);

% reconstruct
Gpart = [G011, G012; eye(ndelta), zeros(ndelta,n_in); G021, G022];
end
