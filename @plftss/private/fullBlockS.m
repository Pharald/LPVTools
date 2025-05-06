function [Gpart,Pi,ndec,cnt,info] = fullBlockS(G,cnt)

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
[G0,delta,blkstruc,~] = lftdata(G.Data);


% get dimensions
nin = size(G,2);
nout = size(G,1);
ndelta = size(delta,1);

% partition G0
G011 = G0(1:ndelta,1:ndelta);
G021 = G0(ndelta+1:end,1:ndelta);
G022 = G0(ndelta+1:end,ndelta+1:end);
G012 = G0(1:ndelta,ndelta+1:end);

% reconstruct
Gpart = [G011, G012; eye(ndelta), zeros(ndelta,nin); G021, G022];


info.ndelta = ndelta;

% Define multiplier
% ---- Can add options and more types of multipliers

% convert the blk structure from lftdata to a mat structure for setting up
% LMI variable
blk = zeros(length(blkstruc),2);

for ii = 1:length(blkstruc)
    if ~strcmp('ureal',blkstruc(ii).Type)
        error('Parameter need to be tvreal')
    end
    blk(ii,1) = blkstruc(ii).Occurrences;
    blk(ii,2) = 1;
end

[R,ndec,sR] = lmivar(1,blk); % block diagonal multiplier to commute with Delta

sS = [];
for ii = 1:size(blk,1)
    sS = blkdiag(sS,skewdec(blk(ii,1),ndec));
    ndec = ndec + blk(ii,1)*(blk(ii,1)-1)/2;
end

[S,ndec,sS] = lmivar(3,sS);         % skew symmetric block diagonal multiplier


Pi = lmivar(3,[-sR, sS; sS', sR]);

info.Pidim = [ndelta*2, ndelta*2]; % dimensions of multiplier

% lmi condition
lmiterm([cnt 1 1 R],1,1); 
% only one condition so cnt doesn't change
% -------------------------------------------------



end





