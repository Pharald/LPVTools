function [K,CL,GAM,INFO] = lpvmixsyn(G,varargin)

% LPVMIXSYN Mixed-sensitivity synthesis for PSS
%
% lft(G,K): [z1] = [W1 0 ][S  -SG ][W2 0 ][w1]
%           [z2]   [0  W3][KS -KSG][0  W4][w2]
%
% CL: with performance ouptuts, inputs={r,d}, outputs={e,u,y}
% 
% [K,CL,GAM,INFO]=mixsyn(G,W1,W2,W3,W4,...) performs an induced L-2 mixed
% sensitivity synthesis for the gerneralised plant constructed of G and the
% weigths W.
%
% [K,CL,GAM,INFO]=mixsyn(G,W1,W2,W3,W4,Xb,Yb)
% Xb and Yb are BASIS objects, which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for K.
%
% See also: lpvsyn, mixsyn, loopsyn.


nin = nargin;
nout = nargout;
narginchk(2, inf)
nargoutchk(0, 4)

if ~isa(varargin(end),'basis')
Xb = [];
Yb = Xb;
else
Xb = varargin(5);
Yb = varargin(6);
end

% assign weights
We = varargin(1);
Wr = varargin(2);
Wu = varargin(3);
Wd = varagin(4);
nmeas = size(G,1); % e
ncon = size(G,2); % u

% all output signals are in the feeback channel
if ~((size(Wr,1) == size(We,2)) == nmeas)
error('the output of W2 and input to W1 must have the same dimensions as each other and as the output from G')
end

if ~((size(Wd,1) == size(Wu,2))== ncon)
error('the output of W4 and input to W3 must have the same dimensions as each other and as the input to G')
end

% checks
if ~isa(Wd,'double')
    error('W4 must be a scalar (type double)')
end

% construct generalised plant
systemnames = 'G Wr We Wu Wd';
inputvar = strcat('[w1{',num2str(size(Wd,2)),'}; w2{',num2str(size(Wd,2)),'}; u{',num2str(ncon),'}]');
outputvar = strcat('[We; Wu; Wr-G]');
input_to_G = '[u + Wd]';
input_to_Wr = '[w1]';
input_to_We = '[Wr - G]';
input_to_Wu = '[u]';
input_to_Wd = '[w2]';
cleanupsysic = 'yes';

Pgen = sysic;

% perform synthesis
[K,GAM,INFO] = lpvsyn(Pgen,nmeas,ncon,Xb,Yb);

% close loop
systemnames = 'G K';
inputvar = strcat('[r{',num2str(nmeas),'}; d{',num2str(ncon),'}]');
outputvar = '[r-G; K; G]';
input_to_G = '[K + d]';
input_to_K = '[r - G]';
cleanupsysic = 'yes';

CL = sysic;

end