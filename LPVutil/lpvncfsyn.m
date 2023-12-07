function [K,CL,GAM,INFO]=lpvncfsyn(G,W1,W2,varargin)
% LPVNCFSYN   Parameter-dependent Glover-McFarlane loop-shaping for PSS.
%
% [K,CL,GAM,INFO]=LPVNCFSYN(G,W1,W2) synthesizes a parameter-dependent 
% Glover-McFarlane loop-shaping controller K for the shaped plant Gs=W2*G*W1.
%
%                    ^ z1       ^ z2
%         ____       |   ____   |        ____    ____
%      ->| W2 |-> + ---->| Ks |--> + -->| W1 |--| G  |--- 
%   + |   ----    ^       ----     ^     ----    ----    |
%     |           |                |                     |
%     |        w1 |             w2 |                     |
%     |__________________________________________________|
%
% G is a PSS, while W1 and W2 can be PSS, PMAT or DOUBLE. K is a 
% parameter-varying controller which minimizes the induced L2 norm 
% of the loop-shaping interconnection defined by G, W1 and W2. K is a PSS 
% defined on same domain as P. K has as many inputs as G has outputs, and 
% as many outputs as G has inputs. GAM is the induced L2 norm of the 
% loop-shaping interconnection. CL is the system taking in [w1; w2] and
% outputting [z1; z2]. INFO is a structure containing data from the Linear 
% Matrix Inequalities that are solved to obtain K. A call to LPVNCFSYN 
% without a BASIS function argument generates a controller assuming no 
% bounds on the parameter rate of variation. 
%
% [K,CL,GAM,INFO]=LPVNCFSYN(G,W1,W2,'ref') synthesizes a parameter-dependent 
% Glover-McFarlane loop-shaping controller K for the shaped plant Gs=W2*G*W1,
% assuming  a reference command. CL is the system taking in [w1; w2; ref] 
% and outputting [z1; z2], and GAM is its induced L2 norm.
%
%                          ^ z1         ^ z2
%             ____         |   ____     |     ____    ____
%      ----->| W2 |-> + ----->|    |--> + -->| W1 |--| G  |--- 
%   + |       ----    ^       | Ks }      ^   ----    ----    |
%     |               |   --->|    |      |                   |
%     |            w1 |   ref  ----    w2 |                   |
%     |                                                       |
%     |_______________________________________________________|
%
%
% [K,CL,GAM,INFO]=LPVNCFSYN(G,W1,W2,Xb,Yb,...) computes the rate-bounded 
% Glover-McFarlane loop-shaping controller K where the rate-bounds of the 
% independent variables of the shaped plant Gs are incuded in the synthesis 
% conditions. Xb and Yb are BASIS objects, which describe the assumed 
% parameter dependence of the lyapunov matrices used in solving for K.
%
% [K,CL,GAM,INFO]=LPVNCFSYN(G,...,OPT) allows the user to pass in
% a LPVSYNOPTIONS object. 
%
% See also: ncfsyn, lpvsynOptions, lpvsyn, lpvsfsyn, lpvestsyn, lpvmixsyn, lpvloopshape.

% Note that help in NCFSYN has the incorrect diagrams MATLAB 2013a.

% Parse inputs
nin = nargin;
if nin==1
    W1=[];
    W2=[];
    rflag = false;
    varargin = cell(0,0);
elseif nin==2
    W2=[];
    rflag = false;
    varargin = cell(0,0);
elseif nin==3
    rflag = false;
    varargin = cell(0,0);
elseif nin>=4
    v1 = varargin{1};
    if ischar(v1)
        % Check if user specified 'ref' option. Remaining inputs must be
        % basis functions and lpvsynOptions. 
        len = min(3,length(v1));
        rflag =  strncmpi(v1,'ref',len);
        varargin = varargin(2:end);
    else 
        % Handle call with rate-bounds
        rflag = false;
    end          
end

% Default weights
[ny,nu] = size(G);
if isempty(W1)
    W1 = ss(eye(nu));
end
if isempty(W2)
    W2 = ss(eye(ny));
end

% Build generalized interconnection structure
% XXX Note that the Hinf Loopshaping problem has special structure
% that is not being exploited in this code.
Gs = W2*G*W1;
if rflag
    % [z1;z2;y] = Pic*[win1; win2; ref; u]
    nref = size(W2,2)-ny;
    
    systemnames = 'Gs';
    inputvar = ['[win1{' int2str(ny) '}; win2{' int2str(nu) ...
        '}; ref{' int2str(nref) '}; u{' int2str(nu) '}]']
    outputvar = '[Gs+win1; u; Gs+win1; ref]';
    input_to_Gs = '[u+win2]';
    Pic = sysic;
else 
    % [z1;z2;y] = Pic*[win1; win2;u]
    nref = 0;
        
    systemnames = 'Gs';
    inputvar = ['[win1{' int2str(ny) '}; win2{' int2str(nu) '}; u{' int2str(nu) '}]'];
    outputvar = '[Gs+win1; u; Gs+win1]';
    input_to_Gs = '[u+win2]';
    Pic = sysic;    
end


% Synthesize controller Ks for shaped plant

% opt = lpvsynOptions('BackOffFactor',1.01);
[Ks,GAM,INFO] = lpvsyn(Pic,ny,nu,varargin{:}); 
K=W1*Ks*W2;

% Store outputs if requested by user
nout = nargout;
if nout==4
    INFO.Ks = Ks;
    INFO.emax = 1/GAM;
    INFO.Gs = Gs;
end
if nout>=2
    CL = lft(Pic,Ks);
end
