function [Pred,info] = lpvbalancmr(P,varargin)
% LPVBALANCMR  Balanced truncation of quadratically stable PSS models.
%
% [PRED,INFO] = LPVBALANCMR(SYS,N) performs balanced truncation model 
% reduction on a PSS SYS. A PSS SYS, with Nx states, is balanced by 
% computing a single balancing transformation for SYS and applying it at 
% every point in its domain. The output PRED has N states and is obtained 
% by truncating from the balanced system the (Nx-N) states which contribute 
% least to its input-output mapping. INFO contains two fields 'StabSV' and 
% 'ErrorBound'.  INFO.StabSV is a vector of singular values describing the 
% input-output mapping of SYSB (comparable to Hankel singular values for 
% LTI systems). INFO.ErrorBound contains the L_2 norm of the difference 
% between SYS and PRED: INFO.ErrorBound = ||SYS-PRED||_2. 
%
% Note that LPVBALANCMR only works for quadratically stable systems. For
% unstable PSS models use LPVNCFMR.
%
% [PRED,INFO] = LPVBALANCMR(SYS,N,OPTION1,VAL1,OPTION2,VAL2,...) provides
% additional options for the balanced truncation model reduction. 
% The current implementation supports the following options:
%
%   OPTION    |       VAL          |           Explanation
%-------------------------------------------------------------------------
%  'weight'   | {Wout,Win}         | LTI weights on input (Win) and output
%             |                    | (Wout). Used to emphasize accuracy in 
%             |                    | different I/O and frequency ranges. 
%             |                    | Must be invertable if method 'invgram'
%             |                    | is used. 
%-------------------------------------------------------------------------
%  'method'   |'gram' or 'invgram' | Solve using either gramians or the
%             |                    | inverse gramians.
%-------------------------------------------------------------------------
%
% See also:	lpvncfmr, lpvbalreal, lpvgram, modred.

% TODO - AH 9/11/14
%      - Implement a call where the funcion loops over different
%        orders of reduction. This is needed to implement 'MaxError'
%        call. It is also needed to replicate functionality of standard
%        MATALB ncfmr which if called as e.g. balancmr(G,[10:2:18]);
%        will return an array of reduced order models 
%        with 10,12,14,16, and 18 states respectively. 
% TODO - AH 9/11/14
%      - Implement 'MaxError' and 'Display' input arguments.

% Parse system data:
[A,B,C,D] = ssdata(P);
nX = size(A,1);
nU = size(B,2);
nY = size(C,1);

solflag =true;
flagout = true;

if nargin == 1
    W = pss(eye(nY));    
    order = nX;
elseif nargin ==2
    W = pss(eye(nY));   
    order = varargin{1};
else
    if isa(varargin{1},'double')
        order = varargin{1};
        varargin(1) =[];
    end
    
    % XXX do more input parsing here- for now assume inputs are weights
    nkey = numel(varargin)/2;
    nvpair = reshape(varargin,[2 nkey]);
    for i=1:nkey
        name = nvpair{1,i};
        val = nvpair{2,i};
        
        if strcmpi(name,'weights')
            Weights = val;
            Wo = Weights{1};
            Wi = Weights{2};
            if isempty(Wi)
                W = Wo;
                if isscalar(W)
                    W = W*eye(nY);
                end
                flagout = true;
            elseif isempty(Wo)
                W = Wi;
                if isscalar(W)
                    W = W*eye(nU);
                end
                flagout = false;
            else
                error('Must specify only one weight at a time')
            end
            
        elseif strcmpi(name,'method')
            if strcmpi(val,'invgram')
                solflag = false;
            end
        else
            error('XXX handle later')
        end
    end
end

% BALANCE MODEL:

[Aw,Bw,Cw,Dw] = ssdata(W);
nW = size(Aw,1);

% XXX check definition of these matrices - eq. 7.22 onward in Wood thesis.
if flagout 
    % output weight
    Abar = [A zeros(nX,nW);Bw*C Aw];
    Bbar = [B ; zeros(nW,nU)];
    Cbar = [C Cw];
else
    % input weight
    Abar = [A B*Cw;zeros(nW,nX) Aw];
    Bbar = [B;Bw];
    Cbar = [C zeros(nY,nW)];
end
Dbar = zeros(nY,nU);

% Define weighted LPV Model
Pw = ss(Abar,Bbar,Cbar,Dbar);
% Balance weighted LPV Model
[Pwbal,G,T,Ti,Wow,Wcw] = lpvbalreal(Pw,solflag);

% Define Weighted gramians
Wo = Wow(1:nX,1:nX);
Wc = Wcw(1:nX,1:nX);


% Balance System ( Moore's algorithm, using notation in section 9.5 in
% Robust Systems by Sanchez-Pena and Sznaier)
[Vc,Sc] = LOCALschursort(Wc);
Scrt = sqrt(diag(Sc));
T1 = lrscale(Vc',1./Scrt,[]); %T1 = inv(Scrt)*Vc'
T1i = lrscale(Vc,[],Scrt); %T1i = Vc*Scrt
Wotil = (T1')\(Wo/T1);
[Vo,So] = LOCALschursort(Wotil);
Sort = diag(So).^(1/4);
T2 = lrscale(Vo',Sort,[]); %T2 = Scrt*Vo'
T2i = lrscale(Vo,[],1./Sort); %T2i = Vo*inv(Scrt)

% Form final transformation
T = T2*T1;
Ti = T1i*T2i;

% Compute balanced original system:
Pbal = ss2ss(P,T);

% Compute Hankel singular values:
G = sqrt(diag(So));

% XXX do a loop on order, if order is a vector creating multiple reduced 
% models for different levels of order reduction 
% info.ErrorBound becomes a vector
% Pred becomes an array of pss

% MODEL REDUCTION
elim = (order+1):nX;
Pred = modred(Pbal,elim,'Truncate');
info.ErrorBound = 2*sum(G(elim));
info.StabSV = G;

function [V,S] = LOCALschursort(W)
    [V,S] = schur(W);
    [Ssort,idx] = sort(diag(S),'descend');
    S = diag(Ssort);
    V = V(:,idx);      





            
    
