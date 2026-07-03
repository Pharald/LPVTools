function [F,Gamma,Info] = lpvsfsyn(P,ncont,Xb,alg,opt)
% LPVSFSYN  Parameter-dependent state feedback controller synthesis for PSS
%
% [F,GAMMA,INFO] = LPVSYN(P,NCONT) computes a parameter-varying
% state-feedback gain F which minimizes the induced L2 norm of the 
% interconnection defined by lft(P,F). F is a PMAT with NCON outputs, 
% defined on same domain as P. GAMMA is the induced L2 norm of lft(P,K).
% This three argument call assumes that the rate-bounds of the independent
% variables in P are [-inf,inf]. INFO is a structure containing data from
% the Linear Matrix Inequalities that are solved to obtain F.
%
% [F,GAMMA,INFO] = LPVSYN(P,NCONT,Xb) computes the rate-bounded 
% parameter-varying state-feedback gain F for a system P. F is the 
% controller which minimizes the induced L2 norm of lft(P,F) when the 
% rate-bounds of the  independent variables of P are incorporated into
% the synthesis. Xb is a BASIS object, which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for F.
%
% [F,GAM,INFO] = LPVSFSYN(P,NCONT,Xb,ALG) computes a parameter-varying 
% state-feedback gain F, which minimizes the induced L2 norm if ALG = 'L2' 
% or the stochastic LPV bound if ALG = 'LQG'. The stochastic LPV bound is 
% defined as the expected value of the average instantaneous power of the 
% output of P, assuming its inputs are zero mean, white-noise processes 
% with unit intensity.
%
% [F,GAMMA,INFO] = LPVSYN(P,NCONT,Xb,ALG,OPT) allows the user to pass in
% a LPVSYNOPTIONS object. 
%
% The default algorithm for LPVSFSYN will solve the given synthesis problem
% twice. The first iteration attempts to find a solution that minimizes the
% induced L2 norm of lft(P,K). The second iteration will solve the 
% optimization problem again, with the caveat that any solution that is 
% found to lie within 15% of the optimal induced L2 norm of lft(P,K) from 
% the first iteration, is satisfactory. This formulation has been found to 
% yield controllers that are better numerically conditioned. The back-off 
% factor of 15% can be reset to a different value in LPVSYNOPTIONS.
%
%
% See also: lpvsynOptions, lpvestsyn, lpvsyn, lpvncfyn, lpvmixsyn, lpvloopshape.



% XXX - OPTIONS are currently only used if 'L2' is selected. 
% [L,GAM,INFO] = LPVESTSYN(P,NMEAS,...,opt) permits the user to pass in 
% See also: lpvsynOptions.

% TODO PJS 5/17/2011:  Combine rate bounds into a single niv-by-3 cell array?
% RateUB, (<=niv)-by-2, cell array, {PMAT/DOUBLE variableName}
% RateLB, (<=niv)-by-2, cell array, {PMAT/DOUBLE variableName}

% Parse Inputs
% XXX Check later
% XXX OPT is not being used.
nin = nargin;
narginchk(2, 5)
nout = nargout;
if nin==2
    opt = lpvsynOptions;
    Xb = [];
    alg = 'L2';
elseif nin==3
    if isa(Xb,'lpvsynOptions')
        opt = Xb;
        Xb = [];
        alg = 'L2';
    elseif isa(Xb,'basis')
        alg = 'L2';
        opt = lpvsynOptions;
    elseif isa(Xb,'char')
        alg = Xb;
        Xb = [];
        opt = lpvsynOptions;
    else
        error(['The third argment in a call to lpvsfsyn must be '...
            'either a BASIS object, a CHAR, or a lpvsynOptions'])
    end
elseif nin==4
    if isa(alg,'lpvsynOptions')
        opt = alg;
        alg = 'L2';
    elseif isa(alg,'char')
        opt = lpvsynOptions;
    else
        error(['The fourth argment in a call to lpvsfsyn must be '...
            'either a a CHAR, or a lpvsynOptions'])
    end
elseif nin==5
    if ~isa(opt,'lpvsynOptions')
        error(['The fifth argment in a call to lpvsfsyn must be '...
            'a lpvsynOptions'])
    end
end

% Always assume one basis function: constant
if isempty(Xb)
   Xb = basis(1,0); 
end

% XXX Does not work for arrays of LPV systems.
% Check for this and error out if there are array (non-LPV) dimensions
% if hasArray(sys.Domain), error, end

% Map the input data into engine data form:
[Pdata,RateBounds,Fbasis,Fgrad] = basis2data(P,Xb);

% Single balancing transformation
% TODO PJS 8/30/2013: Do this before or after orthogonalization?
% (orthogonalization occurs in the engine)
% Need to undo the coordinate transformation in the state-feedback gain.
% nd = size(P,2) - ncont;
% ne = size(P,1);
% blk = [nd ne; ncont nx];
% P = lpvbalance(P);


if strcmpi(alg,'L2')
    % Run the L2 synthesis engine:
    [F,Gamma,Info] = lpvL2sfsynengine(Pdata,ncont,Fbasis,Fgrad,RateBounds,opt);

    % TODO HP 07/12/2023: Add lpvsynoptions for stochastic synthesis engine
    % The options are only used for L2 synthesis at the moment
elseif strcmpi(alg,'LQG')
    % Run the Stochastic synthesis engine:
    [F,Gamma,Info] = lpvLQGsfsynengine(Pdata,ncont,Fbasis,Fgrad,RateBounds,opt);
end



nx = size(Pdata.A,1);
Dom = P.Domain;
F = pmat(reshape(F,[ncont; nx; Dom.LIVData]'),Dom);




