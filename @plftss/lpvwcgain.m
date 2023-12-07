function [WCG,WCU,Info] = lpvwcgain(P,omega,polelist,flag)
% LPVWCGAIN  Worst-case bound on induced L2 norm for PLFTSS systems.
% 
% [WCG,WCU,INFO] = LPVWCGAIN(P) computes the upper-bound on the worst-case
% induced L2 norm of the PLFTSS P. "Worst-case" refers to all the modeled
% uncertainty and parameter-dependence. WCG is the  upper bound on the 
% worst-case induced L2 norm. INFO is a structure containing additional 
% information about the solution, including an estimate of the lower bound
% on the worst-case induced L2 norm, based on LTI worst-case gain analysis.
% WCU is value of the uncertainty associated with the lower-bound of the
% induced L2 norm, which is based on LTI worst-case gain analysis.
%
% [WCG,WCU,INFO] = LPVWCGAIN(P,OMEGA) allows the user to specify a custom 
% frequency vector OMEGA for the analysis.
%
% [WCG,WCU,INFO] = LPVWCGAIN(P,...,POLELIST) allows the user to define
% basis functions for the Integral Quadratic Constraints (IQC) used to 
% bound the uncertainty when formulating the LMI to be solved. POLELIST is 
% a 1xN DOUBLE row vector, of negative values. Each value in POLELIST 
% corresponds to a pole of a stable transfer function that is used as a 
% weight on all signals in the IQCs. A default POLELIST, with five pole 
% values, is used when a POLELIST is not supplied by the user. The five 
% pole values are selected automatically from the frequency range of the 
% system dynamics.
%
% See also: lpvnorm, norm, wcgain.

nin = nargin;
narginchk(1,4);

if nin ==1
    omega = [];
    polelist = [];
    flag = true;
elseif nin==2
    polelist = [];
    flag = true;
elseif nin ==3
    flag = true;
end

[m,delta,blkstruct,normunc] = lftdata(P,[],'All');

S.type = '.';
S.subs = 'Blocks';
unames = fieldnames(subsref(P,S));
[omega,defaultpl] = polefreqauto(m,unames,polelist,omega);

for k = 1:numel(normunc)
    nunc = normunc{k};
    if isa(nunc,'tvreal')
        
        rb = max(abs(nunc.Ratebounds));
        
        if isinf(rb)
            % No rate bound IQC analysis
            S = IQCtvrealnrb;
        else
            % Rate bounded IQC analysis
            vnam = nunc.Name(1:end-numel('Normalized'));
            polelist_k = defaultpl.(vnam);
            idx  = isinf(polelist_k);
            polelist_k(idx) = [];
            S = IQCtvreal(polelist_k,rb);
        end
        % Replace the tvreal with a udyn with IQCs
        normunc{k} = udyn(nunc.Name,[1 1],'UserData',S);
    end
end
% Recreate uncertain system with IQCs replacing tvreals.
Piqc = lft(blkdiag(normunc{:}),m);

% Compute LPV norm (L2 gain)
[UpperBound,Sopt,xopt] = iqcperf4(Piqc,omega);

WCU = [];
% Package the output:
if flag == true
    % Estimate the lowerbound
    pfr = ufrd(P.Data,omega);
    [ltiWCG,WCU,wcinfo] = wcgain(pfr);
    WCG.LowerBound = ltiWCG.LowerBound;    
    Info.Sensitivity = wcinfo.Sensitivity;
    Info.Frequency = omega;
    Info.BadUncertainValues = wcinfo.BadUncertainValues;
    Info.ArrayIndex = wcinfo.ArrayIndex;
    Info.CriticalFrequency = ltiWCG.CriticalFrequency;
end

WCG.UpperBound = UpperBound;
Info.Polelist = defaultpl;
Info.Sopt = Sopt;
Info.xopt = xopt;
