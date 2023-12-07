function [gam,X,Info] = lpvwcgain(P,Xb,polelist)
% LPVWCGAIN  Worst-case bound on induced L2 norm for PSS systems.
% 
% [GAM,X,INFO] = LPVWCGAIN(P,Xb) computes the upper-bound on the worst-case
% induced L2 norm of the PSS P. "Worst-case" refers to all the modeled
% uncertainty and parameter-dependence (including parameter rate-bounds). 
% The upper bound GAM and a parameter-varying matrix X are computed to 
% satisfy the induced L2 norm linear matrix inequality (LMI) condition. 
% Xb is a BASIS object that defines the basis functions which describe the 
% assumed parameter dependence of X. INFO is a structure containing 
% additional information about the solution to the LMI.
%
% [GAM,X,INFO] = LPVWCGAIN(P) computes the upper-bound on the worst-case
% induced L2 norm of the PSS P, assuming no rate-bounds on the independent 
% variables of P.  The upper bound GAM and a constant (parameter independent) 
% matrix X are computed to satisfy the induced L2 norm linear matrix 
% inequality (LMI) condition. 
%
% [GAM,X,INFO] = LPVWCGAIN(P,...,POLELIST) allows the user to define 
% weighting functions for the Integral Quadratic Constraints (IQC) used 
% to bound the uncertainty when formulating the LMI to be solved. POLELIST 
% is a 1xN DOUBLE row vector, of negative values. Each value in POLELIST 
% corresponds to a pole of a stable transfer function that is used as a 
% weight on all signals in the IQCs. A default POLELIST, with three pole 
% values, is used when a POLELIST is not supplied by the user. The three 
% pole values are selected automatically from the frequency range of the 
% system dynamics.
%
% See also: lpvnorm, norm, wcgain.


% Parse Inputs
nin = nargin;
narginchk(1,3)
if nin==1
    Xb = [];
    polelist = [];
elseif nin ==2
    polelist = [];
end

if isempty(Xb)
   Xb = basis(1,0); 
end

% SS array of data. G is the M of the M-Delta form:
% G: Array of LTI systems.
if isuncertain(P)
    [G,unc,blkstruct,~] = lftdata(P);
    fn = fieldnames(unc.Uncertainty);
    % mu block structure:
    % Ublksize: Nblk-by-2 matrix of positive integers. Each row is [R,C] where
    %           R and C are the number of rows and columns of the uncertainty.
    Ublksize = muBlkHelper(P);
    % Process Ublksize
    uidx = find(Ublksize(:,2)==0);
    Ublksize(uidx,:) = abs(Ublksize(uidx,1))*[1 1];
    
    % Get automatically generated frequency grid and pole list:    
    [~,defaultpl] = polefreqauto(G.Data,fn,polelist,'nogenfresp');
else
    blkstruct = [];
    Ublksize  = [];
    defaultpl = [];
    fn = cell(0,1);
end

% Create structure of IQCs:
% IQCS:  Nblk-by-1 structured array specifying the IQCs in factorized
%        form as (psi,M).  The fields of IQCS are:
%        IQCS(i).psi: Ni-by-1 cell array of LTI systems
%        IQCS(i).M: Ni-by-1 cell array of real, symmetric matrices.

IQCS = struct('psi',[],'M',[]);

for k = 1:numel(fn)
    uname = fn{k};
    polenum = numel(defaultpl.(uname));
    
    if isequal(blkstruct(k).Type,'ureal')
        % UREAL
        r = abs(Ublksize(k,1));
        iqcpsi = cell(r*polenum,1);
        iqcM   = cell(r*polenum,1);
        cnt = 1;
        for i = 1:polenum
            pol = abs(defaultpl.(uname)(i));
            %XXX Assumes only first order stable poles.
            if isinf(pol)
                d = 1;
            else
                d = tf(pol,[1 pol]);
            end
            %XXX Assuming Y's in eq 27 in Megretski and Ranzter paper are zero
            for rk = 1:r
                dmat = ss(eye(r));
                dmat(rk,rk) = d;
                iqcpsi{cnt} = blkdiag(dmat,dmat);
                iqcM{cnt}   = blkdiag(eye(r),-eye(r));
                cnt = cnt+1;
            end
        end
    elseif isequal(blkstruct(k).Type,'ultidyn')
        % ULTI
        iqcpsi = cell(polenum,1);
        iqcM   = cell(polenum,1);
        r = Ublksize(k,1);
        c = Ublksize(k,2);
        for i = 1:polenum           
            pol = abs(defaultpl.(uname)(i));
            if isinf(pol)
                d = 1;
            else
                %XXX Assumes only first order stable poles.
                d = tf(pol,[1 pol]);                
            end
            iqcpsi{i} = blkdiag(d*eye(r),d*eye(c));
            iqcM{i}   = blkdiag(eye(r),-eye(c));
        end
    else
        error('Not implemented as of yet')
    end
    
    IQCS(k).psi = iqcpsi;
    IQCS(k).M   = iqcM;
end

[Gdata,RBx,BFdatax,Pdatax] = basis2data(G,Xb);

gam = iqcgridengine(Gdata,BFdatax,Pdatax,RBx,Ublksize,IQCS);
X = [];
Info = [];