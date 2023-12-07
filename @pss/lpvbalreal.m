function [sysb,G,T,Ti,Wo,Wc] = lpvbalreal(sys,solflag)
% LPVBALREAL   Gramian-based balancing for PSS objects
%
% SYSB = LPVBALREAL(SYS) computes a balanced realization of the 
% parameter varying system SYS. SYSB is dervied by computing a single 
% balancing transformation for SYS and applying it at every point in its 
% domain. 
% 
% [SYSB,G] = LPVBALREAL(SYS) returns G, a vector of singular values 
% describing the input-output mapping of SYSB (comparable to Hankel 
% singular values for LTI systems).
%
% [SYSB,G,T,Ti] = lpvbalreal(SYS) returnes the balancing transformation T, and
% its inverse Ti. 
%
% [SYSB,G,T,Ti] = LPVBALREAL(SYS,...,INVERSE) provides an option to use a 
% alternative implementation of the algorithm which computes the balancing 
% transformation. If INVERT is True the alternative formulation is used. It
% can improve the accuracy of the solution for certain systems. The default 
% implementation assumes INVERT=FALSE.
%
% See also: balreal, lpvgram, lpvbalancmr, lpvncfmr.


% XXX Redo help? input stopping conditions and solvingoptions

if nargin == 1
    solflag = true;
end

[A,B,C,D] = ssdata(sys);
nX = size(A,1);
% Define stopping conditions:
maxiter = 10;
reltol = 1e-3;

% Implements variation of the algorithm (7.7.1) presented in G. Woods Thesis
Wo = eye(nX);
% % Compute Initial Weight to use:
% sz = privatesize(sys);
% npts = prod(sz(3:end));
% WoLTI = gram(sys,'o');
% Wo = sum(WoLTI.Data(:,:,:),3);
go = 1;
k = 1;
J = zeros(maxiter,1);
while go == 1
    
    if solflag
        % Compute Gramians
        Wc = lpvgram(sys,'c',Wo,solflag); 
        Wo = lpvgram(sys,'o',Wc,solflag);
    else
        % Note Wood's remarks on page 169 of his PhD thesis: By inverting the
        % matrices in the solution, you can improve the numerical tractability    
        % Compute Gramians
        Wc = lpvgram(sys,'c',inv(Wo),solflag); 
        Wo = lpvgram(sys,'o',inv(Wc),solflag);
    end
    
    % Compute cost function
    J(k) = trace(Wc*Wo);
    
    % Check stopping conditions  
    if k == 1
        reldiff = 2*reltol;
    else
        reldiff = abs( (J(k)-J(k-1))/J(k) );
    end
    if (k>=maxiter) || (reldiff<reltol)
        go = 0;        
    end
    k = k+1; 
end

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

% Compute balanced system:
sysb = ss2ss(sys,T);

% Compute Hankel singular values:
G = sqrt(diag(So));


function [V,S] = LOCALschursort(W)
    [V,S] = schur(W);
    [Ssort,idx] = sort(diag(S),'descend');
    S = diag(Ssort);
    V = V(:,idx);        


