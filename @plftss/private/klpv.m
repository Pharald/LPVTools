function K = klpv(P,nmeas,ncont,nri,R,S,J3,L3,gam)
% Construct "central" LPV controller from LMI solution (R,S,J,L,gam)
%
% Ref:
%  1. Explicit Controller Formulas for LMI-based Hinfty Synthesis,
%       by Gahinet, Automatica, 1996.
%  2. A Convex Characterization of Gain-Scheduled Hinfty Controllers,
%       by Apkarian and Gahinet, IEEE TAC, 1995.

% Tolerances
% XXX How to set this?
tol = sqrt(eps);

% Uncertainty block data
nr = sum(nri);
Nblk = length(nri);

% Plant dimensions
% nd is dim of true generalized disturbances (not including param inputs)
% ne is dim of true generalized errors (not including param outputs)
szP = size(P);
ny = nmeas;
nu = ncont;
ne = szP(1) - ny - nr;
nd = szP(2) - nu - nr;

% Form Augmented Plant (See Ref 2)
% Pa contains extra channels for the parameters going to/from K
Pa = [zeros(nr,nr+szP(2)) eye(nr); ...
    zeros(szP(1),nr) P zeros(szP(1),nr); ...
    eye(nr) zeros(nr,nr+szP(2))];

% Pull out state matrices from Pa
% nda = augmented disturbances including channels for parameters
% nea = augmented errors including channels for parameters
% XXX: Assumes controller uses the same # of copies of param as the plant
[A,B,C,D] = ssdata(Pa);
nx = size(A,1);
nda = nd + 2*nr;
nea = ne + 2*nr;
nya = ny + nr;
nua = nu + nr;

B1 = B(:,1:nda);
B2 = B(:,nda+1:end);

C1 = C(1:nea,:);
C2 = C(nea+1:end,:);

D11 = D(1:nea,1:nda);
D12 = D(1:nea,nda+1:end);
D21 = D(nea+1:end,1:nda);

% Step 1: Compute M,N such that MN' = I-RS
% XXX: Assumes controller order is the same as the plant
[U,Sig,V] = svd( eye(nx) - R*S );
Srt = Sig.^0.5;
M = U*Srt;
N = V*Srt;

% Step 2: Compute scaling L
% XXX: Assumes controller uses the same # of copies of param as the plant
% XXX: Can this be implemented without explicitly constructing L?
L1 = cell(Nblk,1);
L2 = cell(Nblk,1);
ptr = 0;
for i=1:Nblk
    % Index for ith block
    idx = ptr+(1:nri(i));
    ptr = idx(end);
    
    % Compute W,Z such that WZ' = I-J3*L3 for i^th block        
    [U,Sig,V] = svd( eye(nri(i)) - J3(idx,idx)*L3(idx,idx) );
    Srt = diag( diag(Sig).^0.5 );
    W = U*Srt;
    Z = V*Srt;
        
    % Solve for i^th block of L1 and L2
    L2{i} = Z';
    L1{i} = -Z'*J3(idx,idx)/(W');
end
if nr>0
    L1 = blkdiag( L1{:} );
    L2 = blkdiag( L2{:} );
    L = [L1 L2; L2' L3];
    Lrt = chol(L);
else
    L = [];
    Lrt = [];
end

% Step 3: Solve static Parrott equation to compute DK such that
%         Ld  - Dcl'*Le*Dcl >0
% where Dcl = D11+D12*DK*D21
gamrt = sqrt(gam);
Ldrt = blkdiag(Lrt,gamrt*eye(nd));
Lert = blkdiag(Lrt,eye(ne)/gamrt);
D11til = (Lert*D11)/Ldrt;
D12til = Lert*D12;
D21til = D21/Ldrt;

[U12,S12,V12] = svd(D12til);
S12d = diag(S12);
r12 = length(find( S12d>tol*S12(1) ));

[U21,S21,V21] = svd(D21til);
S21d = diag(S21);
r21 = length(find( S21d>tol*S21(1) ));

AA = U12(:,r12+1:end)'*D11til*V21(:,r21+1:end);
BB = U12(:,r12+1:end)'*D11til*V21(:,1:r21);
CC = U12(:,1:r12)'*D11til*V21(:,r21+1:end);
DKtil = LocalParrott(AA,BB,CC);
DK = DKtil-U12(:,1:r12)'*D11til*V21(:,1:r21);
DK = V12*diag(1./S12d(1:r12))*DK*diag(1./S21d(1:r21))*U21';

% Step 4: Solve for BKtil
% XXX This assumes the problem is regular: (I-D21*pinv(D21)) C2 = 0
Dcl = D11+D12*DK*D21;
Ld = blkdiag(L,gam*eye(nd));
Lei = blkdiag(inv(L),gam*eye(ne));
DEL = [ Ld -Dcl'; -Dcl Lei];
tmpB = [zeros(nya) D21 zeros(nya,nea); [D21'; zeros(nea,nya)] -DEL ];
vB = -tmpB\[C2; B1'*S; C1+D12*DK*C2];
BKtil = vB(1:nya,:)';

% Step 5: Solve for CKtil
% XXX This assumes the problem is regular: (I-pinv(D12)*D12) B2' = 0
tmpC = [zeros(nua) zeros(nua,nda) D12'; [zeros(nda,nua); D12] -DEL];
vC = -tmpC\[B2'; (B1+B2*DK*D21)'; C1*R];
CKtil = vC(1:nua,:);

% Step 6: Solve for AKtil 
LL1 = [B1+B2*DK*D21 (C1*R+D12*CKtil)'];
LL2 = [S*B1+BKtil*D21 (C1+D12*DK*C2)'];

AKtil = -(A+B2*DK*C2)'-(LL2/DEL)*LL1';

% Step 7: Solve for (AK,BK,CK)
BK = N\(BKtil-S*B2*DK);
CK = (CKtil-DK*C2*R)/(M');
AK = AKtil-S*B2*CK*M'-N*BK*C2*R-S*(A+B2*DK*C2)*R;
AK = (N\AK)/(M');

% Step 6: Pack controller as a PSSLFT (Eq 6.2)
% XXX: Assumes controller uses the same # of copies of param as the plant
K = ss(AK,BK,CK,DK);
%K = lft( K, DELTA );

% Step 7: Undo loopshift to add back D22 term from plant
%K = feedback( K, D22 );

end



%%
% Solve Static Parrott Problem
%  inf_X || [A B; C X] || = max( ||[A B]||, ||[A;C]|| ) := g*
%
% For any g>g*, X = -C* inv(g^2I-A'A) * A'B gives a cost < g.
function [X,gopt] = LocalParrott(A,B,C,g)

if nargin==3
    g = 1;
end

gopt = max( norm([A,B],2), norm([A;C],2) );
if g<=gopt
    error('Increase cost in Parrott Problem.');
end

na = size(A,2);
X = (g^2*eye(na)-A'*A)\(A'*B);
X = -C*X;

end
