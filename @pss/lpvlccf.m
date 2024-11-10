function [fact,Ml,Nl,L] = lpvlccf(P,Zb)
% LPVLCCF  Compute left contractive coprime factorization for a PSS
%
% [FACT,Ml,Nl,L] = LPVLCCF(P,Zb) computes a contractive
% coprime factorization of the LPV system P (analogous to a normalized
% coprime factorization for LTI systems).Zb is a BASIS object, which 
% describe the assumed parameter dependence of the solution matrix of
% generalised filtering inequality.
% FACT = LPVLCCF(P,Zb) returns a state-space realization FACT of [Ml,Nl].
% [FACT,Ml,Nl,L] = LPVLCCF(P,Zb) also returns the coprime factors Ml, Nl and
% the gain L separately.

nin = nargin;
narginchk(1, 2);
if nin==1
    Zb = [];
end


% Parse system data:
if isempty(Zb)
    Zb = basis(1,0);
elseif ~isa(Zb,'basis')
    error('Zb must be a BASIS object')
end


% Single balancing transformation
% TODO HP 10/11/24: in the original code we did not balance.
% P0 = P;
% P = lpvbalance(P0);

% Transform system and basis to non-LPVTools format.
[Pdata,RBz,BFz,Pz] = basis2data(P,Zb);

% Get state space data of PMAT system P
[a,b,c,d] = ssdata(P);
nX = size(a,1);


% Dimensions
Pdata = Pdata(:,:,:);
szP = size(Pdata);
if numel(szP) ==2
    szP = [szP 1];
end
ny = szP(1);    % # of measurements
nu = szP(2);    % # of inputs 
nmod = szP(3); % # of model in the Pdata array.
nparz = size(RBz,1); % # of parameters in the basis functions
nbasisz = size(BFz,1); % # of basis functions
Zbivn = Zb.IVName;


R = eye(ny)+d*d';
S = eye(nu)+d'*d; 
Atil = a-b*inv(S)*d'*c;
Ctil = c'*inv(R)*c;


% Check for non-rate bounded case
RateBndFlag = 1;
if (nbasisz==1 && nparz==0) 
    RateBndFlag = 0; 
else
    % User has input non-constant basis functions, intends Rate-bounded syn
    % Make sure that all system parameters used in basis functios have 
    % finite ratebounds.
    RBcheck = isfinite(RBz);
    RBresult = all(RBcheck(:));
    if ~RBresult
        error(['While attempting rate-bounded coprime factorization for system P. '...
               'BASIS functions Zb depend on a parameter that '...
               'has non-finite rate-bounds in P.'])
    end    
end

% Create LMI variables
setlmis([]);
Zdec = zeros(nbasisz,1);
for k=1:nbasisz
    Zdec(k) = lmivar(1,[nX 1]);
end
[Wdec,~,Ws] = lmivar(1,[nX 1]);


% Define LMIs
cnt = 1;
for i=1:nmod
    % Transformed Data
    Si = S.Data(:,:,i);   
    Atili = Atil.Data(:,:,i);
    Bi = b.Data(:,:,i);
    Ctili = Ctil.Data(:,:,i);
    
    if RateBndFlag
        Pzi = Pz(:,:,i);          
    end
    BFzi = BFz(:,:,i); 

    for q = 1:2^nparz  % number of variables appearing in Z's basis fcns
        RBV = RBz(:,1);
        tmp = dec2bin(q-1,nparz);
        idx = find(tmp=='1');
        RBV(idx) = RBz(idx,2);

        for k=1:nbasisz
            for j=1:nparz
                % TODO HP 10/11/24: check for the sign 
                lmiterm([cnt 1 1 Zdec(k)],RBV(j)*Pzi(k,j),1);
            end
            lmiterm([cnt 1 1 Zdec(k)],1,Atili*BFzi(k),'s');
            lmiterm([cnt 1 2 Zdec(k)],1,Bi*BFzi(k));
        end

        lmiterm([cnt 1 1 0],-Ctili) 
        lmiterm([cnt 2 2 0],-Si)
        
        cnt = cnt+1;
    end
    if RateBndFlag || i == 1
        lmiterm([-cnt 1 1 Wdec],1,1)
        lmiterm([-cnt 1 2 0],eye(nX))
        for k=1:nbasisz
            lmiterm([-cnt 2 2 Zdec(k)],1,BFzi(k))  
        end
         ZMIN = 1e-9; 
         lmiterm([cnt 1 1 0],ZMIN*eye(nX)) 
         lmiterm([cnt 2 2 0],ZMIN*eye(nX))
        cnt = cnt+1;
    end

end



% Get LMI Options
% TODO HP 10/11/2024: Shall we include options object in the code?

% Default settings for LMI Lab
LMIopt = zeros(5,1);
LMIopt(2) = 350;  % Max # of iters for rate bounded syn
LMIopt(5) = 1;        % Toggle display

% Solve LMI
lmisys = getlmis;
FeasFlag = 1;

% Create objective function:
ndec = decnbr(lmisys); 
cobj = zeros(ndec,1);
cobj(diag(Ws)) = 1; % minimize Trace of W

[copt,xopt] = mincx(lmisys,cobj,LMIopt);
if isempty(copt)
    FeasFlag = 0;
end

% Handle Infeasible LMI Case
if ~FeasFlag
    fact = [];
    Ml =[];
    Nl = [];
    L =[];
    return;
end

Z = LOCALdec2pmat(nX,lmisys,xopt,Zdec,P.Domain,BFz);


% Coprime Factorization State Matrices reconstruction
L = -(b*d'+Z\c')/R; 

Afact = a+L*c;
Bfact = [L, b+L*d];
Cfact = R^(-0.5)*c;
Dfact = [R^(-0.5), R^(-0.5)*d];

Ml = ss(Afact,Bfact(:,1:ny),Cfact,Dfact(:,1:ny));
Nl = ss(Afact,Bfact(:,ny+1:end),Cfact,Dfact(:,ny+1:end));
fact =ss(Afact,Bfact,Cfact,Dfact);

end 


function Y = LOCALdec2pmat(nX,lmisys,xopt,Ydec,Domain,BFy)
% helper function that transforms the decision variables into PMAT objects.
% Y = zeros(nX,nX);
% for k=1:nbasisy
%     Y = Y + BFypmat(k)*dec2mat(lmisys,xopt,Ydec(k));
% end
% Y = lpvsplit(Y,sys.DomainPrivate); XXX Why this?
        
        

nbasisy = size(BFy,1);
nmod = size(BFy,3);

Ymat = zeros(nX,nX,nbasisy);
for k=1:nbasisy
    Ymat(:,:,k) = dec2mat(lmisys,xopt,Ydec(k));
end
    
Yopt = zeros(nX,nX,nmod);
for i1 =1:nmod
    for i2 = 1:nbasisy
        Yopt(:,:,i1) = Yopt(:,:,i1)+BFy(i2,1,i1)*Ymat(:,:,i2);
    end
end
Yopt = reshape(Yopt,[nX;nX;Domain.LIVData]');
if isempty(Yopt)
    Y = pmat([]);
else
    Y = pmat(Yopt,Domain);
end

end
