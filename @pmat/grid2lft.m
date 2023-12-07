function [B,C,RMSerr] = grid2lft(A,varargin)
% GRID2LFT Transform a PMAT into a PLFTMAT model with polynomial dependence
%
% Transform a grid-based LPV model into a LFT model with polynomial 
% parameter dependence. Use linear regression to fit a polynomial in the 
% parameters to the grid based data.
%
% L = grid2lft(G) fits the elements of the data in G with a linear
% parameter dependence.  G is a PMAT and L is an PLFTMAT.
%
% L = grid2lft(G,N) fits the elements of the data in G with a polynomial
% parameter dependence of order N.
% 
% L = grid2lft(G,NAMES,DEGREEMAT) fits the data in G with a polynomial,
% using a linear combination of monomials specified by the data in
% DEGREEMAT.  NAMES is a 1-by-P cell array of chars, consisting of the
% P names of every independent variable in G.  DEGREEMAT is a D-by-P matrix
% of nonnegative integers, each 1-by-P row corresponding to a monomial,
% defined by the nonnegative exponents associated with each independent
% variable. 
% 
% [L,C] = grid2lft(G,...) returns C, the matrix of polynominal 
% coefficients used in the transformation from PMAT to PLFTMAT. If G is a 
% M-by-N matrix that is being fit with a polynominal with B terms, then C 
% is a M-by-N-by-B double array, in which elements (i,k,:) correspond to 
% (i,k)-th matrix element in G, and elements (:,:,r) correspond to the r-th 
% basis function.
% 
% [L,C,E] = grid2lft(G,...) returns E, the RMS error of the linear fit.
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create PMATs M and M2 with two independent variables x and y. 
%   x = pgrid('x',linspace(-2,2,12),[-1 1]);
%   y = pgrid('y',1:5,[-4 8] );
%   M = [x+y-x*y x;3*y -2*x*y];
%   M2 = sqrt(1+x.^2)*y;
% 
%   % Transform both M and M2 into LFT based LPV objects. Use a polynomial 
%   % containing the factors (1,x,y,x*y) to perform the fitting for M, and 
%   % a polynomial (1,x,y,x^2,x*y,x^2*y) to perform the fitting for M2.
% 
%   % Call grid2lft and specify that the fitting of M should use the 
%   % polynomial (1,x,y,x*y)
%   [Mlft,C,E] = grid2lft(M,{'x','y'},[0 0;1 0;0 1;1 1])
% 
%   % Call grid2lft and specify that the fitting of M2 should use the 
%   % polynomial (1,x,y,x^2,x*y,x^2*y)
%   [M2lft,C2,E2] = grid2lft(M2,{'x','y'},[0 0;1 0;0 1;2 0;1 1;2 1])
% 
% See also pss/grid2lft, lft2grid.

nin = nargin;
narginchk(1,3)

if isa(A,'pgrid')
    A = pmat(A);
elseif ~isa(A,'pmat')
    error('First argument must be a PMAT')    
end
if A.Domain.NumIV == 0
    B = plftss(A.Data);
    RMSerr = 0;
    C = [];
    return
end

IVN = A.Domain.IVName;
LN = length(IVN);
if nin ==1
    Deg = [zeros(1,LN);eye(LN)];
    [B,C,RMSerr] = grid2lft(A,IVN,Deg);
    return
elseif nin ==2
    Ndeg = varargin{1};
    Deg = [];
    for i = 0:Ndeg
        Deg = [Deg;GetDegMat(i,LN)];
    end
    [B,C,RMSerr] = grid2lft(A,IVN,Deg);
    return
end

NameIn = varargin{1};
Deg    = varargin{2};
% XXX Add check to make sure rows in Deg are unique.

[~,idx] = setdiff( NameIn,IVN );
NameIn(idx) = [];
Deg(:,idx) = [];

% Make sure the input contains the names of every parameter in A
[~,idxOrder] = ismember(IVN,NameIn);
if any(idxOrder==0)
    error([' Input must include every parameter that is contained in ' inputname(A)])
end

if isa(varargin{2},'double')
    
    Deg = Deg(:,idxOrder);
    
    % Reorder as [row col AD IV]
    niv = A.Domain.NumIV;    
    nad = numel(size(A))-2;
    AData = permute(A.Data,[1 2 (3+niv:2+niv+nad) (3:2+niv)]);
    
    % Reshape data to put grid point matrix data into individual rows in AData
    sza = size(AData);
    Npts = prod(sza(3+nad:end));
    sza(3+nad:end) = [];
    Nelem = prod(sza);
    AData = reshape(AData,[Nelem,Npts])';   
    
    % Reorder the parameter grid into a Npts x niv matrix, each row
    % defines a grid point, i.e. a single parameter combination on the grid.
    GridCell = cell(niv,1);
    IVD = A.Domain.IVData;
    [GridCell{:}] = ndgrid(IVD{:});
    GridMat = zeros(Npts,niv);
    for i =1:niv
        GridMat(:,i) = GridCell{i}(:);
    end
    
    % Define the monomial matrix 'MonMat'
    Nbas = size(Deg,1);
    MonMat = zeros(Npts,Nbas);
    for i = 1:Nbas
        degi = repmat(Deg(i,:),[Npts,1]);
        MonMat(:,i) = prod(GridMat.^degi,2);
    end
    
    % Solve for the least squares fit
    C = MonMat\AData;
    RMSerr = norm(MonMat*C-AData)/sqrt(numel(AData));
    
    % Reshape C back to correct dimension:
    C = reshape(C',[sza,Nbas]);
    
    % Define TVREAL objects for each parameter in A.
    TVout = cell(niv,1);
    IVRB = A.Domain.IVRateBounds;
    for i = 1:niv
        namei = IVN{i};
        rangei = [IVD{i}(1) IVD{i}(end)];
        %nomi = (rangei(1)+rangei(2))/2;
        ratei = IVRB(i,:);
        %TVout{i} = tvreal(namei,nomi,'Range',rangei,'RateBounds',ratei);
        TVout{i} = tvreal(namei,rangei,ratei);
    end
    
    B = directlft(C,Deg,TVout);
    
elseif iscell(varargin{2})
    %XXX fix Deg = Deg(:,idxOrder); for cell array
    
end

end

function degmat = GetDegMat(deg,nvar)

if isa(deg,'double') && isscalar(deg) && ...
        all(floor(deg)==ceil(deg)) && deg>=0
    
    if nvar ==0
        degmat = zeros(1,0);
    elseif nvar ==1
        % Base case for recursion
        degmat = deg(:);
    else
        % Recursive call for nvar>1
        maxd = max(deg);
        
        % Compute all degmats for nvar-1 vars and degs from 0 up to maxd
        degmatcell = cell(maxd+1,1);
        for i1 =0:maxd
            degmatcell{i1+1} = GetDegMat(i1,nvar-1);
        end
        
        % Stack on degrees for last variable
        degmat = [];
        for i1 = 1:length(deg)
            degi = deg(i1);
            degmati = [];
            for i3=0:degi
                temp = degmatcell{degi-i3+1};
                temp = [temp repmat(i3,[size(temp,1),1])];
                degmati = [degmati;temp];
            end
            
            degmat = [degmat;degmati];
        end
    end
    
else
    error('Degree must be a non-negative integer')
end
end


function L = directlft(C,Deg,TVout)
%XXX Array dimensions are not allowed in this first implementation.

niv = length(TVout);
Nterm = size(Deg,1);
szc  = size(C);
tol = 1e-8; %XXX what should we choose?

% Find constant term in C
constidx = find(sum(Deg,2)==0);
if isempty(constidx)
    C0 = zeros(szc(1:2));
else
    C0 = C(:,:,constidx);
end

% Project C into non-constant rows and colums
E = C;
E(:,:,constidx) = [];
absE = sum(abs(E),3);
ridx = sum(absE,2) > tol;
cidx = sum(absE,1) > tol;
Cproj = C(ridx,cidx,:);
szcproj = size(Cproj(:,:,1));

% Handle constant term
Lproj = zeros(szcproj);
for i = 1:Nterm
    % Grab term data
    idegrow = Deg(i,:);
    ci = Cproj(:,:,i);
    d = sum(idegrow);
    
    if d==0
        % Handle Constant Term
        Lterm = ci;
    else
        % Split matrix into LR through svd based left right decomposition
        [Ui,Si,Vi] = svd(ci);
        r = length(find(diag(Si)>tol));
        Sroot = sqrt(Si(1:r,1:r));
        Li = Ui(:,1:r)*Sroot;
        Ri = Sroot*Vi(:,1:r)';
        
        % Construct a LFT description of the monomial
        sd = diag(ones(d-1,1),-1);
        e1 = [1;zeros(d-1,1)];
        ed = [zeros(d-1,1);1];
        m1 = [sd e1;ed' 0];
        D = [];
        for j = 1:niv
            D = blkdiag(D,TVout{j}*eye(idegrow(j)));
        end
        Lmonom = lft(D,m1);
        
        % Construct LFT description of Li*[monomial*Irxr]*Ri
        Lmonom = Lmonom*eye(r);
        m2 = [zeros(r) Ri;Li zeros(szcproj)];
        Lterm = lft(Lmonom,m2);
    end
    
    % Sum term into projected LFT
    Lproj = Lproj + Lterm;
    
end

% inverse projection
L = plftmat(C0);
L(ridx,cidx) = Lproj;

end





