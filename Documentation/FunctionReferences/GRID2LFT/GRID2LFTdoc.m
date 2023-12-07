%% GRID2LFT - Transform a grid-based LPV model into LFT
% 
%  
%% Syntax
% 
%    [L,C,E] = grid2lft(G)
%    [L,C,E] = grid2lft(G,N)
%    [L,C,E] = grid2lft(G,NAMES,DEGREEMAT) 
% 
%% Description
%
% Transform a grid-based LPV model into a LFT model  
% with polynomial parameter dependence. Use linear regression to fit a 
% polynomial in the parameters to the grid based data.
%
% |L = grid2lft(G)| fits the elements of the matrix or state-space data in 
% |G| with a linear parameter dependence. |G| is a grid based LPV model 
% (e.g. pmat or upss) and |L| is an LFT based LPV model (e.g. plftmat, plftss).
% 
% |L = grid2lft(G,N)| fits the elements of the matrix or state-space data 
% in |G| with polynomial parameter dependence of degree |N|.
% 
% |L = grid2lft(G,NAMES,DEGREEMAT)| fits the matrix or state-space data in 
% |G| with polynomials, using a linear combination of monomials specified by the
% data in |DEGREEMAT|.  |NAMES| is a 1-by-P cell array of chars, consisting of
% the P names of every independent variable in |G|.  |DEGREEMAT| is a D-by-P
% matrix of nonnegative integers, each 1-by-P row corresponding to a
% monomial, defined by the nonnegative exponents associated with each
% independent variable.
%
% |[L,C] = grid2lft(G,...)| returns|C|, the matrix of polynominal 
% coefficients used in the transformation from grid-LPV to LFT. If |G| is a 
% M-by-N matrix that is being fit with a polynominal with B terms, then |C| 
% is a M-by-N-by-B |double| array, in which elements (i,k,:) correspond to 
% (i,k)-th matrix element in |G|, and elements (:,:,r) correspond to the r-th 
% basis function.
% 
% |[L,C,E] = grid2lft(G,...)| returns |E|, the root mean square error 
% of the linear fit.

%% Example: Transform a PMAT to a PLFTMAT

% Create PMATs M and M2 with two independent variables x and y. 
x = pgrid('x',linspace(-2,2,12),[-1 1]);
y = pgrid('y',1:5,[-4 8] );
M = [x+y-x*y x;3*y -2*x*y];
M2 = sqrt(1+x.^2)*y;

% Transform both M and M2 into LFT based LPV objects. Use a polynomial 
% containing the factors (1,x,y,x*y) to perform the fitting for M, and 
% a polynomial (1,x,y,x^2,x*y,x^2*y) to perform the fitting for M2.

% Call grid2lft and specify that the fitting of M should use the 
% polynomial (1,x,y,x*y)
[Mlft,C,E] = grid2lft(M,{'x','y'},[0 0;1 0;0 1;1 1]);

%% 
Mlft

%%
C

%%
E

% Call grid2lft and specify that the fitting of M2 should use the 
% polynomial (1,x,y,x^2,x*y,x^2*y)
[M2lft,C2,E2] = grid2lft(M2,{'x','y'},[0 0;1 0;0 1;2 0;1 1;2 1]);

%% 
M2lft

%%
C2

%% 
E2

%% Example: Transform a UPSS to a PLFTSS

% Create UPSS M that depends on two independent variables x and y. 
x = pgrid('x',linspace(-2,2,12),[-1 1]);
y = pgrid('y',1:5,[-4 8] );
u = ureal('u',1);
M = ss(x+y-x*y*u, x+3*y,-2*x*y,pmat(0));

% Transform M into a LFT based LPV object. Use a polynomial containing
% the factors (1,x,y,x^2,x*y,x^2*y) to perform the fitting for M.

% Call grid2lft and specify that the fitting of M should use the 
% polynomial (1,x,y,x^2,x*y,x^2*y)
[Mlft,C] = grid2lft(M,{'x','y'},[0 0;1 0;0 1;2 0;1 1;2 1]);

%%
Mlft

%%
C


