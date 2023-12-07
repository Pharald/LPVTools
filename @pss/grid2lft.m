function [G,F] = grid2lft(S,varargin)
% GRID2LFT Transform a PSS into a PLFTSS model.
%
% Transform a grid-based LPV model into a LFT model with polynomial 
% parameter dependence. Use linear regression to fit a polynomial in the 
% parameters to the grid based data.
%
% L = grid2lft(G) fits the elements of the state-space data in G with a
% linear parameter dependence. G is a PSS and L is an PLFTSS.
% 
% L = grid2lft(G,N) fits the elements of the state-space data in G with
% polynomial parameter dependence of degree N.
% 
% L = grid2lft(G,NAMES,DEGREEMAT) fits the state-space data in G with
% polynomials, using a linear combination of monomials specified by the
% data in DEGREEMAT.  NAMES is a 1-by-P cell array of chars, consisting of
% the P names of every independent variable in G.  DEGREEMAT is a D-by-P
% matrix of nonnegative integers, each 1-by-P row corresponding to a
% monomial, defined by the nonnegative exponents associated with each
% independent variable.
%
% [L,C] = grid2lft(G,...) returns C, the matrix of polynominal 
% coefficients used in the transformation from PSS to PLFTSS. If G is a 
% M-by-N matrix that is being fit with a polynominal with B terms, then C 
% is a M-by-N-by-B double array, in which elements (i,k,:) correspond to 
% (i,k)-th matrix element in G, and elements (:,:,r) correspond to the r-th 
% basis function.
% 
%   % EXAMPLE: (CUT/PASTE)
%   % Create PSS M that depends on two independent variables x and y. 
%   x = pgrid('x',linspace(-2,2,12),[-1 1]);
%   y = pgrid('y',1:5,[-4 8] );
%   M = ss(x+y-x*y, x+3*y,-2*x*y,0);
% 
%   % Transform M into a LFT based LPV object. Use a polynomial containing
%   % the factors (1,x,y,x^2,x*y,x^2*y) to perform the fitting for M.
% 
%   % Call grid2lft and specify that the fitting of M should use the 
%   % polynomial (1,x,y,x^2,x*y,x^2*y)
%   [Mlft,C] = grid2lft(M,{'x','y'},[0 0;1 0;0 1;2 0;1 1;2 1])
% 
% See also pmat/grid2lft, lft2grid.

nin = nargin;

narginchk(1,3)

if ~isa(S,'pss')
    error('First argument must be a PSS')
end
if S.Domain.NumIV == 0
    G = plftss(S.Data);
    F = [];
    return
end

[A,B,C,D] = ssdata(S);
M = [A B;C D];
[Gmat,F] = grid2lft(M,varargin{:});

% Convert matrix back to a system
Nx = size(A,1);
Ag = Gmat(1:Nx,1:Nx);
Bg = Gmat(1:Nx,Nx+1:end);
Cg = Gmat(Nx+1:end,1:Nx);
Dg = Gmat(Nx+1:end,Nx+1:end);

G = ss(Ag,Bg,Cg,Dg);