function [G,F] = grid2lft(S,varargin)
% GRID2LFT Transform a UPSS into a uncertain PLFTSS model.
%
% Transform a grid-based uncertain LPV model into a uncertain LFT model  
% with polynomial parameter dependence. Use linear regression to fit a 
% polynomial in the parameters to the grid based data.
%
% L = grid2lft(G) fits the elements of the state-space data in G with a
% linear parameter dependence. G is a UPSS and L is an PLFTSS.
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
% coefficients used in the transformation from UPSS to PLFTSS. If G is a 
% M-by-N matrix that is being fit with a polynominal with B terms, then C 
% is a M-by-N-by-B double array, in which elements (i,k,:) correspond to 
% (i,k)-th matrix element in G, and elements (:,:,r) correspond to the r-th 
% basis function.
% 
%   % EXAMPLE: (CUT/PASTE)
%   % Create UPSS M that depends on two independent variables x and y. 
%   x = pgrid('x',linspace(-2,2,12),[-1 1]);
%   y = pgrid('y',1:5,[-4 8] );
%   u = ureal('u',1);
%   M = ss(x+y-x*y*u, x+3*y,-2*x*y,pmat(0));
% 
%   % Transform M into a LFT based LPV object. Use a polynomial containing
%   % the factors (1,x,y,x^2,x*y,x^2*y) to perform the fitting for M.
% 
%   % Call grid2lft and specify that the fitting of M should use the 
%   % polynomial (1,x,y,x^2,x*y,x^2*y)
%   [Mlft,C] = grid2lft(M,{'x','y'},[0 0;1 0;0 1;2 0;1 1;2 1])
%
% See also upmat/grid2lft, lft2grid

nin = nargin;

narginchk(1,3)

if ~isa(S,'upss')
    error('First argument must be a UPSS')
end
if S.Domain.NumIV == 0
    G = plftss(S.Data);
    F = [];
    return
end

[M,Delta] = lftdata(S);
[Gmat,F] = grid2lft(M,varargin{:});
G = lft(Delta,Gmat);







