%% DOMUNION - Map LPV objects on a common domain
%
%  
%% Syntax
%
%    [Aext,Bext]=domunion(A,B)
%    [A1ext,...,ANext] = domunion(A1,...,AN)
%
%% Description
%
% Let |A| and |B| be grid-based LPV objects.  If |A| depends on independent 
% variables (X,Y) and |B| depends on independent variables (X,Z) then
% |[Aext,Bext]=domunion(A,B)| returns grid-based LPV objects  Aext and Bext, 
% of the same class, that have a common domain with independent variables 
% (X,Y,Z). |Aext| evaluated at point (x,y,z) is given by |A| evaluated at 
% (x,y). |Bext| evaluated at point (x,y,z) is given by |B| evaluated at (x,z).
%
% Given grid-based LPV objects |A1|,|A2|,...,|AN|, the syntax
%   |[A1ext,...,ANext] = domunion(A1,...,AN)|
% constructs |A1ext|, |A2ext|,..., |ANext| that are defined on a common domain.           

