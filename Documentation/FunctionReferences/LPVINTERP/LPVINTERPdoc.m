%% LPVINTERP - Interpolate a grid-based LPV object
%
%  
%% Syntax
%
%    B = lpvinterp(A,NAME,VALUES)
%    B = lpvinterp(A,NAME1,VALUES1,NAME2,VALUES2,...)
%    B = lpvinterp(A,NAME1,VALUES1,NAME2,VALUES2,....,METHOD)
%
%% Description  
%
% |B = lpvinterp(A,NAME1,VALUES1,NAME2,VALUES2,...)| interpolates a 
% grid-based LPV system |A| on the domain specified by the |NAME| / |VALUES| 
% pairs. Each |VALUE| is a vector of values at which to interpolate |A| 
% along the domain direction of the corresponding |NAME|. If an independent 
% variable of |A| is not listed in the inputs, then all values along this 
% domain direction are retained in the output |B|.
%
% |B = lpvinterp(A,NAME,VALUES)| is an alternative syntax. |NAME| is an N-by-1
% cell array of characters and |VALUES| is an N-by-1 cell array of values.
%
% |B = lpvinterp(A,NAME1,VALUES1,NAME2,VALUES2,....,METHOD)| includes a final 
% input argument called |METHOD|, which specifies the interpolation method to 
% be used. |METHOD| can be: |'nearest'|, |'linear'|, |'spline'|, or |'cubic'|. 
% The default is |'linear'|. 






