%% LPVSUBS - Evaluate LPV object and demote to LTI
%
%  
%% Syntax
% 
%    B = lpvsubs(A,NAME,VALUES)
%    B = lpvsubs(A,NAME,VALUES,METHOD)
%
%% Description
%
% |lpvsubs|  evaluates a grid-based LPV object at points in the domain, 
% and demotes the results to a standard LTI object |double|, |ss|, or |frd|.
%
% |B = lpvsubs(A,NAME,VALUES,METHOD)| evaluates a grid-based LPV object |A| 
% at the domain points specified by the |NAME| / |VALUES| pair.  |NAME| is 
% an N-by-1 cell array of characters. |NAME| must contain all of the 
% independent variable names in |A|, namely |A.IVName|, but may also contain 
% others but they will be ignored. |VALUES| is an N-by-NPTS double array of 
% the corresponding values.  |B| is a double array, with row and column 
% dimensions from |A|.  The 3rd dimension of |B| is NPTS.  |METHOD| is an 
% optional input that specifies the interpolation method and can be: 
% |'nearest'|, |'linear'|, |'spline'|, or 
% |'cubic'|. The default is |'linear'|.