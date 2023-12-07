%% LFT2GRID - Transfrom LFT into grid-based LPV object
%
%  
%% Syntax
%
%   M = lft2grid(L)
%   M = lft2grid(L,N)
%   M = lft2grid(L,DOMAIN)
%
%% Description
%
% Transform a |tvreal|, |plftmat| or |plftss| object into a grid based 
% LPV object |pmat|, |pss|, |upmat|, or |upss|. The transformation is 
% performed by evaluating the PLFT object at a grid of parameter values.
%
% |M = lft2grid(L)| evaluates |L| at 10 values of each independent parameter,
% sampled uniformly from the range of each parameter.
%
% |M = lft2grid(L,N)| evaluates |L| at |N| values of each independent parameter,
% sampled uniformly from the range of each parameter.
%
% |M = lft2grid(L,DOMAIN)| evaluates |L| at each point containted in |DOMAIN|. 
% |DOMAIN| is an |rgrid| object that must contain the same independent 
% variables as |L|.
