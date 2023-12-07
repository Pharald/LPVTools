%% PLFTSS - Parameter-varying state-space model in LFT framework
%
%  
%% Syntax
%
% |M = plftss(Data,RateBounds)|
%
%% Description
%
% |M = plftss(Data,RateBounds)| creates a parameter-varying state-space model.
% |Data| is a |uss|. |RateBounds| is a N-by-2 cell array listing the
% rate bound information for each independent variable in the |plftss|.  
% |RateBounds{i,1}| is the character string name of the i-th independent 
% variable and |RateBounds{i,2}| is a sorted real vector of form [Low, High] 
% specifying its rate bounds.  |RateBounds| must only contain names of |ureal| 
% objects that exist in |Data| and this indicates that these |ureal| are actually 
% |tvreal| representing the independent variables.
% 
%% Example


% Create a 1-by-1 state-space model S that depends on tvreal b.
b = tvreal('b',[2 20]);
S = ss(-b,b,1,0)
