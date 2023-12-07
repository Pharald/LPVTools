%% PLFTMAT - Parameter-varying matrix in LFT framework
%
%  
%% Syntax
%
% |M = plft(Data,RateBounds)|
%
%% Description
%
% |M = plft(Data,RateBounds)| creates a parameter-varying matrix.
% |Data| is a |umat|. |RateBounds| is a N-by-2 cell array listing the
% rate bound information for each independent variable in the |plftmat|.  
% |RateBounds{i,1}| is the character string name of the i-th independent 
% variable and |RateBounds{i,2}| is a sorted real vector of form [Low, High] 
% specifying its rate bounds.  |RateBounds| must only contain names of |ureal| 
% objects that exist in |Data| and this indicates that those |ureal| are actually 
% |tvreal| representing the independent variables.
% 
%% Example


% Create a 2-by-2 matrix M that depends on TVREAL a.
a = tvreal('a',[-2 2],[-1 1]);   
M = [1, a;a^2, -a]