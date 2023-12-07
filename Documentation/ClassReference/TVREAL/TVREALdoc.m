%% TVREAL - Time-varying real parameter
%
%  
%% Syntax
%
%   A = tvreal(Name,Range)
%   A = tvreal(Name,Range,RateBounds)
%
%% Description
%
% |A = tvreal(Name,Range,RateBounds)| specifies the name, range, 
% and ratebounds for a time-varying real parameter, where
% * |Name| is a character string specifying the name
% * |Range| is a 1-by-2 row vector specifying the lower and upper limits
%   for the |tvreal|.
% * |RateBounds| is a 1-by-2 row vector specifying lower and upper
% bounds on the derivative of the parameter with respect to time.
% Set |RateBounds(1)=-inf| and/or  |RateBounds(2)=+inf| to denote an
% unbounded rate of change. |RateBounds| are optional, and a two argument 
% call: |A = tvreal(Name,GridData)| 
% will set them to a default value of |[-inf,+inf]|.
% 
%% Example
% 

% Create a tvreal "a" which has range [-2 2] and ratebounds [-1 1].
a = tvreal('a',[-2 2],[-1 1])

%%

% Create a tvreal "b" which has range [2 20] and default 
% ratebounds [-inf inf].
b = tvreal('b',[2 20])

%%

% Use tvreal as a building block: Create a 2-by-2 matrix that depends 
% on the tvreal "a".
M = [1, a;a^2, -a]

%%

% Use tvreal as a building block to build a plftss: Create a 1-by-1 
% state-space model that depends on tvreals "a" and "b".
S = ss(-a,b,1,0)