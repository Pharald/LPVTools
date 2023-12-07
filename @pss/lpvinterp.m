function Si = lpvinterp(G,varargin)
% LPVINTERP  Interpolate a PSS
%
% B = LPVINTERP(A,NAME1,VALUES1,NAME2,VALUES2,...) interpolates a PSS A on
% the domain specified by the NAME/VALUES pairs. Each VALUE is a vector of 
% values at which to interpolate A along the domain direction of the
% corresponding NAME. If an independent variable of A is not listed
% in the inputs, then all values along this domain direction are retained
% in the output B.
%
% B = LPVINTERP(A,NAME,VALUES) is an alternative syntax. NAME is an N-by-1
% cell array of characters and VALUES is an N-by-1 cell array of values.
%
% B = LPVINTERP(A,NAME1,VALUES1,NAME2,VALUES2,....,METHOD) includes a final 
% input argument called METHOD, which specifies the interpolation method to 
% be used. METHOD can be: 'nearest', 'linear', 'spline', or 'cubic'. The
% default is 'linear'. 
%
% See also: lpvsubs, lpvsplit, lpvsample.

% Get state-space data (state-dim is constant)
[A,B,C,D] = ssdata(G.DataPrivate);
szG = size(G);
nS = size(A,1);

% Stacked state-space matrices and use PMAT/LPVINTERP
M = pmat([A B;C D],G.DomainPrivate);
Mi = lpvinterp(M,varargin{:});

% Unstack PMAT state-space matrices, and call constructor
Ai = reshape(Mi(1:nS,1:nS,:),[nS nS szG(3:end)]);
Bi = reshape(Mi(1:nS,nS+1:end,:),[nS szG(2) szG(3:end)]);
Ci = reshape(Mi(nS+1:end,1:nS,:),[szG(1) nS szG(3:end)]);
Di = reshape(Mi(nS+1:end,nS+1:end,:),[szG(1) szG(2) szG(3:end)]);

% Pack up data for output (currently delays are not handled)
% SS Syntax inherits the following generic properties from Data:
%   InputDelay, InputGroup, InputName, Notes, OutputDelay,
%   OutputGroup, OutputName, Ts, UserData
% The following properties are not documented as 'generic' but
% also seem to be set in the SS syntax below:
%   Name, TimeUnit, InputUnit, OutputUnit
Si = pss(Ai,Bi,Ci,Di,G.DataPrivate);

% Set non-generic properties
% TODO PJS 5/1/2011: Are there other properties to be set?
% Si = set(Si,'Scaled',Data.Scaled, ...
%     'StateName',Data.StateName,'StateUnit',Data.StateUnit);

% % Create output variable
% out = m;
% out.DataPrivate = Si;
% out.DomainPrivate = Ai.DomainPrivate;


