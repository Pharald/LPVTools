function out = lpvinterp(m,varargin)
% LPVINTERP  Interpolate a UPFRD
%
% B = LPVINTERP(A,NAME1,VALUES1,NAME2,VALUES2,...) interpolates a UPFRD A on
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

% Get response data
Data = m.Data;
Response = Data.ResponseData;
Freqs = Data.Frequency;
Domain = m.Domain;

% Response is ordered as [row col Freqs IV AD]. Must reorder it 
% as [row col IV AD Freqs]

% ord is the desired order for permute
sz = size(m.DataPrivate.ResponseData);
ord = 1:numel(sz)-2;
ord(1) = [];      
ord = [ord 1];

% Reorder data as [row col IV AD Freqs]
Response = permute(Response,ord);

% Promote Response to UPMAT
Response = upmat(Response,Domain);

% Interpolate response data using UPMAT/LPVINTERP
Ri = lpvinterp(Response,varargin{:});


% Pack up data for output (TODO: currently delays are not handled)

% Reorder as [row col Freqs IV AD].
RDi = Ri.Data;
RDi = ipermute(RDi,ord);
% TODO PJS 5/1/2011: Si inherits 'generic' LTI properties from Data.
%   Do any other properties need to be set?

Si = frd(RDi,Freqs,Data.Ts);
Si.Unit = Data.Unit;

out = upfrd(Si,Ri.Domain);



% OLD CODE
% 
% % Get response data data
% Data = m.DataPrivate;
% Response = Data.ResponseData;
% Freqs = Data.Frequency
% Domain = m.DomainPrivate;
% 
% % Promote Response PMAT
% Response = pmat(Response,Domain);
% 
% % Interpolate response data using PMAT/LPVINTERP
% Ri = lpvinterp(Response,varargin{:});
% 
% % Pack up data for output (currently delays are not handled)
% % TODO PJS 5/1/2011: Si inherits 'generic' LTI properties from Data.
% %   Do any other properties need to be set?
% Si = ufrd(Ri.DataPrivate,Freqs,Data);
% 
% % Create output variable
% out = m;
% out.DataPrivate = Si;
% out.DomainPrivate = Ri.DomainPrivate;

