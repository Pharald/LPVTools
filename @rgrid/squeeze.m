function [E,keepidx,rmidx]= squeeze(R)
% SQUEEZE  Removes singleton array dimensions
%
% SQUEEZE(R) removes all singleton array dimensions from R.
%
% See also: lpvelimiv.

% UNDOCUMENTED FOR USERS

%out = lpvelimiv(m);

% Get rgrid data
ivD = R.IVData;
ivN = R.IVName;
ivRB = R.IVRateBounds;

% Find locations of array dimensions and indep vars
adidx =strncmp(R.ArrayName,R.IVName,length(R.ArrayName));
ividx = ~adidx;

% Remove singleton array dimensions
rmidx = find(adidx & R.LIVData==1);
adidx = find(adidx);

adidx( R.LIVData(adidx)==1 ) = [];

% Construct rgrid object without single array dimensions
E = rgrid( ivN(ividx), ivD(ividx) , ivRB(ividx,:) );
E = rgrid( E, R.LIVData(adidx) );

% Construct numeric indices of retained indep vars and array dims
keepidx = [find(ividx); adidx];

