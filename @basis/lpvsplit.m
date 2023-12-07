function out = lpvsplit(m,varargin)
% LPVSPLIT  Extract BASIS data based on independent variable range
%
% B = LPVSPLIT(A,NAME1,RANGE1,NAME2,RANGE2,....) extracts data from the BASIS
% A on the domain specified by the NAME/RANGE pairs. Each RANGE is a 1-by-2 
% row vector [min, max], that specifies the values of the independent 
% variable to be extracted along the domain direction NAME. The data at all
% points in A.Parameter."NAME".GridData which lies inside RANGE is extracted. 
% If RANGE is a scalar then the data is extracted where the variable NAME 
% is exactly equal to RANGE. If an independent variable of A is not listed 
% in the inputs then all values along this domain direction are retained in B.
%
% B = LPVSPLIT(A,NAME1,INDEX1,NAME2,INDEX2,....,'index') extracts data from
% the BASIS A on the domain specified by the NAME/INDEX pairs. Each INDEX
% is a vector of integers or a logical array that specifies the indices
% of A.BasisFunction.Parameter."NAME".GridData to be extracted.
%
% B = LPVSPLIT(A,NAME1,VALUES1,NAME2,VALUES2,....,'value') extracts data 
% from the BASIS A on the domain specified by the NAME/VALUES pairs. VALUES 
% specifies the A.BasisFunction.Parameter."NAME".GridData to be extracted.
%
% B = LPVSPLIT(A,NAME,RANGE) is an alternative syntax. NAME is an N-by-1
% cell array of characters and RANGE is an N-by-1 cell array of ranges.
% B = LPVSPLIT(A,NAME,INDEX,'index') and B = LPVSPLIT(A,NAME,VALUES,'value')
% also apply for INDEX or VALUES as a cell array. 
%
% B = LPVSPLIT(A,DOMAIN) is an another alternative syntax. DOMAIN is
% an RGRID object. This extracts data from the BASIS A based on the
% independent variables and data ranges in DOMAIN.
%
% See also: lpvinterp.

% XXX verify error catching here.
try
    BF = lpvsplit(m.BasisFunction,varargin{:});
    P = lpvsplit(m.Partials,varargin{:});
    IVN = m.IVName;
catch 
    rethrow(lasterr);
end
out = basis(BF,P,IVN);



