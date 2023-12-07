%% LPVSPLIT - Extract data based on independent variable range
%
%  
%% Syntax
%
%    B = lpvsplit(A,NAME1,RANGE1,NAME2,RANGE2,....)
%    B = lpvsplit(A,NAME1,INDEX1,NAME2,INDEX2,....,'index')
%    B = lpvsplit(A,NAME1,VALUES1,NAME2,VALUES2,....,'value')
%    B = lpvsplit(A,NAME,RANGE)
%    B = lpvsplit(A,NAME,INDEX,'index')
%    B = lpvsplit(A,NAME,VALUES,'value')
%    B = lpvsplit(A,DOMAIN)
%
%% Description
% 
% 
% |B = lpvsplit(A,NAME1,RANGE1,NAME2,RANGE2,....)| extracts data from the
% grid-based LPV object |A| on the domain specified by the |NAME| / |RANGE| pairs. Each |RANGE| is a 1-by-2 
% row vector [min, max], that specifies the values of the independent 
% variable to be extracted along the domain direction |NAME|. The data at all
% points in |A.Parameter."NAME".GridData| which lies inside |RANGE| is extracted. 
% If |RANGE| is a scalar then the data is extracted where the variable |NAME| 
% is exactly equal to |RANGE|. If an independent variable of |A| is not listed 
% in the inputs then all values along this domain direction are retained in |B|.
%
% |B = lpvsplit(A,NAME1,INDEX1,NAME2,INDEX2,....,'index')| extracts data from
% the |pmat| |A| on the domain specified by the |NAME| / |INDEX| pairs.
% Each |INDEX| is a vector of integers or a logical array that specifies 
% the indices of |A.Parameter."NAME".GridData| to be extracted.
%
% |B = lpvsplit(A,NAME1,VALUES1,NAME2,VALUES2,....,'value')| extracts data 
% from the |pmat| |A| on the domain specified by the |NAME| / |VALUES| pairs. 
% |VALUES| specifies the |A.Parameter."NAME".GridData| to be extracted.
%
% |B = lpvsplit(A,NAME,RANGE)| is an alternative syntax. |NAME| is an N-by-1
% cell array of characters and |RANGE| is an N-by-1 cell array of ranges.
% |B = lpvsplit(A,NAME,INDEX,'index')| and |B = lpvsplit(A,NAME,VALUES,'value')|
% also apply for |INDEX| or |VALUES| as a cell array. 
%
% |B = lpvsplit(A,DOMAIN)| is an another alternative syntax. |DOMAIN| is
% an |rgrid| object. This extracts data from the |pmat| |A| based on the
% independent variables and data ranges in |DOMAIN|.