function B = adscalarexp(AData,Domain)
% This function takes array data and an rgrid domain and does scalar
% expansion, if necessary, on the array data to match up array dimensions
% in Domain.


szA = size(AData);
szA = szA(3:end);
% TODO: Size of A is 10x1, but size of the Domain is 1x10. This is
% NOT CORRECT as the size of Domain should be 10....
% Should this be handled in the calling function or here?
szD = size(Domain);
nA = numel(szA);
nD = numel(szD);
szA = [szA ones(1,nD-nA)];
szD = [szD ones(1,nA-nD)];

idx = find( szA==1 & szD>1 );

repsz = ones(1,length(szD));
repsz(idx) = szD(idx);

% repmat works on double/ss/frd
B = repmat( AData, [1 1 repsz]); 

