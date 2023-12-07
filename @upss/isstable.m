function s = isstable(sys)
% ISSTABLE  Pointwise stability test for the nominal value of a UPSS.
%
% ISSTABLE(SYS) tests the stability of the nominal value of SYS at each 
% point in the domain of SYS. Returns a PMAT on the domain of SYS, 
% containing True at those  points where SYS is stable, and False otherwise.
%  
% See also: isstable, pole.

try
    szs = privatesize(sys);
    s = isstable(sys.DataPrivate,'elem');
    s = pmat( reshape(s,[1 1 szs(3:end)]) , sys.DomainPrivate);
catch ME
    throw(ME);
end
    
