function s = isstable(sys)
% ISSTABLE  Pointwise stability test for a PSS object.
%
% ISSTABLE(SYS) tests the stability of SYS at each point in the domain 
% of SYS. Returns a PMAT on the domain of SYS, containing True at those 
% points where SYS is stable, and False otherwise.
%  
% See also: isstable, pole.

try
    szs = privatesize(sys);
    s = isstable(sys.DataPrivate,'-v7.8');
    s = pmat( reshape(s,[1 1 szs(3:end)]) , sys.DomainPrivate);
catch ME
    throw(ME);
end
    
