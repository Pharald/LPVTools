function H = surf(a)

% TODO PJS 4/4/2011: Revisit.
% I implemented the 2d plot functions (plot, semilogx, etc) by generating
% the plots at each point in the domain on a single. This seems consistent
% with our basic approach to overloading DOUBLE/SS/FRD functions. 
% This SURF function deviates from this approach --> Given a 1-by-1
% PMAT A with 2IVs it generates a surf plot of the value as a function
% of the IVs.  Should we rename this function since it does not overload
% the DOUBLE/SURF in the same way that we are overloading other DOUBLE
% functions?

sza = privatesize(a);
IvData = a.DomainPrivate.IVData;
niv = a.DomainPrivate.NumIV;

if niv==2 && sza(1)==1 && sza(2)==1
    h = surf(IvData{2}, IvData{1}, reshape(a.DataPrivate,sza(3:end)));
    xlabel(a.DomainPrivate.IVName{2});
    ylabel(a.DomainPrivate.IVName{1});
else
    error('PMAT should have 2 independent variables')
end
if nargout==1
   H = h;
end
