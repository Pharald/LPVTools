function s = lpvsample(r,Npts,opt)
% LPVSAMPLE   Sample a rectangular grid object RGRID
%
% S=LPVSAMPLE(R,N) returns a sample of N points from the rectangular 
% grid R. S is an Niv-by-N matrix with each column containing one sample of the
% rectangular grid.  The rows of S correspond to the ordering of 
% independent variables in R.IVName.  
%
% S=LPVSAMPLE(R,N,OPT) allows the user to specify the sampling algorithm
% to be used. OPT is a CHAR specifies the type of sampling:
%    -'grid': Draws points drawn randomly (possibly with repeats) from 
%         the rectangular grid.
%    -'uniform' (default): Draws points uniformly from the hypercube
%         specified by the limits of the rectangular grid.
%    -'LHC': Does a Latin Hypercube sample of the domain hypercube.
% For 'uniform' and 'LHC', the samples are not, in general, elements
% of the rectangular grid.
%
%   % EXAMPLE: (CUT/PASTE)
%   % Sample a 2-dimensional RGRID object
%   R = rgrid( {'a', 'b'}, {linspace(-2,2,12), 1:5} );
%   Su = lpvsample(R,15);   % Uniform sample
%   Sg = lpvsample(R,15,'grid');   % Sample from grid
%   plot(Su(1,:),Su(2,:),'bx',Sg(1,:),Sg(2,:),'ro');
%   legend('Uniform','Grid','Location','Best')

% Check # of input/output arguments
nin = nargin;
narginchk(2, 3)
if nin==2
    opt = 'uniform';
end

r = rmprivateiv(r); % Remove Any Extra Array Dimensions Names 

% Get LB/UB Range data
Niv = r.NumIV;
IVData = r.IVData;
LIVData = r.LIVData;
Range = zeros(Niv,2);
for i=1:Niv
    Range(i,:)=[IVData{i}(1) IVData{i}(end)];
end

switch opt
    case 'grid'
        s = zeros(Niv,Npts);
        for i=1:Niv
            idx = ceil( LIVData(i)*rand(1,Npts) );
            s(i,:) = IVData{i}(idx);
        end        
    case 'uniform'
        LB = Range(:,1);
        LEN = Range(:,2)-Range(:,1);        
        s = repmat(LB,[1 Npts]) + repmat(LEN,[1 Npts]).*rand(Niv,Npts);
    case 'LHC'
        s = LOCALlhcdesign(Niv,Npts,Range);
    otherwise
        error(['Sample type ' opt ' not implemented.']);
end
    



function [D,DN] = LOCALlhcdesign(n,p,Range)
% Latin Hypercube Design
% [0 1]^n, p points total
%
% % Example
% % Inputs:
% n = 4;  % 4-dimensional space
% p = 50; % 50 samples
% Range = [-1 1;0 5;-2 0;4 5];
% D = LOCALlhcdesign(n,p,Range);
% D  % 4-by-50 sample matrix in Range

V = rand(p,n);
ILength = 1/p;
IBegin = (0:p-1)'*ILength;

V = repmat(IBegin,[1 n]) + ILength*V;
for i=1:n
   [~,idx] = sort(rand(1,p));
   V(:,i) = V(idx,i);
end
DN = V';
% Range n-by-2
Left = Range(:,1);
Width = Range(:,2) - Range(:,1);
D = repmat(Left,[1 p])+DN.*repmat(Width,[1 p]);

