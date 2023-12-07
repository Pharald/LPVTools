function out = ss(a,b,c,d,varargin)
% Call PSS constructor
%
% SYS = ss(A,B,C,D) creates a parameter-varying state-space model on  
%     the domain of A, B, C and D. SYS will be a PSS.
%
% SYS = ss(A)  creates a parameter dependent static gain matrix on the 
%     domain of A.
% 
% See also: ss, pss.
 

if nargin==1
    out = pss(a);
elseif nargin==4
    out = pss(a,b,c,d);
else
    out = pss(a,b,c,d,varargin{:});
end




% ----- OLD CODE
%
%
% if nargin==0
%    disp('usage: out = ss(a,b,c,d)');
%    return
% end
%
% % if nargin==4
% %    Ts = 0;
% % end
%
% sza = size(a);
% szb = size(b);
% szc = size(c);
% szd = size(d);
%
% if sza(1)~=sza(2)
%    error('A-matrix needs to be square');
% end
% if sza(1)~=szb(1)
%    error('A and B matrices need the same number of rows');
% end
% if sza(2)~=szc(2)
%    error('A and C matrices need the same number of columns');
% end
% if szd(1)~=szc(1)
%    error('C and D matrices need the same number of rows');
% end
% if szd(2)~=szb(2)
%    error('B and D matrices need the same number of columns');
% end
%
% sn = cell(sza(1),1);
% for i=1:sza(1)
%    sn{i} = '';
% end
% ltifield = lti(szc(1),szb(2),Ts);
% mat = [a b;c d];
% switch class(mat)
% case 'double'
%    Nx = sza(1);
%    out = ss(mat(1:Nx,1:Nx),mat(1:Nx,Nx+1:end),...
%          mat(Nx+1:end,1:Nx),mat(Nx+1:end,Nx+1:end),ltifield);
% case 'umat'
%    out = umdl(mat,sn,ltifield);
% case 'pmat'
%    out = pmdl(mat,sn,ltifield);
% case 'pumat'
%    out = pumdl(mat,sn,ltifield);
% otherwise
%    error('Incorrect state-space data')
% end

