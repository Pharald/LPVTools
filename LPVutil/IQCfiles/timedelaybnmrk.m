% Time-Delay System
% Weehon Tan, Berkeley
% Automatica vol33, 1997, pp411-419

% plant model
%clear all

taud = 0.1;

P = tf(10,[1 1]);
C = ss(1);

D = ss(1);
D.InputDelay = taud;

omega= logspace(-2,2,1000);

%Specify approximation
ordp = 0; ordw = 1;
% ordp = 1; ordw = 3;
%ordp = 2; ordw = 3;

[N,weight,Erf]=covertd(taud,ordp,ordw);

M = blkdiag(1,-1);
s = tf('s');
w1 = 1;
w2 = 1/(s+1);
w3 = 1/(s+.01);
w4 = 1/(s+.1);

psi1 = blkdiag(w1*weight,w1);
psi2 = blkdiag(w2*weight,w2);
psi3 = blkdiag(w3*weight,w3);
psi4 = blkdiag(w4*weight,w4);

IQClist = {psi1'*M*psi1,psi4'*M*psi4};

deltadel = udyn('deltadel',[1 1],'UserData',IQCcell(IQClist));

CL = feedback(1,P*(N+deltadel)*C);

[G,Delta] = lftdata(CL);
blksize = [1 1];

iqcs.psi = {psi1}; 
iqcs.M = {M};
iqcs.psi{2} = psi4; 
iqcs.M{2} = M;

Fbasis = 1;
Fgrad = zeros(1,0);
RateBounds = [];

tic
gam = iqcgridengine(G,Fbasis,Fgrad,RateBounds,blksize,iqcs);
%gam = iqcgridengine(G,blksize,iqcs);
toc

disp(['IQC LMILab iqcgain:       ', num2str(gam)]);

% Test an array G
G = cat(3,G,G,G,G);

Fbasis = ones(1,1,4);
Fgrad = zeros(1,0);
RateBounds = [];

tic
gam = iqcgridengine(G,Fbasis,Fgrad,RateBounds,blksize,iqcs);
%gam = iqcgridengine(G,blksize,iqcs);
toc

disp(['IQC LMILab iqcgain:       ', num2str(gam)]);


