% test script for functionality of connect()

addpath(genpath(pwd))

%%
% generate pss object
rho = rgrid('rho',[-1 0 1]);

for ii = 1:3
Pdata(:,:,ii) = tf(1,[1 ii]);

% Pdata(:,:,ii) = ss(1,ii,2,0);
end

P = pss(Pdata,rho);
P.outputname = 'y';
P.inputname = 'u';

K = tf(1,[1 2 3]); % lti system
K.inputname = 'e';
K.outputname = 'u';

for ii = 1:3
Kdata(:,:,ii) = tf(1,[1 ii 3]);
end

Kpss = pss(Kdata,rho);
Kpss.inputname = 'e';
Kpss.outputname = 'u';

% testing functions
Pfeedback = feedback(P*K,1);

%%

clc

try
sum1 = sumblk('u = r-y');
Pconnect = connect(sum1,P,'r','y');
catch 
sprintf('connect not working with pss connected to sumblk')
end

%%
try
Psys = append(K,P);
connections = [2 1; 1 -2];
connect(Psys,connections,1,2);
catch
sprintf('connect is not working with pss and in/out defined by matrix of connections')
end
%%

try
Psys = append(Kpss,P);
connections = [2 1; 1 -2];
connect(P,Kpss,connections,1,2);
catch
sprintf('connect is not working for combining pss objects defined by matrix of connections')
end

try
Pconnect = connect(sum1,P,Kpss,'r','y');
catch
sprintf('connect is not working for combining pss objects')
end


