%% IQC test on a slowly TV real parameter 
% (page 102 of IQCbeta_2004 manual)
IQCbetaResult = 25.7744;

RB = 0.1;
IQClist = rbtvgain(RB,[-1 -3 -5]);
delta = udyn('delta',[1 1],'UserData',IQCcell(IQClist));
G = ss(tf(0.8,[1 0.21 1]));
H = feedback(G,delta);
wgrid = logspace(-3,3,90);
iqcstab4(H,wgrid)
tic
gainbnd = iqcperf4(H,wgrid,'gain')
toc

%% Using new code
% This is where we're having problems.  We get an answer (11-12), but it
% seems largely unrelated to RateBound.  Cuold be OK, but since we don't
% know answer, we don't really know.  MR/IQCtools bound is about 25.  Bound
% with UREAL is about 8.74 (WCGAIN gets 8.3 - difference is due entirely to
% freq grid - i checked this)
RB = 0.1;
S = IQCtvreal([-1 -3 -5],0.1);
delta = udyn('delta',[1 1],'UserData',S);
G = ss(tf(0.8,[1 0.21 1]));
H = feedback(delta*G,1);
wgrid = logspace(-4,4,200);
gainbnd = iqcperf4(H,wgrid,'gain')

G = ss(tf(0.8,[1 0.21 1]));
H = feedback(delta*G,1);
wgrid = logspace(-4,4,200);
tic
gainbnd = iqcperf4(H,wgrid,'gain')
toc

%% Two IQC functions
RB = 0.1;
S(1) = IQCcell(rbtvgain(RB,[-1 -3 -5]));
S(2) = IQCtvreal([-1 -3 -5],0.1);
S(3) = S(2);

delta = udyn('delta',[1 1],'UserData',S);
wgrid = logspace(-3,3,90);
H = feedback(delta*G,1);
tic
gainbnd = iqcperf4(H,wgrid,'gain')
toc

%% pg 102, no rate bounds, from rate-bound example
S = IQCtvrealnrb;
delta = udyn('delta',[1 1],'UserData',S);
%delta = ureal('delta',0);
G = ss(tf(0.2,[1 0.21 1]));
H = feedback(G,delta);
wgrid = logspace(-3,3,90);
%iqcstab2(H,wgrid)
gainbnd = iqcperf3(H,wgrid,'gain')

%% TV real, no rate bound (IQC-Beta iqc tvscalar, Sec 4.31)
a=-0.3;
b=0.8;
k=2.5;
Td=0.86;
A=[-(a+b*k*Td), -b*k; 1 0];
B=[-2, -2*b; 0, 0];
C=[a+b*k*Td, b*k; k*Td, k];
D=[1, 2*b;0, 1];
G=ss(A,B,C,D);
S = IQCtvrealnrb;
delta = udyn('delta',[1 1],'UserData',S);
Delta = [delta 0;0 delta];
deltaMag = 0.225;
H = feedback(deltaMag*G,Delta);
HH = feedback(Delta*deltaMag*G,eye(2));

wgrid = [0 logspace(-4,4,1600)];
gainbnd = iqcperf4(H,wgrid,'gain')
gainbnd = iqcperf4(HH,wgrid,'gain')




%% TVREAL alone
% freq grid - i checked this)
RB = 1000;
clear UD
UD.IQCfunction = @IQCtvreal;
UD.IQCparams.polelist = [-0.1 -1 -5 -20];
UD.IQCparams.RB = RB;
G = tf(1,[1 1]);
delta = udyn('delta',[1 1],'UserData',UD);
H = G*delta - delta*G;
wgrid = logspace(-4,4,200);
%iqcstab2(H,wgrid)
gainbnd = iqcperf3(H,wgrid,'gain')

delta = idf('delta',linspace(-1,1,40),[-RB RB]);
Xb(1).BasisFunction = pmat(1);
Xb(2).BasisFunction = delta;
Xb(3).BasisFunction = delta^2;
Xb(4).BasisFunction = cos(delta);
Xb(5).BasisFunction = sin(delta);
Xb(1).Partial = {pmat(0) , 'delta'};
Xb(2).Partial = {pmat(1) , 'delta'};
Xb(3).Partial = {2*delta, 'delta'};
Xb(4).Partial = {-sin(delta), 'delta'};
Xb(5).Partial = {cos(delta), 'delta'};

G = tf(1,[1 1]);
H = G*delta - delta*G;
H = G*delta;
xp = lpvnorm(H,Xb)


%% Zames-Falb example (1 uncertainty, perf works, all different ways)
G1 = tf([2 1 2],[1 20 100]);
G2 = tf([1 100],[1 5 20]);
G = ss(G1*G2);
f = 10;
SLP = 100;
IQClist = zamesfalb([-1],[0 f*SLP]);
IQClist = [IQClist, popov];
delta = udyn('delta',[1 1],'UserData',IQClist);
H = feedback(G/f,delta)*f;
wgrid = [0 logspace(-4,3,300)];
iqcstab2(H,wgrid)
[gainbnd,Sopt] = iqcperf3(H,wgrid,'gain');
gainbnd

%% IQC-B
s=tf([1 0],1);
G=((2*s*s+s+2)*(s+100))/((s+10)*(s+10)*(s*s+5*s+20));
abst_init_iqc;
w=signal;
f=signal;
v=G*(-w+f);
w==iqc_slope_odd(v,-1,0,0.5,100);
w==iqc_popov_vect(v);
gain=iqc_gain_tbx(f,v)


%%
UD.IQCfunction = @IQCsectorbounded;
UD.IQCparams.SB = [0 5];
delta = udyn('delta',[1 1],'UserData',UD);
H = feedback(G,delta);
wgrid = logspace(-3,3,160);
iqcstab2(H,wgrid)
gainbnd = iqcperf3(H,wgrid,'gain')

% direct creation for sector bound [0 5]
delta = udyn('delta',[1 1],'UserData',{[0 5;5 -2]});
H = feedback(G,delta);
wgrid = logspace(-3,3,160);
iqcstab2(H,wgrid)
gainbnd = iqcperf3(H,wgrid,'gain')

%% UREAL
slope = 1;
delta = ureal('delta',slope/2,'Range',[.1 slope]);
%delta = ureal('delta',0);
H = feedback(G,delta);
wgrid = logspace(-3,3,160);
%iqcstab2(H,wgrid)
gainbnd = iqcperf3(H,wgrid,'gain')
%gainbnd = iqcperf(H,wgrid,'gain')
x = wcgain(ufrd(H,wgrid))


%% LPV 
RB = 0.1;
%delta = tvreal('delta',0,'Range',[-1 1],'RateBounds',[-RB RB]);
delta = idf('delta',linspace(-1,1,40),[-RB RB]);
Xb(1).BasisFunction = pmat(1);
Xb(2).BasisFunction = delta;
Xb(3).BasisFunction = delta^2;
Xb(4).BasisFunction = cos(delta);
Xb(5).BasisFunction = sin(delta);
Xb(1).Partial = {pmat(0) , 'delta'};
Xb(2).Partial = {pmat(1) , 'delta'};
Xb(3).Partial = {2*delta, 'delta'};
Xb(4).Partial = {-sin(delta), 'delta'};
Xb(5).Partial = {cos(delta), 'delta'};

G = ss(tf(0.8,[1 0.21 1]));
H = feedback(G,delta);
xp = lpvnorm(H,Xb)


