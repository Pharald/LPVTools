t = 0:0.002:1000;
N = numel(t);
d = (N-1)/4;
scale = 1000;
hleg1= zeros(1,d);
hleg2 = (t(d+1:2*d)-t(d))/120+hleg1(end);
hleg3 = cos(t(1:d)/250)+hleg2(end)-1;
hleg4 = zeros(1,d+1)+hleg3(end);
h = [hleg1,hleg2,hleg3,hleg4]*scale+6000;

% Plot settings
fz =10;
     
f1=figure(1);
clf
subplot(211)
plot(t,h,'b')     
hold on
plot([t(1) t(end)],[10000 10000],'r--')
plot([t(1) t(end)],[5000 5000],'r--')
xlabel('time [sec]')
ylabel('Altitude [ft]')
set(gca,'FontSize',fz)
set(findall(gcf,'type','text'),'FontSize',fz,'fontWeight','bold')
set(findall(gcf,'type','Line'),'LineWidth',2) 
ylim([3000 12000])
legend('Altitude','Limits on altitude','location','best')

hdotleg1 = zeros(1,d);
hdotleg2 = ones(1,d)*scale/120;
hdotleg3 = -sin(t(1:d)/250)*scale/250;
hdotleg4 = zeros(1,d+1);

hdot = [hdotleg1,hdotleg2,hdotleg3,hdotleg4];

subplot(212)
plot(t,hdot,'k:') 
hold on 
plot([t(1) t(end)],[10 10],'r--')
plot([t(1) t(end)],[-10 -10],'r--')
xlabel('time [sec]')
ylabel('Rate of Change in Altitude [ft/s]')
set(gca,'FontSize',fz)
set(findall(gcf,'type','text'),'FontSize',fz,'fontWeight','bold')
set(findall(gcf,'type','Line'),'LineWidth',2) 
ylim([-15 15])
legend('Rate of climb/descent','Limits on climb/descent rate','location','best')


pos = [2,1,16,12];
set(gcf,'Units','centimeters');
set(gcf,'Position',pos);
% print(f1,'html\AltitudeExample','-dpng')

set(gcf,'PaperPositionMode','auto')
print(f1,'html\AltitudeExample','-dpng','-r0')




