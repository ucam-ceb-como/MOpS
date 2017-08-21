clc, clear, close all

set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

pdiags = csvread('Part-split-diagnostics(stage1).csv',1);
cdiags = csvread('Chem-split-diagnostics(stage1).csv',1);

figure(1)
set(gcf,'color','white')
subplot(121)
plot(pdiags(:,1)*1000,pdiags(:,4))
hold on
plot(pdiags(:,1)*1000,pdiags(:,5),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('SV (-)')
subplot(122)
plot(pdiags(:,1)*1000,pdiags(:,6))
hold on
plot(pdiags(:,1)*1000,pdiags(:,7),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('NSP (-)')

figure(2)
set(gcf,'color','white')
subplot(131)
plot(pdiags(:,1)*1000,sum(pdiags(:,8:8+105),2))
xlabel('Time (ms)')
ylabel('Inc. events (-)')
subplot(132)
plot(pdiags(:,1)*1000,pdiags(:,8+106))
xlabel('Time (ms)')
ylabel('Surf. growth events (-)')
subplot(133)
plot(pdiags(:,1)*1000,sum(pdiags(:,8+107:end-2),2))
xlabel('Time (ms)')
ylabel('Coag. events (-)')

% figure(3)
% set(gcf,'color','white')
% subplot(131)
% loglog(pdiags(:,1)*1000,cumsum(sum(pdiags(:,8:8+105),2)))
% xlabel('Time (ms)')
% ylabel('Inc. events (-)')
% subplot(132)
% loglog(pdiags(:,1)*1000,cumsum(pdiags(:,8+106)))
% xlabel('Time (ms)')
% ylabel('Surf. growth events (-)')
% subplot(133)
% loglog(pdiags(:,1)*1000,cumsum(sum(pdiags(:,8+107:end-2),2)))
% xlabel('Time (ms)')
% ylabel('Coag. events (-)')

totevs = sum(pdiags(:,8:end-2),2);
figure(4)
set(gcf,'color','white')
loglog(pdiags(:,1)*1000,sum(pdiags(:,8:8+105)./totevs,2))
hold on
loglog(pdiags(:,1)*1000,pdiags(:,8+106)./totevs,'--')
loglog(pdiags(:,1)*1000,sum(pdiags(:,8+107:end-2),2)./totevs,':')
xlabel('Time (ms)')
ylabel('Fraction of events (-)')
legend('Inc.','Surf. growth','Coag.','location','southwest')

figure(5)
set(gcf,'color','white')
subplot(231)
loglog(cdiags(:,1)*1000,cdiags(:,4)*1e6)
hold on
loglog(cdiags(:,1)*1000,cdiags(:,5)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('TiCl$_4$ conc. (mol$\cdot$m$^{-3}$)')
subplot(232)
loglog(cdiags(:,1)*1000,cdiags(:,40)*1e6)
hold on
loglog(cdiags(:,1)*1000,cdiags(:,41)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('O$_2$ conc. (mol$\cdot$m$^{-3}$)')
subplot(233)
loglog(cdiags(:,1)*1000,cdiags(:,58)*1e6)
hold on
loglog(cdiags(:,1)*1000,cdiags(:,59)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Ar conc. (mol$\cdot$m$^{-3}$)')
subplot(234)
loglog(cdiags(:,1)*1000,cdiags(:,42)*1e6)
hold on
loglog(cdiags(:,1)*1000,cdiags(:,43)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Cl$_2$ conc. (mol$\cdot$m$^{-3}$)')
subplot(235)
loglog(cdiags(:,1)*1000,cdiags(:,20)*1e6)
hold on
loglog(cdiags(:,1)*1000,cdiags(:,21)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('TiO$_2$Cl$_3$ conc. (mol$\cdot$m$^{-3}$)')
subplot(236)
loglog(cdiags(:,1)*1000,cdiags(:,44)*1e6)
hold on
loglog(cdiags(:,1)*1000,cdiags(:,45)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('O conc. (mol$\cdot$m$^{-3}$)')
% 
% cdiags0 = csvread('Chem-split-diagnostics(stage0).csv',1);
% 
% figure(6)
% set(gcf,'color','white')
% subplot(121)
% loglog(cdiags0(:,1)*1000,cdiags0(:,40)*1e6)
% hold on
% loglog(cdiags0(:,1)*1000,cdiags0(:,41)*1e6,':')
% legend('Pre','Post','location','North','orientation','horizontal')
% xlabel('Time (ms)')
% ylabel('O$_2$ conc. (mol$\cdot$m$^{-3}$)')
% subplot(122)
% loglog(cdiags0(:,1)*1000,cdiags0(:,44)*1e6)
% hold on
% loglog(cdiags0(:,1)*1000,cdiags0(:,45)*1e6,':')
% legend('Pre','Post','location','North','orientation','horizontal')
% xlabel('Time (ms)')
% ylabel('O conc. (mol$\cdot$m$^{-3}$)')
