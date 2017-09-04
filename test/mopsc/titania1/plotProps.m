% clc, clear, close all

set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

fold = '';
legd = '';

part = csvread([fold 'Network(stage1)-part.csv'],1);
rate = csvread([fold 'Network(stage1)-part-rates.csv'],1);
pslf = csvread([fold 'Network(stage1)-psl(0.1s).csv'],1);
chem = csvread([fold 'Network(stage1)-chem.csv'],1);
cput = csvread([fold 'Network(stage1)-cput.csv'],1);

%% Part

figure(1)
set(gcf,'color','white')
subplot(231)
plot(part(:,2)*1000,part(:,3))
hold on
xlabel('Time (ms)')
ylabel('NSP (-)')
legend(legd)

subplot(232)
% figure(2)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,5))
hold on
xlabel('Time (ms)')
ylabel('M0 (m$^{-3}$)')
legend(legd)

subplot(233)
% figure(3)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,9)*1e9)
hold on
xlabel('Time (ms)')
ylabel('Collision diamter (nm)')
legend(legd)

subplot(234)
% figure(4)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,21))
hold on
xlabel('Time (ms)')
ylabel('Mass (kg$\cdot$m$^{-3}$)')
legend(legd)

subplot(235)
% figure(5)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,37))
hold on
xlabel('Time (ms)')
ylabel('Number primaries per particle (-)')
legend(legd)

subplot(236)
% figure(6)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,39)*1e9)
hold on
xlabel('Time (ms)')
ylabel('Primary average diamter (nm)')
legend(legd)

%% PSL

nbins=25;

figure(2)
% figure()
set(gcf,'color','white')
[n_d,x_d] = histwc(pslf(:,3),pslf(:,9),nbins);
subplot(231)
plot(x_d,n_d/max(n_d))
hold on
xlabel('Collision diameter (nm)')
ylabel('Count divided by max. count (-)')
legend(legd)

[n_d,x_d] = histwc(pslf(:,5),pslf(:,9),nbins);
subplot(232)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Surface area (cm$^{2}$)')
ylabel('Count divided by max. count (-)')
legend(legd)

[n_d,x_d] = histwc(pslf(:,12),pslf(:,9),nbins);
subplot(233)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Number of primaries (nm)')
ylabel('Count divided by max. count (-)')
legend(legd)

[n_d,x_d] = histwc(pslf(:,13),pslf(:,9),nbins);
subplot(234)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Primary diameter (nm)')
ylabel('Count divided by max. count (-)')
legend(legd)

[n_d,x_d] = histwc(pslf(:,16),pslf(:,9),nbins);
subplot(235)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Std of primary diameter (nm)')
ylabel('Count divided by max. count (-)')
legend(legd)

[n_d,x_d] = histwc(pslf(:,8),pslf(:,9),nbins);
subplot(236)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Particle age (s)')
ylabel('Count divided by max. count (-)')
legend(legd)

figure(20)
set(gcf,'color','white')
[n_d,x_d] = hist(pslf(:,9));
% subplot(231)
plot(x_d,n_d/max(n_d))
hold on
xlabel('Weight (-)')
ylabel('Count divided by max. count (-)')
legend(legd)

%% Chem

figure(3)
set(gcf,'color','white')
subplot(231)
plot(chem(:,2)*1000,chem(:,3)*1e6)
hold on
xlabel('Time (ms)')
ylabel('TiCl$_4$ conc. (mol$\cdot$m$^{-3}$)')
legend(legd)

subplot(232)
% figure(2)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,39)*1e6)
hold on
xlabel('Time (ms)')
ylabel('O$_2$ conc. (mol$\cdot$m$^{-3}$)')
legend(legd)

subplot(233)
% figure(3)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,45)*1e6)
hold on
xlabel('Time (ms)')
ylabel('Cl$_2$ conc. (mol$\cdot$m$^{-3}$)')
legend(legd)

subplot(234)
% figure(4)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,19)*1e6)
hold on
xlabel('Time (ms)')
ylabel('TiO$_2$Cl$_3$ conc. (mol$\cdot$m$^{-3}$)')
legend(legd)

subplot(235)
% figure(5)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,57)*1e6)
hold on
xlabel('Time (ms)')
ylabel('Ar conc. (mol$\cdot$m$^{-3}$)')
legend(legd)

subplot(236)
% figure(6)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,61))
hold on
xlabel('Time (ms)')
ylabel('Temperature (K)')
legend(legd)

%% CPUT

figure(4)
set(gcf,'color','white')
plot(cput(:,2)*1000,cput(:,3)/60)
hold on
xlabel('Time (ms)')
ylabel('Solver time (min)')
legend(legd)

%% Rates

figure(5)
set(gcf,'color','white')
subplot(221)
plot(rate(:,2)*1000,sum(rate(:,3:2:211),2))
hold on
xlabel('Time (ms)')
ylabel('Inc. rate (m$^{-3}\cdot$s$^{-1}$)')
legend(legd)

subplot(222)
% figure(3)
% set(gcf,'color','white')
plot(rate(:,2)*1000,rate(:,215))
hold on
xlabel('Time (ms)')
ylabel('Surf. growth rate (m$^{-3}\cdot$s$^{-1}$)')
legend(legd)

subplot(223)
% figure(3)
% set(gcf,'color','white')
plot(rate(:,2)*1000,sum(rate(:,217:2:end),2))
hold on
xlabel('Time (ms)')
ylabel('Coag. rate (m$^{-3}\cdot$s$^{-1}$)')
legend(legd)

subplot(224)
% figure(2)
% set(gcf,'color','white')
plot(rate(:,2)*1000,rate(:,213))
hold on
xlabel('Time (ms)')
ylabel('No inc. rate (m$^{-3}\cdot$s$^{-1}$)')
legend(legd) 