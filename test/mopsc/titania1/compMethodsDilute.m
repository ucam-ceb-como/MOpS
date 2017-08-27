clc, clear, close all

set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

fold = 'vart-2048SPs/'; 

part_dd = csvread([fold 'dsa-d/PP-part.csv'],1);
part_sd = csvread([fold 'swa-d/PP-part.csv'],1);
part_wd1 = csvread([fold 'swa-dn1/PP-part.csv'],1);
part_wd2 = csvread([fold 'swa-dn2/PP-part.csv'],1);

rate_dd = csvread([fold 'dsa-d/PP-part-rates.csv'],1);
rate_sd = csvread([fold 'swa-d/PP-part-rates.csv'],1);
rate_wd1 = csvread([fold 'swa-dn1/PP-part-rates.csv'],1);
rate_wd2 = csvread([fold 'swa-dn2/PP-part-rates.csv'],1);

psl_dd = csvread([fold 'dsa-d/PP-psl(0.5s).csv'],1);
psl_sd = csvread([fold 'swa-d/PP-psl(0.5s).csv'],1);
psl_wd1 = csvread([fold 'swa-dn1/PP-psl(0.5s).csv'],1);
psl_wd2 = csvread([fold 'swa-dn2/PP-psl(0.5s).csv'],1);

chem_dd = csvread([fold 'dsa-d/PP-chem.csv'],1);
chem_sd = csvread([fold 'swa-d/PP-chem.csv'],1);
chem_wd1 = csvread([fold 'swa-dn1/PP-chem.csv'],1);
chem_wd2 = csvread([fold 'swa-dn2/PP-chem.csv'],1);

cput_dd = csvread([fold 'dsa-d/PP-cput.csv'],1);
cput_sd = csvread([fold 'swa-d/PP-cput.csv'],1);
cput_wd1 = csvread([fold 'swa-dn1/PP-cput.csv'],1);
cput_wd2 = csvread([fold 'swa-dn2/PP-cput.csv'],1);

part_ds = csvread([fold 'dsa-s/PP-part.csv'],1);
part_ss = csvread([fold 'swa-s/PP-part.csv'],1);

rate_ds = csvread([fold 'dsa-s/PP-part-rates.csv'],1);
rate_ss = csvread([fold 'swa-s/PP-part-rates.csv'],1);

psl_ds = csvread([fold 'dsa-s/PP-psl(0.5s).csv'],1);
psl_ss = csvread([fold 'swa-s/PP-psl(0.5s).csv'],1);

chem_ds = csvread([fold 'dsa-s/PP-chem.csv'],1);
chem_ss = csvread([fold 'swa-s/PP-chem.csv'],1);

cput_ds = csvread([fold 'dsa-s/PP-cput.csv'],1);
cput_ss = csvread([fold 'swa-s/PP-cput.csv'],1);

leg = {'DSA det, orig','SWA det, orig',...
       'DSA sph, orig','SWA sph, orig',...
       'SWA det, wt const','SWA det, wt var'};

legs = {'DSA det, orig','SWA det, orig',...
       'SWA det, wt const','SWA det, wt var'};

%% Part

figure(1)
set(gcf,'color','white')
subplot(231)
plot(part_dd(:,2)*1000,part_dd(:,3))
hold on
plot(part_sd(:,2)*1000,part_sd(:,3),':')
plot(part_ds(:,2)*1000,part_ds(:,3),'--')
plot(part_ss(:,2)*1000,part_ss(:,3),'-.')
plot(part_wd1(:,2)*1000,part_wd1(:,3),'--')
plot(part_wd2(:,2)*1000,part_wd2(:,3),'-.')
xlabel('Time (ms)')
ylabel('NSP (-)')
legend(leg)

subplot(232)
% figure(2)
% set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,5))
hold on
plot(part_sd(:,2)*1000,part_sd(:,5),':')
plot(part_ds(:,2)*1000,part_ds(:,5),'--')
plot(part_ss(:,2)*1000,part_ss(:,5),'-.')
plot(part_wd1(:,2)*1000,part_wd1(:,5),'--')
plot(part_wd2(:,2)*1000,part_wd2(:,5),'-.')
xlabel('Time (ms)')
ylabel('M0 (m$^{-3}$)')
legend(leg)

subplot(233)
% figure(3)
% set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,9)*1e9)
hold on
plot(part_sd(:,2)*1000,part_sd(:,9)*1e9,':')
plot(part_ds(:,2)*1000,part_ds(:,9)*1e9,'--')
plot(part_ss(:,2)*1000,part_ss(:,9)*1e9,'-.')
plot(part_wd1(:,2)*1000,part_wd1(:,9)*1e9,'--')
plot(part_wd2(:,2)*1000,part_wd2(:,9)*1e9,'-.')
xlabel('Time (ms)')
ylabel('Collision diamter (nm)')
legend(leg)

subplot(234)
% figure(4)
% set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,21))
hold on
plot(part_sd(:,2)*1000,part_sd(:,21),':')
plot(part_ds(:,2)*1000,part_ds(:,21),'--')
plot(part_ss(:,2)*1000,part_ss(:,21),'-.')
plot(part_wd1(:,2)*1000,part_wd1(:,21),'--')
plot(part_wd2(:,2)*1000,part_wd2(:,21),'-.')
xlabel('Time (ms)')
ylabel('Mass (kg$\cdot$m$^{-3}$)')
legend(leg)

subplot(235)
% figure(5)
% set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,37))
hold on
plot(part_sd(:,2)*1000,part_sd(:,37),':')
plot(part_wd1(:,2)*1000,part_wd1(:,37),'--')
plot(part_wd2(:,2)*1000,part_wd2(:,37),'-.')
xlabel('Time (ms)')
ylabel('Number primaries per particle (-)')
legend(legs)

subplot(236)
% figure(6)
% set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,39)*1e9)
hold on
plot(part_sd(:,2)*1000,part_sd(:,39)*1e9,':')
plot(part_wd1(:,2)*1000,part_wd1(:,39)*1e9,'--')
plot(part_wd2(:,2)*1000,part_wd2(:,39)*1e9,'-.')
xlabel('Time (ms)')
ylabel('Primary average diamter (nm)')
legend(legs)

%% PSL

figure(2)
% figure()
set(gcf,'color','white')
[n_d,x_d] = hist(psl_dd(:,3));
[n_s,x_s] = hist(psl_sd(:,3));
[n_ds,x_ds] = hist(psl_ds(:,3));
[n_ss,x_ss] = hist(psl_ss(:,3));
[n_d2,x_d2] = hist(psl_wd1(:,3));
[n_s2,x_s2] = hist(psl_wd2(:,3));
subplot(231)
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_ds,n_ds/max(n_ds),'--',x_ss,n_ss/max(n_ss),'-.')
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Collision diameter (nm)')
ylabel('Count divided by max. count (-)')
legend(leg)

[n_d,x_d] = hist(psl_dd(:,5));
[n_s,x_s] = hist(psl_sd(:,5));
[n_d2,x_d2] = hist(psl_wd1(:,5));
[n_s2,x_s2] = hist(psl_wd2(:,5));
subplot(232)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Surface area (cm$^{2}$)')
ylabel('Count divided by max. count (-)')
legend(legs)

[n_d,x_d] = hist(psl_dd(:,12));
[n_s,x_s] = hist(psl_sd(:,12));
[n_d2,x_d2] = hist(psl_wd1(:,12));
[n_s2,x_s2] = hist(psl_wd2(:,12));
subplot(233)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Number of primaries (nm)')
ylabel('Count divided by max. count (-)')
legend(legs)

[n_d,x_d] = hist(psl_dd(:,13));
[n_s,x_s] = hist(psl_sd(:,13));
[n_d2,x_d2] = hist(psl_wd1(:,13));
[n_s2,x_s2] = hist(psl_wd2(:,13));
subplot(234)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Primary diameter (nm)')
ylabel('Count divided by max. count (-)')
legend(legs)

[n_d,x_d] = hist(psl_dd(:,16));
[n_s,x_s] = hist(psl_sd(:,16));
[n_d2,x_d2] = hist(psl_wd1(:,16));
[n_s2,x_s2] = hist(psl_wd2(:,16));
subplot(235)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Std of primary diameter (nm)')
ylabel('Count divided by max. count (-)')
legend(legs)

[n_d,x_d] = hist(psl_dd(:,8));
[n_s,x_s] = hist(psl_sd(:,8));
[n_d2,x_d2] = hist(psl_wd1(:,8));
[n_s2,x_s2] = hist(psl_wd2(:,8));
subplot(236)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Particle age (s)')
ylabel('Count divided by max. count (-)')
legend(legs)

%% Chem

figure(3)
set(gcf,'color','white')
subplot(231)
plot(chem_dd(:,2)*1000,chem_dd(:,3)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,3)*1e6,':')
plot(chem_ds(:,2)*1000,chem_ds(:,3)*1e6,'--')
plot(chem_ss(:,2)*1000,chem_ss(:,3)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,3)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,3)*1e6,'-.')
xlabel('Time (ms)')
ylabel('TiCl$_4$ conc. (mol$\cdot$m$^{-3}$)')
legend(leg)

subplot(232)
% figure(2)
% set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,39)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,39)*1e6,':')
plot(chem_ds(:,2)*1000,chem_ds(:,39)*1e6,'--')
plot(chem_ss(:,2)*1000,chem_ss(:,39)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,39)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,39)*1e6,'-.')
xlabel('Time (ms)')
ylabel('O$_2$ conc. (mol$\cdot$m$^{-3}$)')
legend(leg)

subplot(233)
% figure(3)
% set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,45)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,45)*1e6,':')
plot(chem_ds(:,2)*1000,chem_ds(:,45)*1e6,'--')
plot(chem_ss(:,2)*1000,chem_ss(:,45)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,45)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,45)*1e6,'-.')
xlabel('Time (ms)')
ylabel('Cl$_2$ conc. (mol$\cdot$m$^{-3}$)')
legend(leg)

subplot(234)
% figure(4)
% set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,19)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,19)*1e6,':')
plot(chem_ds(:,2)*1000,chem_ds(:,19)*1e6,'--')
plot(chem_ss(:,2)*1000,chem_ss(:,19)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,19)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,19)*1e6,'-.')
xlabel('Time (ms)')
ylabel('TiO$_2$Cl$_3$ conc. (mol$\cdot$m$^{-3}$)')
legend(leg)

subplot(235)
% figure(5)
% set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,57)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,57)*1e6,':')
plot(chem_ds(:,2)*1000,chem_ds(:,57)*1e6,'--')
plot(chem_ss(:,2)*1000,chem_ss(:,57)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,57)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,57)*1e6,'-.')
xlabel('Time (ms)')
ylabel('Ar conc. (mol$\cdot$m$^{-3}$)')
legend(leg)

subplot(236)
% figure(6)
% set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,61))
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,61),':')
plot(chem_ds(:,2)*1000,chem_ds(:,61),'--')
plot(chem_ss(:,2)*1000,chem_ss(:,61),'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,61),'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,61),'-.')
xlabel('Time (ms)')
ylabel('Temperature (K)')
legend(leg)

%% CPUT

figure(4)
set(gcf,'color','white')
plot(cput_dd(:,2)*1000,cput_dd(:,3)/60)
hold on
plot(cput_sd(:,2)*1000,cput_sd(:,3)/60,':')
plot(cput_ds(:,2)*1000,cput_ds(:,3)/60,'--')
plot(cput_ss(:,2)*1000,cput_ss(:,3)/60,'-.')
plot(cput_wd1(:,2)*1000,cput_wd1(:,3)/60,'--')
plot(cput_wd2(:,2)*1000,cput_wd2(:,3)/60,'-.')
xlabel('Time (ms)')
ylabel('Solver time (min)')
legend(leg)

%% Rates

figure(5)
set(gcf,'color','white')
subplot(221)
plot(rate_dd(:,2)*1000,sum(rate_dd(:,3:2:211),2))
hold on
plot(rate_sd(:,2)*1000,sum(rate_sd(:,3:2:211),2),':')
plot(rate_ds(:,2)*1000,sum(rate_ds(:,3:2:211),2),'--')
plot(rate_ss(:,2)*1000,sum(rate_ss(:,3:2:211),2),'-.')
plot(rate_wd1(:,2)*1000,sum(rate_wd1(:,3:2:211),2),'--')
plot(rate_wd2(:,2)*1000,sum(rate_wd2(:,3:2:211),2),'-.')
xlabel('Time (ms)')
ylabel('Inc. rate (m$^{-3}\cdot$s$^{-1}$)')
legend(leg)

subplot(222)
% figure(3)
% set(gcf,'color','white')
plot(rate_dd(:,2)*1000,rate_dd(:,215))
hold on
plot(rate_sd(:,2)*1000,rate_sd(:,215),':')
plot(rate_ds(:,2)*1000,rate_ds(:,215),'--')
plot(rate_ss(:,2)*1000,rate_ss(:,215),'-.')
plot(rate_wd1(:,2)*1000,rate_wd1(:,215),'--')
plot(rate_wd2(:,2)*1000,rate_wd2(:,215),'-.')
xlabel('Time (ms)')
ylabel('Surf. growth rate (m$^{-3}\cdot$s$^{-1}$)')
legend(leg)

subplot(223)
% figure(3)
% set(gcf,'color','white')
plot(rate_dd(:,2)*1000,sum(rate_dd(:,217:2:end),2))
hold on
plot(rate_sd(:,2)*1000,sum(rate_sd(:,217:2:end),2),':')
plot(rate_ds(:,2)*1000,sum(rate_ds(:,217:2:end),2),'--')
plot(rate_ss(:,2)*1000,sum(rate_ss(:,217:2:end),2),'-.')
plot(rate_wd1(:,2)*1000,sum(rate_wd1(:,217:2:end),2),'--')
plot(rate_wd2(:,2)*1000,sum(rate_wd2(:,217:2:end),2),'-.')
xlabel('Time (ms)')
ylabel('Coag. rate (m$^{-3}\cdot$s$^{-1}$)')
legend(leg)

subplot(224)
% figure(2)
% set(gcf,'color','white')
plot(rate_dd(:,2)*1000,rate_dd(:,213))
hold on
plot(rate_sd(:,2)*1000,rate_sd(:,213),':')
plot(rate_ds(:,2)*1000,rate_ds(:,213),'--')
plot(rate_ss(:,2)*1000,rate_ss(:,213),'-.')
plot(rate_wd1(:,2)*1000,rate_wd1(:,213),'--')
plot(rate_wd2(:,2)*1000,rate_wd2(:,213),'-.')
xlabel('Time (ms)')
ylabel('No inc. rate (m$^{-3}\cdot$s$^{-1}$)')
legend(leg)
