clc, clear, close all

set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

nbins=50;

fold = 'local/cstr_2048SPs_2runs/'; 

part_dd = csvread([fold 'dsa-d/Network(stage1)-part.csv'],1);
part_sd = csvread([fold 'swa-d/Network(stage1)-part.csv'],1);
part_wd1 = csvread([fold 'swa-dn1/Network(stage1)-part.csv'],1);
part_wd2 = csvread([fold 'swa-dn2e2/Network(stage1)-part.csv'],1);

rate_dd = csvread([fold 'dsa-d/Network(stage1)-part-rates.csv'],1);
rate_sd = csvread([fold 'swa-d/Network(stage1)-part-rates.csv'],1);
rate_wd1 = csvread([fold 'swa-dn1/Network(stage1)-part-rates.csv'],1);
rate_wd2 = csvread([fold 'swa-dn2e2/Network(stage1)-part-rates.csv'],1);

psl_dd = csvread([fold 'dsa-d/Network(stage1)-psl(0.5s).csv'],1);
psl_sd = csvread([fold 'swa-d/Network(stage1)-psl(0.5s).csv'],1);
psl_wd1 = csvread([fold 'swa-dn1/Network(stage1)-psl(0.5s).csv'],1);
psl_wd2 = csvread([fold 'swa-dn2e2/Network(stage1)-psl(0.5s).csv'],1);

chem_dd = csvread([fold 'dsa-d/Network(stage1)-chem.csv'],1);
chem_sd = csvread([fold 'swa-d/Network(stage1)-chem.csv'],1);
chem_wd1 = csvread([fold 'swa-dn1/Network(stage1)-chem.csv'],1);
chem_wd2 = csvread([fold 'swa-dn2e2/Network(stage1)-chem.csv'],1);

cput_dd = csvread([fold 'dsa-d/Network(stage1)-cput.csv'],1);
cput_sd = csvread([fold 'swa-d/Network(stage1)-cput.csv'],1);
cput_wd1 = csvread([fold 'swa-dn1/Network(stage1)-cput.csv'],1);
cput_wd2 = csvread([fold 'swa-dn2e2/Network(stage1)-cput.csv'],1);

% part_ds = csvread([fold 'dsa-s/Network(stage1)-part.csv'],1);
% part_ss = csvread([fold 'swa-s/Network(stage1)-part.csv'],1);
% 
% rate_ds = csvread([fold 'dsa-s/Network(stage1)-part-rates.csv'],1);
% rate_ss = csvread([fold 'swa-s/Network(stage1)-part-rates.csv'],1);
% 
% psl_ds = csvread([fold 'dsa-s/Network(stage1)-psl(0.5s).csv'],1);
% psl_ss = csvread([fold 'swa-s/Network(stage1)-psl(0.5s).csv'],1);
% 
% chem_ds = csvread([fold 'dsa-s/Network(stage1)-chem.csv'],1);
% chem_ss = csvread([fold 'swa-s/Network(stage1)-chem.csv'],1);
% 
% cput_ds = csvread([fold 'dsa-s/Network(stage1)-cput.csv'],1);
% cput_ss = csvread([fold 'swa-s/Network(stage1)-cput.csv'],1);

% leg = {'DSA det, orig','SWA det, orig',...
%        'DSA sph, orig','SWA sph, orig',...
%        'SWA det, wt const','SWA det, wt var'};      
% 
leg = {'DSA','CIW SWA, $w_{\textrm{inc}}=1$',...
       'CIW SWA, $w_{\textrm{inc}}=100$',...
       'AIW SWA, $w_{\textrm{inc}}\in\left[1,500\right]$'};

% leg = {'DSA','LAIW SWA, $w_{\textrm{inc}}\in\left[1,500\right]$',...
%        'QAIW SWA, $w_{\textrm{inc}}\in\left[1,500\right]$',...
%        'EAIW SWA, $w_{\textrm{inc}}\in\left[1,500\right]$'};
   
studyid = '_real_1tau_1024SP';
savfigs = 0;
fdir    = 'C:\Users\Astrid\Documents\Projects\Network temperature dependence\IdeasFromHMMeeting\figures\AIWSWA\';
if savfigs
    disp('Saving figures - continuing may overwrite existing!')
    pause 
    disp('Continuing...')
    saveas = @(str)(print([fdir str],'-depsc'));
else
    saveas = @(str)(disp('Not saving'));
end

%% Part

figure(1)
set(gcf,'color','white')
% subplot(231)
plot(part_dd(:,2)*1000,part_dd(:,3))
hold on
plot(part_sd(:,2)*1000,part_sd(:,3),':')
% plot(part_ds(:,2)*1000,part_ds(:,3),'--')
% plot(part_ss(:,2)*1000,part_ss(:,3),'-.')
plot(part_wd1(:,2)*1000,part_wd1(:,3),'--')
plot(part_wd2(:,2)*1000,part_wd2(:,3),'-.')
xlabel('Time (ms)')
ylabel('NSP (-)')
l = legend(leg);
l.Location = 'SouthEast';
l.Interpreter = 'latex';
saveas(['nsp' studyid])

% subplot(232)
figure(2)
set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,5))
hold on
plot(part_sd(:,2)*1000,part_sd(:,5),':')
% plot(part_ds(:,2)*1000,part_ds(:,5),'--')
% plot(part_ss(:,2)*1000,part_ss(:,5),'-.')
plot(part_wd1(:,2)*1000,part_wd1(:,5),'--')
plot(part_wd2(:,2)*1000,part_wd2(:,5),'-.')
xlabel('Time (ms)')
ylabel('M0 (m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['m0' studyid])

% subplot(233)
figure(3)
set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,9)*1e9)
hold on
plot(part_sd(:,2)*1000,part_sd(:,9)*1e9,':')
% plot(part_ds(:,2)*1000,part_ds(:,9)*1e9,'--')
% plot(part_ss(:,2)*1000,part_ss(:,9)*1e9,'-.')
plot(part_wd1(:,2)*1000,part_wd1(:,9)*1e9,'--')
plot(part_wd2(:,2)*1000,part_wd2(:,9)*1e9,'-.')
xlabel('Time (ms)')
ylabel('Collision diamter (nm)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['dcol_ave' studyid])

% subplot(234)
figure(4)
set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,21))
hold on
plot(part_sd(:,2)*1000,part_sd(:,21),':')
% plot(part_ds(:,2)*1000,part_ds(:,21),'--')
% plot(part_ss(:,2)*1000,part_ss(:,21),'-.')
plot(part_wd1(:,2)*1000,part_wd1(:,21),'--')
plot(part_wd2(:,2)*1000,part_wd2(:,21),'-.')
xlabel('Time (ms)')
ylabel('Mass (kg$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['mass_ave' studyid])

% subplot(235)
figure(5)
set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,37))
hold on
plot(part_sd(:,2)*1000,part_sd(:,37),':')
plot(part_wd1(:,2)*1000,part_wd1(:,37),'--')
plot(part_wd2(:,2)*1000,part_wd2(:,37),'-.')
xlabel('Time (ms)')
ylabel('Number primaries per particle (-)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['npri_ave' studyid])

% subplot(236)
figure(6)
set(gcf,'color','white')
plot(part_dd(:,2)*1000,part_dd(:,39)*1e9)
hold on
plot(part_sd(:,2)*1000,part_sd(:,39)*1e9,':')
plot(part_wd1(:,2)*1000,part_wd1(:,39)*1e9,'--')
plot(part_wd2(:,2)*1000,part_wd2(:,39)*1e9,'-.')
xlabel('Time (ms)')
ylabel('Primary average diamter (nm)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['dpri_ave' studyid])

%% PSL

figure(7)
set(gcf,'color','white')
[n_d,x_d] = histwc(psl_dd(:,3),psl_dd(:,9),nbins);
[n_s,x_s] = histwc(psl_sd(:,3),psl_sd(:,9),nbins);
% [n_ds,x_ds] = hist(psl_ds(:,3));
% [n_ss,x_ss] = hist(psl_ss(:,3));
[n_d2,x_d2] = histwc(psl_wd1(:,3),psl_wd1(:,9),nbins);
[n_s2,x_s2] = histwc(psl_wd2(:,3),psl_wd2(:,9),nbins);
% subplot(231)
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
% plot(x_ds,n_ds/max(n_ds),'--',x_ss,n_ss/max(n_ss),'-.')
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Collision diameter (nm)')
ylabel('Count divided by max. count (-)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['dcol' studyid])

[n_d,x_d] = histwc(psl_dd(:,5),psl_dd(:,9),nbins);
[n_s,x_s] = histwc(psl_sd(:,5),psl_sd(:,9),nbins);
[n_d2,x_d2] = histwc(psl_wd1(:,5),psl_wd1(:,9),nbins);
[n_s2,x_s2] = histwc(psl_wd2(:,5),psl_wd2(:,9),nbins);
% subplot(232)
figure(8)
set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Surface area (cm$^{2}$)')
ylabel('Count divided by max. count (-)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['sa' studyid])

[n_d,x_d] = histwc(psl_dd(:,12),psl_dd(:,9),nbins);
[n_s,x_s] = histwc(psl_sd(:,12),psl_sd(:,9),nbins);
[n_d2,x_d2] = histwc(psl_wd1(:,12),psl_wd1(:,9),nbins);
[n_s2,x_s2] = histwc(psl_wd2(:,12),psl_wd2(:,9),nbins);
% subplot(233)
figure(9)
set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Number of primaries (nm)')
ylabel('Count divided by max. count (-)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['npri' studyid])

[n_d,x_d] = histwc(psl_dd(:,13),psl_dd(:,9),nbins);
[n_s,x_s] = histwc(psl_sd(:,13),psl_sd(:,9),nbins);
[n_d2,x_d2] = histwc(psl_wd1(:,13),psl_wd1(:,9),nbins);
[n_s2,x_s2] = histwc(psl_wd2(:,13),psl_wd2(:,9),nbins);
% subplot(234)
figure(10)
set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Primary diameter (nm)')
ylabel('Count divided by max. count (-)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['dprim' studyid])

[n_d,x_d] = histwc(psl_dd(:,16),psl_dd(:,9),nbins);
[n_s,x_s] = histwc(psl_sd(:,16),psl_sd(:,9),nbins);
[n_d2,x_d2] = histwc(psl_wd1(:,16),psl_wd1(:,9),nbins);
[n_s2,x_s2] = histwc(psl_wd2(:,16),psl_wd2(:,9),nbins);
% subplot(235)
figure(11)
set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Std of primary diameter (nm)')
ylabel('Count divided by max. count (-)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['std_part_diam' studyid])

[n_d,x_d] = histwc(psl_dd(:,8),psl_dd(:,9),nbins);
[n_s,x_s] = histwc(psl_sd(:,8),psl_sd(:,9),nbins);
[n_d2,x_d2] = histwc(psl_wd1(:,8),psl_wd1(:,9),nbins);
[n_s2,x_s2] = histwc(psl_wd2(:,8),psl_wd2(:,9),nbins);
% subplot(236)
figure(12)
set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Particle age (s)')
ylabel('Count divided by max. count (-)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
ax = gca;
ax.XLim(1) = 0;
saveas(['part_age' studyid])

[n_d,x_d] = hist(psl_dd(:,9));
[n_s,x_s] = hist(psl_sd(:,9));
[n_d2,x_d2] = hist(psl_wd1(:,9));
[n_s2,x_s2] = hist(psl_wd2(:,9));
% subplot(236)
figure(121)
set(gcf,'color','white')
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Statistical weight (-)')
ylabel('Count divided by max. count (-)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
ax = gca;
ax.XLim(1) = 0;
saveas(['stat_weight' studyid])

%% Chem

figure(13)
set(gcf,'color','white')
% subplot(231)
plot(chem_dd(:,2)*1000,chem_dd(:,3)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,3)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,3)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,3)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,3)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,3)*1e6,'-.')
xlabel('Time (ms)')
ylabel('TiCl$_4$ conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['ticl4' studyid])

% subplot(232)
figure(14)
set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,39)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,39)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,39)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,39)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,39)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,39)*1e6,'-.')
xlabel('Time (ms)')
ylabel('O$_2$ conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['o2' studyid])

% subplot(233)
figure(15)
set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,45)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,45)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,45)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,45)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,45)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,45)*1e6,'-.')
xlabel('Time (ms)')
ylabel('Cl$_2$ conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['cl2' studyid])

% subplot(234)
figure(16)
set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,19)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,19)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,19)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,19)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,19)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,19)*1e6,'-.')
xlabel('Time (ms)')
ylabel('TiO$_2$Cl$_3$ conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['tio2cl3' studyid])

% subplot(235)
figure(17)
set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,57)*1e6)
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,57)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,57)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,57)*1e6,'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,57)*1e6,'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,57)*1e6,'-.')
xlabel('Time (ms)')
ylabel('Ar conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['ar' studyid])

% subplot(236)
figure(18)
set(gcf,'color','white')
plot(chem_dd(:,2)*1000,chem_dd(:,61))
hold on
plot(chem_sd(:,2)*1000,chem_sd(:,61),':')
% plot(chem_ds(:,2)*1000,chem_ds(:,61),'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,61),'-.')
plot(chem_wd1(:,2)*1000,chem_wd1(:,61),'--')
plot(chem_wd2(:,2)*1000,chem_wd2(:,61),'-.')
xlabel('Time (ms)')
ylabel('Temperature (K)')
l = legend(leg);
l.Location = 'SouthEast';
l.Interpreter = 'latex';
saveas(['temp' studyid])

%% CPUT

figure(19)
set(gcf,'color','white')
plot(cput_dd(:,2)*1000,cput_dd(:,3)/60)
hold on
plot(cput_sd(:,2)*1000,cput_sd(:,3)/60,':')
% plot(cput_ds(:,2)*1000,cput_ds(:,3)/60,'--')
% plot(cput_ss(:,2)*1000,cput_ss(:,3)/60,'-.')
plot(cput_wd1(:,2)*1000,cput_wd1(:,3)/60,'--')
plot(cput_wd2(:,2)*1000,cput_wd2(:,3)/60,'-.')
xlabel('Time (ms)')
ylabel('Solver time (min)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['solver_times' studyid])

%% Rates

figure(20)
set(gcf,'color','white')
% subplot(221)
plot(rate_dd(:,2)*1000,sum(rate_dd(:,3:2:211),2))
hold on
plot(rate_sd(:,2)*1000,sum(rate_sd(:,3:2:211),2),':')
% plot(rate_ds(:,2)*1000,sum(rate_ds(:,3:2:211),2),'--')
% plot(rate_ss(:,2)*1000,sum(rate_ss(:,3:2:211),2),'-.')
plot(rate_wd1(:,2)*1000,sum(rate_wd1(:,3:2:211),2),'--')
plot(rate_wd2(:,2)*1000,sum(rate_wd2(:,3:2:211),2),'-.')
xlabel('Time (ms)')
ylabel('Inc. rate (m$^{-3}\cdot$s$^{-1}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['inc_rate' studyid])

% subplot(222)
figure(21)
set(gcf,'color','white')
plot(rate_dd(:,2)*1000,rate_dd(:,215))
hold on
plot(rate_sd(:,2)*1000,rate_sd(:,215),':')
% plot(rate_ds(:,2)*1000,rate_ds(:,215),'--')
% plot(rate_ss(:,2)*1000,rate_ss(:,215),'-.')
plot(rate_wd1(:,2)*1000,rate_wd1(:,215),'--')
plot(rate_wd2(:,2)*1000,rate_wd2(:,215),'-.')
xlabel('Time (ms)')
ylabel('Surf. growth rate (m$^{-3}\cdot$s$^{-1}$)')
l = legend(leg);
l.Location = 'northWest';
l.Interpreter = 'latex';
saveas(['sg_rate' studyid])

% subplot(223)
figure(22)
set(gcf,'color','white')
plot(rate_dd(:,2)*1000,sum(rate_dd(:,217:2:end),2))
hold on
plot(rate_sd(:,2)*1000,sum(rate_sd(:,217:2:end),2),':')
% plot(rate_ds(:,2)*1000,sum(rate_ds(:,217:2:end),2),'--')
% plot(rate_ss(:,2)*1000,sum(rate_ss(:,217:2:end),2),'-.')
plot(rate_wd1(:,2)*1000,sum(rate_wd1(:,217:2:end),2),'--')
plot(rate_wd2(:,2)*1000,sum(rate_wd2(:,217:2:end),2),'-.')
xlabel('Time (ms)')
ylabel('Coag. rate (m$^{-3}\cdot$s$^{-1}$)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['coag_rate' studyid])

% subplot(224)
figure(23)
set(gcf,'color','white')
plot(rate_dd(:,2)*1000,rate_dd(:,213))
hold on
plot(rate_sd(:,2)*1000,rate_sd(:,213),':')
% plot(rate_ds(:,2)*1000,rate_ds(:,213),'--')
% plot(rate_ss(:,2)*1000,rate_ss(:,213),'-.')
plot(rate_wd1(:,2)*1000,rate_wd1(:,213),'--')
plot(rate_wd2(:,2)*1000,rate_wd2(:,213),'-.')
xlabel('Time (ms)')
ylabel('No inc. rate (m$^{-3}\cdot$s$^{-1}$)')
l = legend(leg);
l.Location = 'SouthWest';
l.Interpreter = 'latex';
saveas(['no_inc_rate' studyid])
