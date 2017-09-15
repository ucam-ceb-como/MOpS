clc, 
clear, close all

%% Setup

% set defaults
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

% for finding files
basedir  = 'vienna/vart_2048SPs_10runs/16384SPs_10runs/'; 
filedir1 = 'dsa-d/';
filedir2 = 'swa-d/';
filedir3 = 'swa-dn1/';
filedir4 = 'swa-dn2/';
filebase = 'PP';

% for plotting psl data
psltime  = '0.5'; 
npslbins = 25;

% for saving images
studyid  = '_batch_10run_16384SP';
projpath = 'C:\Users\Astrid\Documents\Projects\';
projdir  = 'Network temperature dependence\IdeasFromHMMeeting\';
imagedir = 'figures\AIWSWA\';
savefigs = 0;
legtype  = 'algorithm_comparison'; % or weightfun_comparison

% load all data
part_dsa = csvread([basedir filedir1 filebase '-part.csv'],1);
part_swa = csvread([basedir filedir2 filebase '-part.csv'],1);
part_ciw = csvread([basedir filedir3 filebase '-part.csv'],1);
part_aiw = csvread([basedir filedir4 filebase '-part.csv'],1);

rate_dsa = csvread([basedir filedir1 filebase '-part-rates.csv'],1);
rate_swa = csvread([basedir filedir2 filebase '-part-rates.csv'],1);
rate_ciw = csvread([basedir filedir3 filebase '-part-rates.csv'],1);
rate_aiw = csvread([basedir filedir4 filebase '-part-rates.csv'],1);

psl_dsa = csvread([basedir filedir1 filebase '-psl(' psltime 's).csv'],1);
psl_swa = csvread([basedir filedir2 filebase '-psl(' psltime 's).csv'],1);
psl_ciw = csvread([basedir filedir3 filebase '-psl(' psltime 's).csv'],1);
psl_aiw = csvread([basedir filedir4 filebase '-psl(' psltime 's).csv'],1);

chem_dsa = csvread([basedir filedir1 filebase '-chem.csv'],1);
chem_swa = csvread([basedir filedir2 filebase '-chem.csv'],1);
chem_ciw = csvread([basedir filedir3 filebase '-chem.csv'],1);
chem_aiw = csvread([basedir filedir4 filebase '-chem.csv'],1);

cput_dsa = csvread([basedir filedir1 filebase '-cput.csv'],1);
cput_swa = csvread([basedir filedir2 filebase '-cput.csv'],1);
cput_ciw = csvread([basedir filedir3 filebase '-cput.csv'],1);
cput_aiw = csvread([basedir filedir4 filebase '-cput.csv'],1);

% part_ds = csvread([fold 'dsa-s/' fname '-part.csv'],1);
% part_ss = csvread([fold 'swa-s/' fname '-part.csv'],1);
% 
% rate_ds = csvread([fold 'dsa-s/' fname '-part-rates.csv'],1);
% rate_ss = csvread([fold 'swa-s/' fname '-part-rates.csv'],1);
% 
% psl_ds = csvread([fold 'dsa-s/' fname '-psl(0.5s).csv'],1);
% psl_ss = csvread([fold 'swa-s/' fname '-psl(0.5s).csv'],1);
% 
% chem_ds = csvread([fold 'dsa-s/' fname '-chem.csv'],1);
% chem_ss = csvread([fold 'swa-s/' fname '-chem.csv'],1);
% 
% cput_ds = csvread([fold 'dsa-s/' fname '-cput.csv'],1);
% cput_ss = csvread([fold 'swa-s/' fname '-cput.csv'],1);

% choose appropriate legend
if strcmp(legtype,'algorithm_comparison')
    leg = {'DSA','CIW SWA, $w_{\textrm{inc}}=1$',...
           'CIW SWA, $w_{\textrm{inc}}=100$',...
           'EAIW SWA, $w_{\textrm{inc}}\in\left[1,500\right]$'};
elseif strcmp(legtype,'weightfun_comparison')
    leg = {'DSA','LAIW SWA, $w_{\textrm{inc}}\in\left[1,500\right]$',...
           'QAIW SWA, $w_{\textrm{inc}}\in\left[1,500\right]$',...
           'EAIW SWA, $w_{\textrm{inc}}\in\left[1,500\right]$'};
end

% initialise saveas function if saving figures, otherwise do nothing here
if savefigs
    disp('Saving figures - continuing may overwrite existing!')
    pause 
    disp('Continuing...')
    saveas = @(str)(print([projpath projdir imagedir str],'-depsc'));
else
    saveas = @(str)(disp('Not saving'));
end

%% Part

figure(1)
set(gcf,'color','white')
% subplot(231)
plot(part_dsa(:,2)*1000,part_dsa(:,3))
hold on
plot(part_swa(:,2)*1000,part_swa(:,3),':')
% plot(part_ds(:,2)*1000,part_ds(:,3),'--')
% plot(part_ss(:,2)*1000,part_ss(:,3),'-.')
plot(part_ciw(:,2)*1000,part_ciw(:,3),'--')
plot(part_aiw(:,2)*1000,part_aiw(:,3),'-.')
xlabel('Time (ms)')
ylabel('NSP (-)')
l = legend(leg);
l.Location = 'SouthEast';
l.Interpreter = 'latex';
saveas(['nsp' studyid])

% subplot(232)
figure(2)
set(gcf,'color','white')
plot(part_dsa(:,2)*1000,part_dsa(:,5))
hold on
plot(part_swa(:,2)*1000,part_swa(:,5),':')
% plot(part_ds(:,2)*1000,part_ds(:,5),'--')
% plot(part_ss(:,2)*1000,part_ss(:,5),'-.')
plot(part_ciw(:,2)*1000,part_ciw(:,5),'--')
plot(part_aiw(:,2)*1000,part_aiw(:,5),'-.')
xlabel('Time (ms)')
ylabel('M0 (m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['m0' studyid])

% subplot(233)
figure(3)
set(gcf,'color','white')
plot(part_dsa(:,2)*1000,part_dsa(:,9)*1e9)
hold on
plot(part_swa(:,2)*1000,part_swa(:,9)*1e9,':')
% plot(part_ds(:,2)*1000,part_ds(:,9)*1e9,'--')
% plot(part_ss(:,2)*1000,part_ss(:,9)*1e9,'-.')
plot(part_ciw(:,2)*1000,part_ciw(:,9)*1e9,'--')
plot(part_aiw(:,2)*1000,part_aiw(:,9)*1e9,'-.')
xlabel('Time (ms)')
ylabel('Collision diamter (nm)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['dcol_ave' studyid])

% subplot(234)
figure(4)
set(gcf,'color','white')
plot(part_dsa(:,2)*1000,part_dsa(:,21))
hold on
plot(part_swa(:,2)*1000,part_swa(:,21),':')
% plot(part_ds(:,2)*1000,part_ds(:,21),'--')
% plot(part_ss(:,2)*1000,part_ss(:,21),'-.')
plot(part_ciw(:,2)*1000,part_ciw(:,21),'--')
plot(part_aiw(:,2)*1000,part_aiw(:,21),'-.')
xlabel('Time (ms)')
ylabel('Mass (kg$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['mass_ave' studyid])

% subplot(235)
figure(5)
set(gcf,'color','white')
plot(part_dsa(:,2)*1000,part_dsa(:,37))
hold on
plot(part_swa(:,2)*1000,part_swa(:,37),':')
plot(part_ciw(:,2)*1000,part_ciw(:,37),'--')
plot(part_aiw(:,2)*1000,part_aiw(:,37),'-.')
xlabel('Time (ms)')
ylabel('Number primaries per particle (-)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['npri_ave' studyid])

% subplot(236)
figure(6)
set(gcf,'color','white')
plot(part_dsa(:,2)*1000,part_dsa(:,39)*1e9)
hold on
plot(part_swa(:,2)*1000,part_swa(:,39)*1e9,':')
plot(part_ciw(:,2)*1000,part_ciw(:,39)*1e9,'--')
plot(part_aiw(:,2)*1000,part_aiw(:,39)*1e9,'-.')
xlabel('Time (ms)')
ylabel('Primary average diamter (nm)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['dpri_ave' studyid])

%% PSL

[n_d,x_d] = histwc(psl_dsa(:,3),psl_dsa(:,9),npslbins);
[n_s,x_s] = histwc(psl_swa(:,3),psl_swa(:,9),npslbins);
[n_d2,x_d2] = histwc(psl_ciw(:,3),psl_ciw(:,9),npslbins);
[n_s2,x_s2] = histwc(psl_aiw(:,3),psl_aiw(:,9),npslbins);
figure(7)
set(gcf,'color','white')
% subplot(231)
plot(x_d,n_d/max(n_d),x_s,n_s/max(n_s),':')
hold on
plot(x_d2,n_d2/max(n_d2),'--')
plot(x_s2,n_s2/max(n_s2),'-.')
xlabel('Collision diameter (nm)')
ylabel('Count divided by max. count (-)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
ax = gca;
ax.XLim(1) = 0.49;
saveas(['dcol' studyid])

[n_d,x_d] = histwc(psl_dsa(:,5),psl_dsa(:,9),npslbins);
[n_s,x_s] = histwc(psl_swa(:,5),psl_swa(:,9),npslbins);
[n_d2,x_d2] = histwc(psl_ciw(:,5),psl_ciw(:,9),npslbins);
[n_s2,x_s2] = histwc(psl_aiw(:,5),psl_aiw(:,9),npslbins);
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
ax = gca;
ax.XLim(1) = 0;
saveas(['sa' studyid])

[n_d,x_d] = histwc(psl_dsa(:,12),psl_dsa(:,9),npslbins);
[n_s,x_s] = histwc(psl_swa(:,12),psl_swa(:,9),npslbins);
[n_d2,x_d2] = histwc(psl_ciw(:,12),psl_ciw(:,9),npslbins);
[n_s2,x_s2] = histwc(psl_aiw(:,12),psl_aiw(:,9),npslbins);
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
ax = gca;
ax.XLim(1) = 0;
saveas(['npri' studyid])

[n_d,x_d] = histwc(psl_dsa(:,13),psl_dsa(:,9),npslbins);
[n_s,x_s] = histwc(psl_swa(:,13),psl_swa(:,9),npslbins);
[n_d2,x_d2] = histwc(psl_ciw(:,13),psl_ciw(:,9),npslbins);
[n_s2,x_s2] = histwc(psl_aiw(:,13),psl_aiw(:,9),npslbins);
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
ax = gca;
ax.XLim(1) = 0.49;
saveas(['dprim' studyid])

[n_d,x_d] = histwc(psl_dsa(:,16),psl_dsa(:,9),npslbins);
[n_s,x_s] = histwc(psl_swa(:,16),psl_swa(:,9),npslbins);
[n_d2,x_d2] = histwc(psl_ciw(:,16),psl_ciw(:,9),npslbins);
[n_s2,x_s2] = histwc(psl_aiw(:,16),psl_aiw(:,9),npslbins);
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
ax = gca;
ax.XLim(1) = 0;
saveas(['std_part_diam' studyid])

[n_d,x_d] = histwc(psl_dsa(:,8),psl_dsa(:,9),npslbins);
[n_s,x_s] = histwc(psl_swa(:,8),psl_swa(:,9),npslbins);
[n_d2,x_d2] = histwc(psl_ciw(:,8),psl_ciw(:,9),npslbins);
[n_s2,x_s2] = histwc(psl_aiw(:,8),psl_aiw(:,9),npslbins);
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

[n_d,x_d] = hist(psl_dsa(:,9));
[n_s,x_s] = hist(psl_swa(:,9));
[n_d2,x_d2] = hist(psl_ciw(:,9));
[n_s2,x_s2] = hist(psl_aiw(:,9));
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
plot(chem_dsa(:,2)*1000,chem_dsa(:,3)*1e6)
hold on
plot(chem_swa(:,2)*1000,chem_swa(:,3)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,3)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,3)*1e6,'-.')
plot(chem_ciw(:,2)*1000,chem_ciw(:,3)*1e6,'--')
plot(chem_aiw(:,2)*1000,chem_aiw(:,3)*1e6,'-.')
xlabel('Time (ms)')
ylabel('TiCl$_4$ conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['ticl4' studyid])

% subplot(232)
figure(14)
set(gcf,'color','white')
plot(chem_dsa(:,2)*1000,chem_dsa(:,39)*1e6)
hold on
plot(chem_swa(:,2)*1000,chem_swa(:,39)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,39)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,39)*1e6,'-.')
plot(chem_ciw(:,2)*1000,chem_ciw(:,39)*1e6,'--')
plot(chem_aiw(:,2)*1000,chem_aiw(:,39)*1e6,'-.')
xlabel('Time (ms)')
ylabel('O$_2$ conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['o2' studyid])

% subplot(233)
figure(15)
set(gcf,'color','white')
plot(chem_dsa(:,2)*1000,chem_dsa(:,45)*1e6)
hold on
plot(chem_swa(:,2)*1000,chem_swa(:,45)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,45)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,45)*1e6,'-.')
plot(chem_ciw(:,2)*1000,chem_ciw(:,45)*1e6,'--')
plot(chem_aiw(:,2)*1000,chem_aiw(:,45)*1e6,'-.')
xlabel('Time (ms)')
ylabel('Cl$_2$ conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['cl2' studyid])

% subplot(234)
figure(16)
set(gcf,'color','white')
plot(chem_dsa(:,2)*1000,chem_dsa(:,19)*1e6)
hold on
plot(chem_swa(:,2)*1000,chem_swa(:,19)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,19)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,19)*1e6,'-.')
plot(chem_ciw(:,2)*1000,chem_ciw(:,19)*1e6,'--')
plot(chem_aiw(:,2)*1000,chem_aiw(:,19)*1e6,'-.')
xlabel('Time (ms)')
ylabel('TiO$_2$Cl$_3$ conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['tio2cl3' studyid])

% subplot(235)
figure(17)
set(gcf,'color','white')
plot(chem_dsa(:,2)*1000,chem_dsa(:,57)*1e6)
hold on
plot(chem_swa(:,2)*1000,chem_swa(:,57)*1e6,':')
% plot(chem_ds(:,2)*1000,chem_ds(:,57)*1e6,'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,57)*1e6,'-.')
plot(chem_ciw(:,2)*1000,chem_ciw(:,57)*1e6,'--')
plot(chem_aiw(:,2)*1000,chem_aiw(:,57)*1e6,'-.')
xlabel('Time (ms)')
ylabel('Ar conc. (mol$\cdot$m$^{-3}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['ar' studyid])

% subplot(236)
figure(18)
set(gcf,'color','white')
plot(chem_dsa(:,2)*1000,chem_dsa(:,61))
hold on
plot(chem_swa(:,2)*1000,chem_swa(:,61),':')
% plot(chem_ds(:,2)*1000,chem_ds(:,61),'--')
% plot(chem_ss(:,2)*1000,chem_ss(:,61),'-.')
plot(chem_ciw(:,2)*1000,chem_ciw(:,61),'--')
plot(chem_aiw(:,2)*1000,chem_aiw(:,61),'-.')
xlabel('Time (ms)')
ylabel('Temperature (K)')
l = legend(leg);
l.Location = 'SouthEast';
l.Interpreter = 'latex';
saveas(['temp' studyid])

%% CPUT

figure(19)
set(gcf,'color','white')
plot(cput_dsa(:,2)*1000,cput_dsa(:,3)/60)
hold on
plot(cput_swa(:,2)*1000,cput_swa(:,3)/60,':')
% plot(cput_ds(:,2)*1000,cput_ds(:,3)/60,'--')
% plot(cput_ss(:,2)*1000,cput_ss(:,3)/60,'-.')
plot(cput_ciw(:,2)*1000,cput_ciw(:,3)/60,'--')
plot(cput_aiw(:,2)*1000,cput_aiw(:,3)/60,'-.')
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
plot(rate_dsa(:,2)*1000,sum(rate_dsa(:,3:2:211),2))
hold on
plot(rate_swa(:,2)*1000,sum(rate_swa(:,3:2:211),2),':')
% plot(rate_ds(:,2)*1000,sum(rate_ds(:,3:2:211),2),'--')
% plot(rate_ss(:,2)*1000,sum(rate_ss(:,3:2:211),2),'-.')
plot(rate_ciw(:,2)*1000,sum(rate_ciw(:,3:2:211),2),'--')
plot(rate_aiw(:,2)*1000,sum(rate_aiw(:,3:2:211),2),'-.')
xlabel('Time (ms)')
ylabel('Inc. rate (m$^{-3}\cdot$s$^{-1}$)')
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
saveas(['inc_rate' studyid])

% subplot(222)
figure(21)
set(gcf,'color','white')
plot(rate_dsa(:,2)*1000,rate_dsa(:,215))
hold on
plot(rate_swa(:,2)*1000,rate_swa(:,215),':')
% plot(rate_ds(:,2)*1000,rate_ds(:,215),'--')
% plot(rate_ss(:,2)*1000,rate_ss(:,215),'-.')
plot(rate_ciw(:,2)*1000,rate_ciw(:,215),'--')
plot(rate_aiw(:,2)*1000,rate_aiw(:,215),'-.')
xlabel('Time (ms)')
ylabel('Surf. growth rate (m$^{-3}\cdot$s$^{-1}$)')
l = legend(leg);
l.Location = 'northWest';
l.Interpreter = 'latex';
saveas(['sg_rate' studyid])

% subplot(223)
figure(22)
set(gcf,'color','white')
plot(rate_dsa(:,2)*1000,sum(rate_dsa(:,217:2:end),2))
hold on
plot(rate_swa(:,2)*1000,sum(rate_swa(:,217:2:end),2),':')
% plot(rate_ds(:,2)*1000,sum(rate_ds(:,217:2:end),2),'--')
% plot(rate_ss(:,2)*1000,sum(rate_ss(:,217:2:end),2),'-.')
plot(rate_ciw(:,2)*1000,sum(rate_ciw(:,217:2:end),2),'--')
plot(rate_aiw(:,2)*1000,sum(rate_aiw(:,217:2:end),2),'-.')
xlabel('Time (ms)')
ylabel('Coag. rate (m$^{-3}\cdot$s$^{-1}$)')
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
saveas(['coag_rate' studyid])

% subplot(224)
figure(23)
set(gcf,'color','white')
plot(rate_dsa(:,2)*1000,rate_dsa(:,213))
hold on
plot(rate_swa(:,2)*1000,rate_swa(:,213),':')
% plot(rate_ds(:,2)*1000,rate_ds(:,213),'--')
% plot(rate_ss(:,2)*1000,rate_ss(:,213),'-.')
plot(rate_ciw(:,2)*1000,rate_ciw(:,213),'--')
plot(rate_aiw(:,2)*1000,rate_aiw(:,213),'-.')
xlabel('Time (ms)')
ylabel('No inc. rate (m$^{-3}\cdot$s$^{-1}$)')
l = legend(leg);
l.Location = 'SouthWest';
l.Interpreter = 'latex';
saveas(['no_inc_rate' studyid])
