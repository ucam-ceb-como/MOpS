clc, 
% clear, close all

%% Setup
% set defaults
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

% for finding files
basedir  = 'vienna/real/real_inc_factor/justforfun/';
filedir  = 'swa/';
filecmp  = 'swa/';
filebase = 'Network(stage1)';

% for plotting psl data
psltime  = '0.021';
npslbins = 10;

% for saving images
studyid  = '';
projpath = 'C:\Users\Astrid\Documents\Projects\';
projdir  = 'Network temperature dependence\IdeasFromHMMeeting\';
imagedir = 'figures\';
savefigs = 0;
leg_vals = '';

% load all data
m0dat = csvread([basedir filecmp filebase '-part.csv'],1);
part = csvread([basedir filedir filebase '-part.csv'],1);
rate = csvread([basedir filedir filebase '-part-rates.csv'],1);
pslf = csvread([basedir filedir filebase '-psl(' psltime 's).csv'],1);
chem = csvread([basedir filedir filebase '-chem.csv'],1);
cput = csvread([basedir filedir filebase '-cput.csv'],1);

% initialise saveas function if saving figures, otherwise do nothing here
if savefigs
    disp('Saving figures - continuing may overwrite existing!')
    pause 
    disp('Continuing...')
    saveas = @(str)(print([projpath projdir imagedir str],'-depsc'));
else
    saveas = @(str)(disp('Not saving'));
end

% initialise addlegend function 
if leg_vals ~= ''
    addlegend = @(leg_entries)(legend(leg_entries));
else
    addlegend = @(leg_entries)(0);
end


%% Part

c1 = [0.65 0.65 0.65];

figure(1)
set(gcf,'color','white')
subplot(231)
plot(part(:,2)*1000,part(:,3))
hold on
xlabel('Time (ms)')
ylabel('NSP (-)')
addlegend(leg_vals);
saveas(['nsp' studyid])

subplot(232)
% figure(2)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,5))
hold on
% plot(part(:,2)*1000,part(:,5)+part(:,6),'--','Color',c1)
% plot(part(:,2)*1000,part(:,5)-part(:,6),'--','Color',c1)
xlabel('Time (ms)')
ylabel('M0 (m$^{-3}$)')
addlegend(leg_vals);
saveas(['m0' studyid])

subplot(233)
% figure(3)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,9)*1e9)
hold on
xlabel('Time (ms)')
ylabel('Collision diameter (nm)')
addlegend(leg_vals);
saveas(['dcol_ave' studyid])

subplot(234)
% figure(4)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,21))
hold on
% plot(part(:,2)*1000,part(:,21)+part(:,22),'--','Color',c1)
% plot(part(:,2)*1000,part(:,21)-part(:,22),'--','Color',c1)
xlabel('Time (ms)')
ylabel('Mass (kg$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['mass_ave' studyid])

subplot(235)
% figure(5)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,37))
hold on
xlabel('Time (ms)')
ylabel('Number primaries per particle (-)')
addlegend(leg_vals);
saveas(['npri_ave' studyid])

subplot(236)
% figure(6)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,39)*1e9)
hold on
xlabel('Time (ms)')
ylabel('Primary average diameter (nm)')
addlegend(leg_vals);
saveas(['dpri_ave' studyid])

figure(60)
set(gcf,'color','white')
semilogy(part(:,2)*1000,part(:,33)*1e6)
hold on
xlabel('Time (ms)')
ylabel('Total rutile (m$^{-3}$)')
addlegend(leg_vals);
saveas(['rut_ave' studyid])

figure(90)
set(gcf,'color','white')
% subplot(121)
% plot(part(:,2)*1000,part(:,5))
% hold on
% plot(part(:,2)*1000,part(:,5)+part(:,6),'--','Color',c1)
% plot(part(:,2)*1000,part(:,5)-part(:,6),'--','Color',c1)
% xlabel('Time (ms)')
% ylabel('M0 (m$^{-3}$)')
% addlegend(leg_vals);
% subplot(122)
plot(part(:,2)*1000,abs(part(:,5)-m0dat(:,5))./m0dat(:,5))
hold on
xlabel('Time (ms)')
ylabel('Relative difference in m0 (-)')
saveas(['m0_err' studyid])

figure(91)
set(gcf,'color','white')
plot(part(:,2)*1000,abs(part(:,39)-m0dat(:,39))./m0dat(:,39))
hold on
xlabel('Time (ms)')
ylabel('Relative difference in primary diameter (-)')
saveas(['dpri_err' studyid])

figure(92)
set(gcf,'color','white')
plot(part(:,2)*1000,abs(part(:,21)-m0dat(:,21))./m0dat(:,21))
hold on
xlabel('Time (ms)')
ylabel('Relative difference in m1 (-)')
saveas(['m1_err' studyid])

%% PSL

figure(2)
set(gcf,'color','white')
[n_d,x_d] = histwc(pslf(:,3),pslf(:,9),npslbins);
subplot(231)
plot(x_d,n_d/max(n_d))
hold on
xlabel('Collision diameter (nm)')
ylabel('Count divided by max. count (-)')
addlegend(leg_vals);
% ax = gca;
% ax.XLim(1) = 0.49;
saveas(['dcol_psd' studyid])

[n_d,x_d] = histwc(pslf(:,5),pslf(:,9),npslbins);
subplot(232)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Surface area (cm$^{2}$)')
ylabel('Count divided by max. count (-)')
addlegend(leg_vals);
% ax = gca;
% ax.XLim(1) = 0.0;
saveas(['sa_psd' studyid])

[n_d,x_d] = histwc(pslf(:,12),pslf(:,9),npslbins);
subplot(233)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Number of primaries (nm)')
ylabel('Count divided by max. count (-)')
addlegend(leg_vals);
ax = gca;
ax.XLim(1) = 0.0;
saveas(['npri_psd' studyid])

[n_d,x_d] = histwc(pslf(:,13),pslf(:,9),npslbins);
subplot(234)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Primary diameter (nm)')
ylabel('Count divided by max. count (-)')
addlegend(leg_vals);
% ax = gca;
% ax.XLim(1) = 0.49;
saveas(['dpri_psd' studyid])

[n_d,x_d] = histwc(pslf(:,16),pslf(:,9),npslbins);
subplot(235)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Std of primary diameter (nm)')
ylabel('Count divided by max. count (-)')
addlegend(leg_vals);
ax = gca;
ax.XLim(1) = 0.0;
saveas(['std_dpri_psd' studyid])

[n_d,x_d] = histwc(pslf(:,8),pslf(:,9),npslbins);
subplot(236)
% figure()
% set(gcf,'color','white')
plot(x_d,n_d/max(n_d))
hold on
xlabel('Particle age (s)')
ylabel('Count divided by max. count (-)')
addlegend(leg_vals);
ax = gca;
ax.XLim(1) = 0.0;
saveas(['page_psd' studyid])

figure(20)
set(gcf,'color','white')
[n_d,x_d] = hist(pslf(:,9));
% subplot(231)
plot(x_d,n_d/max(n_d))
hold on
xlabel('Weight (-)')
ylabel('Count divided by max. count (-)')
addlegend(leg_vals);
ax = gca;
ax.XLim(1) = 0.0;
saveas(['weight_psd' studyid])

figure(30)
set(gcf,'color','white')
% subplot(131)
[n_d,x_d] = histwc(pslf(:,11),pslf(:,9),npslbins);
plot(x_d,n_d/max(n_d))
hold on
xlabel('Rutile units (-)')
ylabel('Count divided by max. count (-)')
addlegend(leg_vals);
ax = gca;
ax.XLim(1) = 0.0;
saveas(['rutile_psd' studyid])

%% Chem

figure(3)
set(gcf,'color','white')
subplot(231)
plot(chem(:,2)*1000,chem(:,3)*1e6)
hold on
xlabel('Time (ms)')
ylabel('TiCl$_4$ conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['ticl4' studyid])

subplot(232)
% figure(2)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,39)*1e6)
hold on
xlabel('Time (ms)')
ylabel('O$_2$ conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['o2' studyid])

subplot(233)
% figure(3)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,45)*1e6)
hold on
xlabel('Time (ms)')
ylabel('Cl$_2$ conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['cl2' studyid])

subplot(234)
% figure(4)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,19)*1e6)
hold on
xlabel('Time (ms)')
ylabel('TiO$_2$Cl$_3$ conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['tio2cl3' studyid])

subplot(235)
% figure(5)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,57)*1e6)
hold on
xlabel('Time (ms)')
ylabel('Ar conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['ar' studyid])

subplot(236)
% figure(6)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,61))
hold on
xlabel('Time (ms)')
ylabel('Temperature (K)')
addlegend(leg_vals);
saveas(['temp' studyid])

%% CPUT

cadj = [0;(cput(2:end,3)-cput(2,3))/3600];
figure(4)
set(gcf,'color','white')
plot(cput(:,2)*1000,cadj)
hold on
xlabel('Time (ms)')
ylabel('Solver time (min)')
addlegend(leg_vals);
saveas(['cput' studyid])

%% Rates

figure(5)
set(gcf,'color','white')
subplot(221)
plot(rate(:,2)*1000,sum(rate(:,3:2:211),2))
hold on
xlabel('Time (ms)')
ylabel('Inc. rate (m$^{-3}\cdot$s$^{-1}$)')
addlegend(leg_vals);
saveas(['inc_rate' studyid])

subplot(222)
% figure(3)
% set(gcf,'color','white')
plot(rate(:,2)*1000,rate(:,215))
hold on
xlabel('Time (ms)')
ylabel('Surf. growth rate (m$^{-3}\cdot$s$^{-1}$)')
addlegend(leg_vals);
saveas(['sg_rate' studyid])

subplot(223)
% figure(3)
% set(gcf,'color','white')
plot(rate(:,2)*1000,sum(rate(:,217:2:end),2))
hold on
xlabel('Time (ms)')
ylabel('Coag. rate (m$^{-3}\cdot$s$^{-1}$)')
addlegend(leg_vals);
saveas(['coag_rate' studyid])

subplot(224)
% figure(2)
% set(gcf,'color','white')
plot(rate(:,2)*1000,rate(:,213))
hold on
xlabel('Time (ms)')
ylabel('No inc. rate (m$^{-3}\cdot$s$^{-1}$)')
addlegend(leg_vals); 
saveas(['no_inc_rate' studyid])