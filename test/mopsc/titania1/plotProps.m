clc,
% clear, close all

%% Setup
% set defaults
set(0,'defaulttextinterpreter','tex')
set(0,'defaultaxesfontname','Helvetica')
set(0,'defaulttextfontname','Helvetica')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

% for finding files
basedir  = 'C:\Users\aab64\Documents\Data\MF_Comparison\P-no-dp-longer\';
filedir  = '';
filecmp  = '';
filebase = 'Network(stage1)';

% for plotting psl data
psltime  = '0.01';
npslbins = 10;

% for saving images
studyid  = '';
projpath = 'C:\Users\Astrid\Documents\Projects\';
projdir  = 'Network temperature dependence\Recent_work\';
imagedir = 'figures\AIWSWA\batch_comp\';
savefigs = 0;
leg_vals = '';
splots = 0;

part = csvread([basedir filedir filebase '-part.csv'],1);
pnpt = csvread([basedir filedir filebase '-part-temp.csv'],1);
rate = csvread([basedir filedir filebase '-part-rates.csv'],1);
% pslf = csvread([basedir filedir filebase '-psl(' psltime 's).csv'],1);
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

%% Part

figure(1000)
set(gcf,'color','white')
plot(part(:,2)*1000,pnpt(:,3),':')
hold on
xlabel('Time (ms)')
ylabel('Number of PN model particles (-)')
saveas(['nsp_pn_' studyid])

figure(1001)
set(gcf,'color','white')
plot(part(:,2)*1000,part(:,3),':')
hold on
xlabel('Time (ms)')
ylabel('Number of P model particles (-)')
saveas(['nsp_p' studyid])

figure(1002)
set(gcf,'color','white')
plot(part(:,2)*1000,part(:,3)+pnpt(:,3),':')
hold on
xlabel('Time (ms)')
ylabel('Number of particles (-)')
saveas(['nsp_p' studyid])

figure(1)
if splots
    subplot(231)
end
set(gcf,'color','white')
plot(part(:,2)*1000,(part(:,3)+pnpt(:,3))./part(:,5))
hold on
% plot(part(:,2)*1000,pnpt(:,3),':')
xlabel('Time (ms)')
ylabel('Sample volume (m^{-3})')
% legend('Ensemble','PN list')
saveas(['nsp_' studyid])


if ~splots
    figure(2)
    set(gcf,'color','white')
else
    subplot(232)
end
semilogy(part(:,2)*1000,part(:,5))
hold on
if (sum(part(:,6)) ~= 0)
    plot(part(:,2)*1000,part(:,5)+part(:,6),'k--')
    plot(part(:,2)*1000,part(:,5)-part(:,6),'k--')
end
xlabel('Time (ms)')
ylabel('M0 (m^{-3})')
saveas(['m0' studyid])

if ~splots
    figure(3)
    set(gcf,'color','white')
else
    subplot(233)
end
plot(part(:,2)*1000,part(:,9)*1e9)
hold on
if (sum(part(:,10)) ~= 0)
    plot(part(:,2)*1000,(part(:,9)+part(:,10))*1e9,'k--')
    plot(part(:,2)*1000,(part(:,9)-part(:,10))*1e9,'k--')
end
xlabel('Time (ms)')
ylabel('Collision diameter (nm)')
saveas(['dcol_ave' studyid])

if ~splots
    figure(4)
    set(gcf,'color','white')
else
    subplot(234)
end
plot(part(:,2)*1000,part(:,21))
hold on
if (sum(part(:,22)) ~= 0)
    plot(part(:,2)*1000,part(:,21)+part(:,22),'k--')
    plot(part(:,2)*1000,part(:,21)-part(:,22),'k--')
end
xlabel('Time (ms)')
ylabel('Mass (kg/m^3)')
saveas(['mass_ave' studyid])

if ~splots
    figure(5)
    set(gcf,'color','white')
else
    subplot(235)
end
plot(part(:,2)*1000,part(:,37))
hold on
if (sum(part(:,38)) ~= 0)
    plot(part(:,2)*1000,part(:,37)+part(:,38),'k--')
    plot(part(:,2)*1000,part(:,37)-part(:,38),'k--')
end
xlabel('Time (ms)')
ylabel('Number primaries per particle (-)')
saveas(['npri_ave' studyid])

if ~splots
    figure(6)
    set(gcf,'color','white')
else
    subplot(236)
end
plot(part(:,2)*1000,part(:,39)*1e9)
hold on
if (sum(part(:,40)) ~= 0)
    plot(part(:,2)*1000,(part(:,39)+part(:,40))*1e9,'k--')
    plot(part(:,2)*1000,(part(:,39)-part(:,40))*1e9,'k--')
end
xlabel('Time (ms)')
ylabel('Primary average diameter (nm)')
saveas(['dpri_ave' studyid])

figure(7)
set(gcf,'color','white')
subplot(131)
plot(part(:,2)*1000,part(:,33))
hold on
xlabel('Time (ms)')
ylabel('Total rutile (m^{-3})')
subplot(132)
plot(part(:,2)*1000,part(:,35))
hold on
xlabel('Time (ms)')
ylabel('Average rutile')
saveas(['rut_ave' studyid])

%% PSL

% figure(8)
% set(gcf,'color','white')
% [n_d,x_d] = histwc(pslf(:,3),pslf(:,9),npslbins);
% if ~splots
% else
%     subplot(231)
% end
% plot(x_d,n_d/max(n_d))
% hold on
% xlabel('Collision diameter (nm)')
% ylabel('Count divided by max. count (-)')
% ax = gca;
% ax.XLim(1) = 0.49;
% saveas(['dcol_psd' studyid])
% 
% [n_d,x_d] = histwc(pslf(:,5),pslf(:,9),npslbins);
% if ~splots
%     figure(9)
%     set(gcf,'color','white')
% else
%     subplot(232)
% end
% plot(x_d,n_d/max(n_d))
% hold on
% xlabel('Surface area (cm^{2})')
% ylabel('Count divided by max. count (-)')
% ax = gca;
% ax.XLim(1) = 0.0;
% saveas(['sa_psd' studyid])
% 
% [n_d,x_d] = histwc(pslf(:,12),pslf(:,9),npslbins);
% if ~splots
%     figure(10)
%     set(gcf,'color','white')
% else
%     subplot(233)
% end
% plot(x_d,n_d/max(n_d))
% hold on
% xlabel('Number of primaries (nm)')
% ylabel('Count divided by max. count (-)')
% ax = gca;
% ax.XLim(1) = 0.0;
% saveas(['npri_psd' studyid])
% 
% [n_d,x_d] = histwc(pslf(:,13),pslf(:,9),npslbins);
% if ~splots
%     figure(11)
%     set(gcf,'color','white')
% else
%     subplot(234)
% end
% plot(x_d,n_d/max(n_d))
% hold on
% xlabel('Primary diameter (nm)')
% ylabel('Count divided by max. count (-)')
% ax = gca;
% ax.XLim(1) = 0.49;
% saveas(['dpri_psd' studyid])
% 
% [n_d,x_d] = histwc(pslf(:,16),pslf(:,9),npslbins);
% if ~splots
%     figure(12)
%     set(gcf,'color','white')
% else
%     subplot(235)
% end
% plot(x_d,n_d/max(n_d))
% hold on
% xlabel('Std of primary diameter (nm)')
% ylabel('Count divided by max. count (-)')
% ax = gca;
% ax.XLim(1) = 0.0;
% saveas(['std_dpri_psd' studyid])
% 
% [n_d,x_d] = histwc(pslf(:,8),pslf(:,9),npslbins);
% if ~splots
%     figure(13)
%     set(gcf,'color','white')
% else
%     subplot(236)
% end
% plot(x_d,n_d/max(n_d))
% hold on
% xlabel('Particle age (s)')
% ylabel('Count divided by max. count (-)')
% ax = gca;
% ax.XLim(1) = 0.0;
% saveas(['page_psd' studyid])
% 
% figure(7)
% set(gcf,'color','white')
% subplot(133)
% [n_d,x_d] = histwc(pslf(:,11),pslf(:,9),npslbins);
% plot(x_d,n_d/max(n_d))
% hold on
% xlabel('Rutile units (-)')
% ylabel('Count divided by max. count (-)')
% ax = gca;
% ax.XLim(1) = 0.0;
% saveas(['rutile_psd' studyid])

%% Chem

figure(14)
set(gcf,'color','white')
if ~splots
else
    subplot(231)
end
plot(chem(:,2)*1000,chem(:,3)*1e6)
hold on
xlabel('Time (ms)')
ylabel('TiCl_4 conc. (mol/m^3)')
saveas(['ticl4' studyid])

if ~splots
    figure(15)
    set(gcf,'color','white')
else
    subplot(232)
end
plot(chem(:,2)*1000,chem(:,39)*1e6)
hold on
xlabel('Time (ms)')
ylabel('O_2 conc. (mol/m^3)')
saveas(['o2' studyid])

if ~splots
    figure(16)
    set(gcf,'color','white')
else
    subplot(233)
end
plot(chem(:,2)*1000,chem(:,45)*1e6)
hold on
xlabel('Time (ms)')
ylabel('Cl_2 conc. (mol/m^3)')
saveas(['cl2' studyid])

if ~splots
    figure(17)
    set(gcf,'color','white')
else
    subplot(234)
end
plot(chem(:,2)*1000,chem(:,19)*1e6)
hold on
xlabel('Time (ms)')
ylabel('TiO_2Cl_3 conc. (mol/m^3)')
saveas(['tio2cl3' studyid])

if ~splots
    figure(18)
    set(gcf,'color','white')
else
    subplot(235)
end
plot(chem(:,2)*1000,chem(:,57)*1e6)
hold on
xlabel('Time (ms)')
ylabel('Ar conc. (mol/m^3)')
saveas(['ar' studyid])

if ~splots
    figure(19)
    set(gcf,'color','white')
else
    subplot(236)
end
plot(chem(:,2)*1000,chem(:,61))
hold on
xlabel('Time (ms)')
ylabel('Temperature (K)')
saveas(['temp' studyid])

%% CPUT

cadj = [0;(cput(2:end,3)-cput(2,3))/60];
figure(20)
set(gcf,'color','white')
plot(cput(:,2)*1000,cadj)
hold on
xlabel('Time (ms)')
ylabel('Solver time (min)')
saveas(['cput' studyid])

%% Rates

figure(21)
set(gcf,'color','white')
subplot(131)
semilogy(rate(:,2)*1000,sum(rate(:,3:2:212),2))
hold on
xlabel('Time (ms)')
ylabel('Inc. rate (1/(m^3s))')
saveas(['inc_rate' studyid])

% figure(22)
% set(gcf,'color','white')
subplot(132)
semilogy(rate(:,2)*1000,rate(:,215))
hold on
xlabel('Time (ms)')
ylabel('Surf. growth rate (m^2/(m^3s))')
saveas(['sg_rate' studyid])

% figure(23)
% set(gcf,'color','white')
subplot(133)
semilogy(rate(:,2)*1000,sum(rate(:,217:2:end),2))
hold on
xlabel('Time (ms)')
ylabel('Coag. rate (1/(m^3s))')
saveas(['coag_rate' studyid])

xs = chem(end,3:2:59)./sum(chem(end,3:2:59));
ind=find(xs>=1e-5);
disp(ind)
disp(xs(ind))
disp(part(end,5))