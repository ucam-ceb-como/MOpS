clc, 
clear, close all

%% Setup

set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

nsbins = 25;
sigma  = 0.1;
npts   = 1000;
leg    = {'SWA1','SWA4','SWA9','SWA9b','DSA'};
swa_leg = {'SWA1','SWA4','SWA9','SWA9b'};

% for saving images
studyid  = '';
projpath = 'C:\Users\Astrid\Documents\Projects\';
projdir  = 'Network temperature dependence\Recent_work\';
imagedir = 'figures\swa_comparisons\real\largerSV_cdelete\'; 
savefigs = 0;
% initialise saveas function if saving figures, otherwise do nothing here
if (savefigs)
    disp('Saving figures - continuing may overwrite existing!')
    pause 
    disp('Continuing...')
    saveas = @(str)(print([projpath projdir imagedir str],'-depsc'));
else
    saveas = @(str)(disp('Not saving'));
end

%% Data
base0 = 'vienna/real/comparisons/13ms/cdelete/largerSV/5runs/';
base = 'darwin/real/comparisons/13ms/';

s_swa9 = csvread([base 'swa9b/Network(stage1)-psl(0.013s).csv'],1);
s_swa8 = csvread([base 'swa9/Network(stage1)-psl(0.013s).csv'],1);
s_swa4 = csvread([base 'swa4/Network(stage1)-psl(0.013s).csv'],1);
s_swa1 = csvread([base 'swa1/Network(stage1)-psl(0.013s).csv'],1);
s_dsa  = csvread([base 'dsa/Network(stage1)-psl(0.013s).csv'],1);

p_swa9 = csvread([base 'swa9b/Network(stage1)-part.csv'],1);
p_swa8 = csvread([base 'swa9/Network(stage1)-part.csv'],1);
p_swa4 = csvread([base 'swa4/Network(stage1)-part.csv'],1);
p_swa1 = csvread([base 'swa1/Network(stage1)-part.csv'],1);
p_dsa  = csvread([base 'dsa/Network(stage1)-part.csv'],1);

c_swa9 = csvread([base 'swa9b/Network(stage1)-cput.csv'],1);
c_swa8 = csvread([base 'swa9/Network(stage1)-cput.csv'],1);
c_swa4 = csvread([base 'swa4/Network(stage1)-cput.csv'],1);
c_swa1 = csvread([base 'swa1/Network(stage1)-cput.csv'],1);
c_dsa  = csvread([base 'dsa/Network(stage1)-cput.csv'],1);

%% PSL - plot basic distributions

% [n_d,x_d]   = histwc(s_dsa(:,3),s_dsa(:,9),nsbins);
% [n_s1,x_s1] = histwc(s_swa1(:,3),s_swa1(:,9),nsbins);
% 
% figure(1)
% set(gcf,'color','white')
% plot(x_d,n_d/max(n_d),x_s1,n_s1/max(n_s1),'--')
% hold on
% xlabel('Collision diameter (nm)')
% ylabel('Count divided by max. count (-)')
% l = legend(leg);
% l.Location = 'NorthEast';
% l.Interpreter = 'latex';
% ax = gca;
% ax.XLim(1) = 0.49;
% 
% [n_d,x_d]   = histwc(s_dsa(:,12),s_dsa(:,9),nsbins);
% [n_s1,x_s1] = histwc(s_swa1(:,12),s_swa1(:,9),nsbins);
% 
% figure(2)
% set(gcf,'color','white')
% plot(x_d,n_d/max(n_d),x_s1,n_s1/max(n_s1),'--')
% hold on
% xlabel('Primary diameter (nm)')
% ylabel('Count divided by max. count (-)')
% l = legend(leg);
% l.Location = 'NorthEast';
% l.Interpreter = 'latex';
% ax = gca;
% ax.XLim(1) = 0.49;

%% dc - plot lognormal distributions and compute KL divergences from DSA

% Compute lognormal PDF
[D_s1c,DgofD_s1c] = kernelest_w(p_swa1(end,5),s_swa1(:,3),sigma,npts,[],s_swa1(:,9));
[D_s9c,DgofD_s9c] = kernelest_w(p_swa9(end,5),s_swa9(:,3),sigma,npts,[],s_swa9(:,9));
[D_s8c,DgofD_s8c] = kernelest_w(p_swa8(end,5),s_swa8(:,3),sigma,npts,[],s_swa8(:,9));
[D_s4c,DgofD_s4c] = kernelest_w(p_swa4(end,5),s_swa4(:,3),sigma,npts,[],s_swa4(:,9));
[D_dc,DgofD_dc]   = kernelest_w(p_dsa(end,5),s_dsa(:,3),sigma,npts,[],s_dsa(:,9));

figure(3)
set(gcf,'color','white')
semilogx(D_s1c,DgofD_s1c,'-',D_s4c,DgofD_s4c,'--',D_s8c,DgofD_s8c,'-.',D_s9c,DgofD_s9c,':',D_dc,DgofD_dc,'-')
hold on
xlabel('Collision diameter, $d_{\textrm{c}}$ (nm)')
ylabel('$\textrm{d}n/\textrm{dln}(d_{\textrm{c}})$ (m$^{-3}$)')
% set(gca,'YLim',[1 10^17])
ax = gca;
ax.XLim(1) = 0.49;
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
h = gca;
h.Children(1).LineWidth = 3.5;
h.Children(2).LineWidth = 2;
h.Children(3).LineWidth = 2;
h.Children(4).LineWidth = 2;
h.Children(5).LineWidth = 3.5;
h.Children(1).Color = 'k';
h.Children(2).Color = [0      0.7490 0.7490];
h.Children(3).Color = [0.0784 0.1686 0.5490];
h.Children(4).Color = [0.3020 0.7451 0.9333];
h.Children(5).Color = 'r';
% saveas('dcol_5dist_real1')

disp('Collision diameter distributions, D_KL swa 1, 4, 8, 9')
disp(klest(s_dsa(:,3),s_swa1(:,3),s_dsa(:,9),s_swa1(:,9)))
disp(klest(s_dsa(:,3),s_swa4(:,3),s_dsa(:,9),s_swa4(:,9)))
disp(klest(s_dsa(:,3),s_swa8(:,3),s_dsa(:,9),s_swa8(:,9)))
disp(klest(s_dsa(:,3),s_swa9(:,3),s_dsa(:,9),s_swa9(:,9)))

%% dp - plot lognormal distributions and compute KL divergences from DSA

% Compute lognormal PDF
[D_s1p,DgofD_s1p] = kernelest_w(p_swa1(end,5),s_swa1(:,13),sigma,npts,[],s_swa1(:,9));
[D_s9p,DgofD_s9p] = kernelest_w(p_swa9(end,5),s_swa9(:,13),sigma,npts,[],s_swa9(:,9));
[D_s8p,DgofD_s8p] = kernelest_w(p_swa8(end,5),s_swa8(:,13),sigma,npts,[],s_swa8(:,9));
[D_s4p,DgofD_s4p] = kernelest_w(p_swa4(end,5),s_swa4(:,13),sigma,npts,[],s_swa4(:,9));
[D_dp,DgofD_dp]   = kernelest_w(p_dsa(end,5),s_dsa(:,13),sigma,npts,[],s_dsa(:,9));

figure(4)
set(gcf,'color','white')
semilogx(D_s1p,DgofD_s1p,'-',D_s4p,DgofD_s4p,'--',D_s8p,DgofD_s8p,'-.',D_s9p,DgofD_s9p,':',D_dp,DgofD_dp,'-')
hold on
xlabel('Primary diameter, $d_{\textrm{p}}$ (nm)')
ylabel('$\textrm{d}n/\textrm{dln}(d_{\textrm{p}})$ (m$^{-3}$)')
% set(gca,'YLim',[1 10^17])
ax = gca;
ax.XLim(1) = 0.49;
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
h = gca;
h.Children(1).LineWidth = 3.5;
h.Children(2).LineWidth = 2;
h.Children(3).LineWidth = 2;
h.Children(4).LineWidth = 2;
h.Children(5).LineWidth = 3.5;
h.Children(1).Color = 'k';
h.Children(2).Color = [0      0.7490 0.7490];
h.Children(3).Color = [0.0784 0.1686 0.5490];
h.Children(4).Color = [0.3020 0.7451 0.9333];
h.Children(5).Color = 'r';
% saveas('dpri_5dist_real1')

disp('Primary diameter distributions, D_KL swa 1, 4, 8, 9')
disp(klest(s_dsa(:,13),s_swa1(:,13),s_dsa(:,9),s_swa1(:,9)))
disp(klest(s_dsa(:,13),s_swa4(:,13),s_dsa(:,9),s_swa4(:,9)))
disp(klest(s_dsa(:,13),s_swa8(:,13),s_dsa(:,9),s_swa8(:,9)))
disp(klest(s_dsa(:,13),s_swa9(:,13),s_dsa(:,9),s_swa9(:,9)))

disp('Number of primaries distributions, D_KL swa 1, 4, 8, 9')
disp(klest(s_dsa(:,12),s_swa1(:,12),s_dsa(:,9),s_swa1(:,9)))
disp(klest(s_dsa(:,12),s_swa4(:,12),s_dsa(:,9),s_swa4(:,9)))
disp(klest(s_dsa(:,12),s_swa8(:,12),s_dsa(:,9),s_swa8(:,9)))
disp(klest(s_dsa(:,12),s_swa9(:,12),s_dsa(:,9),s_swa9(:,9)))

%% CPUT - plot average final solver time

figure(5)
set(gcf,'color','white')
bar(1,(c_swa1(end,3)-c_swa1(2,3))/3600,'FaceColor','b')
hold on 
bar(2,(c_swa4(end,3)-c_swa4(2,3))/3600,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,(c_swa8(end,3)-c_swa8(2,3))/3600,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,(c_swa9(end,3)-c_swa9(2,3))/3600,'FaceColor',[0      0.7490 0.7490])
bar(5,(c_dsa(end,3)-c_dsa(2,3))/3600,'FaceColor','r')
set(gca,'XTick',1:5,'XTickLabel',leg)
ylabel('Final solver time (hr)')
saveas(['stime_5dist_real1' studyid])

%% REs

rem_s1 = abs(p_dsa(end,5)-p_swa1(end,5))/p_dsa(end,5);
rem_s4 = abs(p_dsa(end,5)-p_swa4(end,5))/p_dsa(end,5);
rem_s8 = abs(p_dsa(end,5)-p_swa8(end,5))/p_dsa(end,5);
rem_s9 = abs(p_dsa(end,5)-p_swa9(end,5))/p_dsa(end,5);

rec_s1 = abs(p_dsa(end,9)-p_swa1(end,9))/p_dsa(end,9);
rec_s4 = abs(p_dsa(end,9)-p_swa4(end,9))/p_dsa(end,9);
rec_s8 = abs(p_dsa(end,9)-p_swa8(end,9))/p_dsa(end,9);
rec_s9 = abs(p_dsa(end,9)-p_swa9(end,9))/p_dsa(end,9);

ren_s1 = abs(p_dsa(end,37)-p_swa1(end,37))/p_dsa(end,37);
ren_s4 = abs(p_dsa(end,37)-p_swa4(end,37))/p_dsa(end,37);
ren_s8 = abs(p_dsa(end,37)-p_swa8(end,37))/p_dsa(end,37);
ren_s9 = abs(p_dsa(end,37)-p_swa9(end,37))/p_dsa(end,37);

rep_s1 = abs(p_dsa(end,39)-p_swa1(end,39))/p_dsa(end,39);
rep_s4 = abs(p_dsa(end,39)-p_swa4(end,39))/p_dsa(end,39);
rep_s8 = abs(p_dsa(end,39)-p_swa8(end,39))/p_dsa(end,39);
rep_s9 = abs(p_dsa(end,39)-p_swa9(end,39))/p_dsa(end,39);

rek_s1 = abs(p_dsa(end,21)-p_swa1(end,21))/p_dsa(end,21);
rek_s4 = abs(p_dsa(end,21)-p_swa4(end,21))/p_dsa(end,21);
rek_s8 = abs(p_dsa(end,21)-p_swa8(end,21))/p_dsa(end,21);
rek_s9 = abs(p_dsa(end,21)-p_swa9(end,21))/p_dsa(end,21);

figure(6)
set(gcf,'color','white')
bar(1,rem_s1,'FaceColor','r')
hold on 
bar(2,rem_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,rem_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,rem_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',swa_leg)
ylabel('Absolute relative error in $M_0$ (-)')
saveas('rem_5dist')

figure(7)
set(gcf,'color','white')
bar(1,rec_s1,'FaceColor','r')
hold on 
bar(2,rec_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,rec_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,rec_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',swa_leg)
ylabel('Absolute relative error in $d_{\textrm{c}}$ (-)')
saveas('rec_5dist')

figure(8)
set(gcf,'color','white')
bar(1,ren_s1,'FaceColor','r')
hold on 
bar(2,ren_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,ren_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,ren_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',swa_leg)
ylabel('Absolute relative error in $n_{\textrm{pri}}$ (-)')
saveas('ren_5dist')

figure(9)
set(gcf,'color','white')
bar(1,rep_s1,'FaceColor','r')
hold on 
bar(2,rep_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,rep_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,rep_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',swa_leg)
ylabel('Absolute relative error in $d_{\textrm{p}}$ (-)')
saveas('rep_5dist')

figure(10)
set(gcf,'color','white')
bar(1,rek_s1,'FaceColor','r')
hold on 
bar(2,rek_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,rek_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,rek_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',swa_leg)
ylabel('Absolute relative error in $M_1$ (-)')
saveas('rek_5dist')

%% RUs

rum_d  = p_dsa(end,6)/p_dsa(end,5);
rum_s1 = p_swa1(end,6)/p_swa1(end,5);
rum_s4 = p_swa4(end,6)/p_swa4(end,5);
rum_s8 = p_swa8(end,6)/p_swa8(end,5);
rum_s9 = p_swa9(end,6)/p_swa9(end,5);

ruc_d  = p_dsa(end,10)/p_dsa(end,9);
ruc_s1 = p_swa1(end,10)/p_swa1(end,9);
ruc_s4 = p_swa4(end,10)/p_swa4(end,9);
ruc_s8 = p_swa8(end,10)/p_swa8(end,9);
ruc_s9 = p_swa9(end,10)/p_swa9(end,9);

run_d  = p_dsa(end,38)/p_dsa(end,37);
run_s1 = p_swa1(end,38)/p_swa1(end,37);
run_s4 = p_swa4(end,38)/p_swa4(end,37);
run_s8 = p_swa8(end,38)/p_swa8(end,37);
run_s9 = p_swa9(end,38)/p_swa9(end,37);

rup_d  = p_dsa(end,40)/p_dsa(end,39);
rup_s1 = p_swa1(end,40)/p_swa1(end,39);
rup_s4 = p_swa4(end,40)/p_swa4(end,39);
rup_s8 = p_swa8(end,40)/p_swa8(end,39);
rup_s9 = p_swa9(end,40)/p_swa9(end,39);

ruk_d  = p_dsa(end,22)/p_dsa(end,21);
ruk_s1 = p_swa1(end,22)/p_swa1(end,21);
ruk_s4 = p_swa4(end,22)/p_swa4(end,21);
ruk_s8 = p_swa8(end,22)/p_swa8(end,21);
ruk_s9 = p_swa9(end,22)/p_swa9(end,21);

% figure(11)
% set(gcf,'color','white')
% % bar(1,rum_d)
% % hold on 
% bar(1,rum_s1,'FaceColor','b')
% hold on
% bar(2,rum_s4,'FaceColor',[0.3020 0.7451 0.9333])
% bar(3,rum_s8,'FaceColor',[0.0784 0.1686 0.5490])
% bar(4,rum_s9,'FaceColor',[0      0.7490 0.7490])
% set(gca,'XTick',1:4,'XTickLabel',swa_leg)
% ylabel('Relative uncertainty in $M_0$ (-)')
% saveas('rum_5dist')
% 
% figure(12)
% set(gcf,'color','white')
% % bar(1,ruc_d)
% % hold on 
% bar(1,ruc_s1,'FaceColor','b')
% hold on
% bar(2,ruc_s4,'FaceColor',[0.3020 0.7451 0.9333])
% bar(3,ruc_s8,'FaceColor',[0.0784 0.1686 0.5490])
% bar(4,ruc_s9,'FaceColor',[0      0.7490 0.7490])
% set(gca,'XTick',1:4,'XTickLabel',swa_leg)
% ylabel('Relative uncertainty in $d_{\textrm{c}}$ (-)')
% saveas('ruc_5dist')
% 
% figure(13)
% set(gcf,'color','white')
% % bar(1,run_d)
% % hold on 
% bar(1,run_s1,'FaceColor','b')
% hold on
% bar(2,run_s4,'FaceColor',[0.3020 0.7451 0.9333])
% bar(3,run_s8,'FaceColor',[0.0784 0.1686 0.5490])
% bar(4,run_s9,'FaceColor',[0      0.7490 0.7490])
% set(gca,'XTick',1:4,'XTickLabel',swa_leg)
% ylabel('Relative uncertainty in $n_{\textrm{pri}}$ (-)')
% saveas('run_5dist')
% 
% figure(14)
% set(gcf,'color','white')
% % bar(1,rup_d)
% % hold on 
% bar(1,rup_s1,'FaceColor','b')
% hold on
% bar(2,rup_s4,'FaceColor',[0.3020 0.7451 0.9333])
% bar(3,rup_s8,'FaceColor',[0.0784 0.1686 0.5490])
% bar(4,rup_s9,'FaceColor',[0      0.7490 0.7490])
% set(gca,'XTick',1:4,'XTickLabel',swa_leg)
% ylabel('Relative uncertainty in $d_{\textrm{p}}$ (-)')
% saveas('rup_5dist')
% 
% figure(15)
% set(gcf,'color','white')
% % bar(1,ruk_d)
% % hold on 
% bar(1,ruk_s1,'FaceColor','b')
% hold on
% bar(2,ruk_s4,'FaceColor',[0.3020 0.7451 0.9333])
% bar(3,ruk_s8,'FaceColor',[0.0784 0.1686 0.5490])
% bar(4,ruk_s9,'FaceColor',[0      0.7490 0.7490])
% set(gca,'XTick',1:4,'XTickLabel',swa_leg)
% ylabel('Relative uncertainty in $M_1$ (-)')
% saveas('ruk_5dist')

%% Display relative errors compared to DSA and relative uncertainties

disp('FINAL RECKONING')
disp('---------------')
disp('* m0 re swa 1, 4, 8, 9')
disp(rem_s1)
disp(rem_s4)
disp(rem_s8)
disp(rem_s9)
disp('* m0 ru dsa, m0 ru swa 1, 4, 8, 9')
disp(rum_d)
disp(rum_s1)
disp(rum_s4)
disp(rum_s8)
disp(rum_s9)

disp('* dc re swa 1, 4, 8, 9')
disp(rec_s1)
disp(rec_s4)
disp(rec_s8)
disp(rec_s9)
disp('* dc ru dsa, dc ru swa 1, 4, 8, 9')
disp(ruc_d)
disp(ruc_s1)
disp(ruc_s4)
disp(ruc_s8)
disp(ruc_s9)

disp('* npri re swa 1, 4, 8, 9')
disp(ren_s1)
disp(ren_s4)
disp(ren_s8)
disp(ren_s9)
disp('* npri ru dsa, npri ru swa 1, 4, 8, 9')
disp(run_d)
disp(run_s1)
disp(run_s4)
disp(run_s8)
disp(run_s9)

disp('* dp re swa 1, 4, 8, 9')
disp(rep_s1)
disp(rep_s4)
disp(rep_s8)
disp(rep_s9)
disp('* dp ru dsa, dp ru swa 1, 4, 8, 9')
disp(rup_d)
disp(rup_s1)
disp(rup_s4)
disp(rup_s8)
disp(rup_s9)

disp('* m1 re swa 1, 4, 8, 9')
disp(rek_s1)
disp(rek_s4)
disp(rek_s8)
disp(rek_s9)
disp('* m1 ru dsa, m1 ru swa 1, 4, 8, 9')
disp(ruk_d)
disp(ruk_s1)
disp(ruk_s4)
disp(ruk_s8)
disp(ruk_s9)

%% REs to SWA

leg    = {'SWA','AIW$_{\textrm{e}}$','PSC$_{\textrm{b}}$','PSC$_{\textrm{w}}$'};
leg2   = {'AIW$_{\textrm{e}}$','PSC$_{\textrm{b}}$','PSC$_{\textrm{w}}$'};

% leg    = {'SWA','AIW$_{\textrm{e}}$','HCI','PSC$_{\textrm{w}}$'};
% leg2   = {'AIW$_{\textrm{e}}$','HCI','PSC$_{\textrm{w}}$'};

rem2_s1 = abs(p_swa1(end,5)-p_swa1(end,5))/p_swa1(end,5);
rem2_s4 = abs(p_swa1(end,5)-p_swa4(end,5))/p_swa1(end,5);
rem2_s8 = abs(p_swa1(end,5)-p_swa8(end,5))/p_swa1(end,5);
rem2_s9 = abs(p_swa1(end,5)-p_swa9(end,5))/p_swa1(end,5);

rec2_s1 = abs(p_swa1(end,9)-p_swa1(end,9))/p_swa1(end,9);
rec2_s4 = abs(p_swa1(end,9)-p_swa4(end,9))/p_swa1(end,9);
rec2_s8 = abs(p_swa1(end,9)-p_swa8(end,9))/p_swa1(end,9);
rec2_s9 = abs(p_swa1(end,9)-p_swa9(end,9))/p_swa1(end,9);

ren2_s1 = abs(p_swa1(end,37)-p_swa1(end,37))/p_swa1(end,37);
ren2_s4 = abs(p_swa1(end,37)-p_swa4(end,37))/p_swa1(end,37);
ren2_s8 = abs(p_swa1(end,37)-p_swa8(end,37))/p_swa1(end,37);
ren2_s9 = abs(p_swa1(end,37)-p_swa9(end,37))/p_swa1(end,37);

rep2_s1 = abs(p_swa1(end,39)-p_swa1(end,39))/p_swa1(end,39);
rep2_s4 = abs(p_swa1(end,39)-p_swa4(end,39))/p_swa1(end,39);
rep2_s8 = abs(p_swa1(end,39)-p_swa8(end,39))/p_swa1(end,39);
rep2_s9 = abs(p_swa1(end,39)-p_swa9(end,39))/p_swa1(end,39);

rek2_s1 = abs(p_swa1(end,21)-p_swa1(end,21))/p_swa1(end,21);
rek2_s4 = abs(p_swa1(end,21)-p_swa4(end,21))/p_swa1(end,21);
rek2_s8 = abs(p_swa1(end,21)-p_swa8(end,21))/p_swa1(end,21);
rek2_s9 = abs(p_swa1(end,21)-p_swa9(end,21))/p_swa1(end,21);

figure(16)
set(gcf,'color','white')
bar(1,rem2_s4,'FaceColor',[0.3020 0.7451 0.9333])
hold on
bar(2,rem2_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(3,rem2_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:3,'XTickLabel',leg2,'TickLabelInterpreter','latex')
ylabel('Absolute relative error in $M_0$ (-)')
saveas(['rem2_5dist' studyid])

figure(17)
set(gcf,'color','white')
bar(1,rec2_s4,'FaceColor',[0.3020 0.7451 0.9333])
hold on
bar(2,rec2_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(3,rec2_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:3,'XTickLabel',leg2,'TickLabelInterpreter','latex')
ylabel('Absolute relative error in $d_{\textrm{c}}$ (-)')
saveas(['rec2_5dist' studyid])

figure(18)
set(gcf,'color','white')
bar(1,ren2_s4,'FaceColor',[0.3020 0.7451 0.9333])
hold on
bar(2,ren2_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(3,ren2_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:3,'XTickLabel',leg2,'TickLabelInterpreter','latex')
ylabel('Absolute relative error in $n_{\textrm{pri}}$ (-)')
saveas(['ren2_5dist' studyid])

figure(19)
set(gcf,'color','white')
bar(1,rep2_s4,'FaceColor',[0.3020 0.7451 0.9333])
hold on
bar(2,rep2_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(3,rep2_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:3,'XTickLabel',leg2,'TickLabelInterpreter','latex')
ylabel('Absolute relative error in $d_{\textrm{p}}$ (-)')
saveas(['rep2_5dist' studyid])

figure(20)
set(gcf,'color','white')
bar(1,rek2_s4,'FaceColor',[0.3020 0.7451 0.9333])
hold on
bar(2,rek2_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(3,rek2_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:3,'XTickLabel',leg2,'TickLabelInterpreter','latex')
ylabel('Absolute relative error in $M_1$ (-)')
saveas(['rek2_5dist' studyid])

%% RUs 2

figure(21)
set(gcf,'color','white')
bar(1,rum_s1,'FaceColor','r')
hold on
bar(2,rum_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,rum_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,rum_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',leg,'TickLabelInterpreter','latex')
ylabel('Relative uncertainty in $M_0$ (-)')
saveas(['rum2_5dist' studyid])

figure(22)
set(gcf,'color','white')
bar(1,ruc_s1,'FaceColor','r')
hold on
bar(2,ruc_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,ruc_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,ruc_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',leg,'TickLabelInterpreter','latex')
ylabel('Relative uncertainty in $d_{\textrm{c}}$ (-)')
saveas(['ruc2_5dist' studyid])

figure(23)
set(gcf,'color','white')
bar(1,run_s1,'FaceColor','r')
hold on
bar(2,run_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,run_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,run_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',leg,'TickLabelInterpreter','latex')
ylabel('Relative uncertainty in $n_{\textrm{pri}}$ (-)')
saveas(['run2_5dist' studyid])

figure(24)
set(gcf,'color','white')
bar(1,rup_s1,'FaceColor','r')
hold on
bar(2,rup_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,rup_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,rup_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',leg,'TickLabelInterpreter','latex')
ylabel('Relative uncertainty in $d_{\textrm{p}}$ (-)')
saveas(['rup2_5dist' studyid])

figure(25)
set(gcf,'color','white')
bar(1,ruk_s1,'FaceColor','r')
hold on
bar(2,ruk_s4,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,ruk_s8,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,ruk_s9,'FaceColor',[0      0.7490 0.7490])
set(gca,'XTick',1:4,'XTickLabel',leg,'TickLabelInterpreter','latex')
ylabel('Relative uncertainty in $M_1$ (-)')
saveas(['ruk2_5dist' studyid])


%% CPUT 2

cmean=(c_swa1(end,3)-c_swa1(2,3))/3600;

figure(26)
set(gcf,'color','white')
bar(1,(c_swa1(end,3)-c_swa1(2,3))/3600,'FaceColor','r')
hold on 
bar(2,(c_swa4(end,3)-c_swa4(2,3))/3600,'FaceColor',[0.3020 0.7451 0.9333])
bar(3,(c_swa8(end,3)-c_swa8(2,3))/3600,'FaceColor',[0.0784 0.1686 0.5490])
bar(4,(c_swa9(end,3)-c_swa9(2,3))/3600,'FaceColor',[0      0.7490 0.7490])
% errorbar(1,(c_swa1(end,3)-c_swa1(2,3))/3600,(c_swa1(end,4)-c_swa1(2,4))/3600,'Color','k')
% errorbar(2,(c_swa4(end,3)-c_swa4(2,3))/3600,(c_swa4(end,4)-c_swa4(2,4))/3600,'Color','k')
% errorbar(3,(c_swa8(end,3)-c_swa8(2,3))/3600,(c_swa8(end,4)-c_swa8(2,4))/3600,'Color','k')
% errorbar(4,(c_swa9(end,3)-c_swa9(2,3))/3600,(c_swa9(end,4)-c_swa9(2,4))/3600,'Color','k')
text(0.65,1.04*(c_swa1(end,3)-c_swa1(2,3))/3600,...
    [num2str(100*(c_swa1(end,3)-c_swa1(2,3))/3600/cmean,'%1.1f') '\%'],...
    'Interpreter','latex','FontSize',14)
text(1.7,1.04*(c_swa4(end,3)-c_swa4(2,3))/3600,...
    [num2str(100*(c_swa4(end,3)-c_swa4(2,3))/3600/cmean,'%1.1f') '\%'],...
    'Interpreter','latex','FontSize',14)
text(2.7,1.04*(c_swa8(end,3)-c_swa8(2,3))/3600,...
    [num2str(100*(c_swa8(end,3)-c_swa8(2,3))/3600/cmean,'%1.1f') '\%'],...
    'Interpreter','latex','FontSize',14)
text(3.7,1.04*(c_swa9(end,3)-c_swa9(2,3))/3600,...
    [num2str(100*(c_swa9(end,3)-c_swa9(2,3))/3600/cmean,'%1.1f') '\%'],...
    'Interpreter','latex','FontSize',14)
set(gca,'XTick',1:4,'XTickLabel',leg,'YTick',0:2:8,'YTickLabelRotation',90,...
    'TickLabelInterpreter','latex')
ylabel('Final solver time (hr)')
saveas(['stime2_5dist_real1' studyid])

%% NSP

pmean = mean(p_swa1(:,3));
     
figure(27)
set(gcf,'color','white')
bar(1,mean(p_swa1(:,3)),'FaceColor','r')
hold on 
bar(2,mean(p_swa4(:,3)),'FaceColor',[0.3020 0.7451 0.9333])
bar(3,mean(p_swa8(:,3)),'FaceColor',[0.0784 0.1686 0.5490])
bar(4,mean(p_swa9(:,3)),'FaceColor',[0      0.7490 0.7490])
% errorbar(1,mean(p_swa1(:,3)),mean(p_swa1(:,4)),'Color','k')
% errorbar(2,mean(p_swa4(:,3)),mean(p_swa4(:,4)),'Color','k')
% errorbar(3,mean(p_swa8(:,3)),mean(p_swa8(:,4)),'Color','k')
% errorbar(4,mean(p_swa9(:,3)),mean(p_swa9(:,4)),'Color','k')
text(0.65,1.04*mean(p_swa1(:,3)),...
    [num2str(100*mean(p_swa1(:,3))/pmean,'%1.1f') '\%'],...
    'Interpreter','latex','FontSize',14)
text(1.7,1.04*mean(p_swa4(:,3)),...
    [num2str(100*mean(p_swa4(:,3))/pmean,'%1.1f') '\%'],...
    'Interpreter','latex','FontSize',14)
text(2.7,1.04*mean(p_swa8(:,3)),...
    [num2str(100*mean(p_swa8(:,3))/pmean,'%1.1f') '\%'],...
    'Interpreter','latex','FontSize',14)
text(3.7,1.04*mean(p_swa9(:,3)),...
    [num2str(100*mean(p_swa9(:,3))/pmean,'%1.1f') '\%'],...
    'Interpreter','latex','FontSize',14)
set(gca,'XTick',1:4,'XTickLabel',leg,'YLim',[0 20000],'YTick',0:4000:16000,...
    'YTickLabel',{'0','4000','6000','8000','12000','16000'},...
    'YTickLabelRotation',90,'TickLabelInterpreter','latex')
ylabel('Average number of particles (-)')
saveas(['nsp2_5dist_real1' studyid])

figure(28)
set(gcf,'color','white')
loglog(p_swa1(:,2),p_swa1(:,3),'-','Color','r')
hold on 
loglog(p_swa1(:,2),p_swa4(:,3),'--','Color',[0.3020 0.7451 0.9333])
loglog(p_swa1(:,2),p_swa8(:,3),'-.','Color',[0.0784 0.1686 0.5490])
loglog(p_swa1(:,2),p_swa9(:,3),':','Color',[0      0.7490 0.7490])
ylabel('Number of particles (-)')
xlabel('Time')
set(gca,'TickLabelInterpreter','latex')
l = legend(leg);
l.Location = 'SouthEast';
l.Interpreter = 'latex';
saveas(['nspt_5dist_real1' studyid])


%% dc 2

% Convert SWA weights to equivalent number of particles
% dcnew_s1 = convert_weighted(s_swa1(:,3),s_swa1(:,9));
% dcnew_s9 = convert_weighted(s_swa9(:,3),s_swa9(:,9));
% dcnew_s8 = convert_weighted(s_swa8(:,3),s_swa8(:,9));
% dcnew_s4 = convert_weighted(s_swa4(:,3),s_swa4(:,9));

% Compute lognormal PDF (m0,dn,sigma,npts,D,w)
[D_s1c,DgofD_s1c] = kernelest_w(p_swa1(end,5),s_swa1(:,3),sigma,npts,[],s_swa1(:,9));
[D_s4c,DgofD_s4c] = kernelest_w(p_swa4(end,5),s_swa4(:,3),sigma,npts,[],s_swa4(:,9));
[D_s8c,DgofD_s8c] = kernelest_w(p_swa8(end,5),s_swa8(:,3),sigma,npts,[],s_swa8(:,9));
[D_s9c,DgofD_s9c] = kernelest_w(p_swa9(end,5),s_swa9(:,3),sigma,npts,[],s_swa9(:,9));
[D_dc,DgofD_dc]   = kernelest(p_dsa(end,5),s_dsa(:,3),sigma,npts,[]);

figure(29)
set(gcf,'color','white')
semilogx(D_s1c,DgofD_s1c,'-',D_s4c,DgofD_s4c,'--',D_s8c,DgofD_s8c,'-.',D_s9c,DgofD_s9c,':')
hold on
xlabel('Collision diameter, $d_{\textrm{c}}$ (nm)')
ylabel('$\textrm{d}n/\textrm{dln}(d_{\textrm{c}})$ (m$^{-3}$)')
ax = gca;
ax.XLim(1) = 0.49;
ax.XLim(2) = 1e4;
ax.TickLabelInterpreter = 'latex';
% ax.YLim=[0 7e16];
l = legend(leg);
l.Location = 'NorthWest';
l.Interpreter = 'latex';
h = gca;
h.Children(1).LineWidth = 2;
h.Children(2).LineWidth = 2;
h.Children(3).LineWidth = 2;
h.Children(4).LineWidth = 3.5;
h.Children(1).Color = [0      0.7490 0.7490];
h.Children(2).Color = [0.0784 0.1686 0.5490];
h.Children(3).Color = [0.3020 0.7451 0.9333];
h.Children(4).Color = 'r';
saveas(['dcol2_5dist_real1' studyid])

disp('Collision diameter distributions, D_KL swa 4, 8, 9')
disp(klest(s_swa1(:,3),s_swa4(:,3),s_swa1(:,9),s_swa4(:,9)))
disp(klest(s_swa1(:,3),s_swa8(:,3),s_swa1(:,9),s_swa8(:,9)))
disp(klest(s_swa1(:,3),s_swa9(:,3),s_swa1(:,9),s_swa9(:,9)))


%% npri - plot lognormal distributions and compute KL divergences from DSA

[D_s1pri,DgofD_s1pri] = kernelest_w(p_swa1(end,5),s_swa1(:,12),sigma,npts,[],s_swa1(:,9));
[D_s4pri,DgofD_s4pri] = kernelest_w(p_swa4(end,5),s_swa4(:,12),sigma,npts,[],s_swa4(:,9));
[D_s8pri,DgofD_s8pri] = kernelest_w(p_swa8(end,5),s_swa8(:,12),sigma,npts,[],s_swa8(:,9));
[D_s9pri,DgofD_s9pri] = kernelest_w(p_swa9(end,5),s_swa9(:,12),sigma,npts,[],s_swa9(:,9));
[D_dpri,DgofD_dpri]   = kernelest(p_dsa(end,5),s_dsa(:,12),sigma,npts,[]);

figure(30)
set(gcf,'color','white')
semilogx(D_s1pri,DgofD_s1pri,'-',D_s4pri,DgofD_s4pri,'--',D_s8pri,DgofD_s8pri,'-.',D_s9pri,DgofD_s9pri,':')
hold on
xlabel('Number of primaries, $n_{\textrm{pri}}$ (-)')
ylabel('$\textrm{d}n/\textrm{dln}(n_{\textrm{pri}})$ (m$^{-3}$)')
ax = gca;
ax.XLim(1) = 0;
ax.XLim(2) = 1e2;
ax.TickLabelInterpreter = 'latex';
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
h = gca;
h.Children(1).LineWidth = 2;
h.Children(2).LineWidth = 2;
h.Children(3).LineWidth = 2;
h.Children(4).LineWidth = 3.5;
h.Children(1).Color = [0      0.7490 0.7490];
h.Children(2).Color = [0.0784 0.1686 0.5490];
h.Children(3).Color = [0.3020 0.7451 0.9333];
h.Children(4).Color = 'r';
saveas(['npri2_5dist_real1' studyid])

disp('Number of primaries distributions, D_KL swa 1, 4, 8, 9')
disp(klest(s_swa1(:,12),s_swa4(:,12),s_swa1(:,9),s_swa4(:,9)))
disp(klest(s_swa1(:,12),s_swa8(:,12),s_swa1(:,9),s_swa8(:,9)))
disp(klest(s_swa1(:,12),s_swa9(:,12),s_swa1(:,9),s_swa9(:,9)))

%% dp - plot lognormal distributions and compute KL divergences from SWA

% Compute lognormal PDF
[D_s1p,DgofD_s1p] = kernelest_w(p_swa1(end,5),s_swa1(:,13),sigma,npts,[],s_swa1(:,9));
[D_s4p,DgofD_s4p] = kernelest_w(p_swa4(end,5),s_swa4(:,13),sigma,npts,[],s_swa4(:,9));
[D_s8p,DgofD_s8p] = kernelest_w(p_swa8(end,5),s_swa8(:,13),sigma,npts,[],s_swa8(:,9));
[D_s9p,DgofD_s9p] = kernelest_w(p_swa9(end,5),s_swa9(:,13),sigma,npts,[],s_swa9(:,9));
[D_dp,DgofD_dp]   = kernelest(p_dsa(end,5),s_dsa(:,13),sigma,npts,[]);

figure(31)
set(gcf,'color','white')
semilogx(D_s1p,DgofD_s1p,'-',D_s4p,DgofD_s4p,'--',D_s8p,DgofD_s8p,'-.',D_s9p,DgofD_s9p,':')
hold on
xlabel('Primary diameter, $d_{\textrm{p}}$ (nm)')
ylabel('$\textrm{d}n/\textrm{dln}(d_{\textrm{p}})$ (m$^{-3}$)')
ax = gca;
ax.XLim(1) = 0.49;
ax.YLim(2) = 8e16;
ax.TickLabelInterpreter = 'latex';
l = legend(leg);
l.Location = 'NorthEast';
l.Interpreter = 'latex';
h = gca;
h.Children(1).LineWidth = 2;
h.Children(2).LineWidth = 2;
h.Children(3).LineWidth = 2;
h.Children(4).LineWidth = 3.5;
h.Children(1).Color = [0      0.7490 0.7490];
h.Children(2).Color = [0.0784 0.1686 0.5490];
h.Children(3).Color = [0.3020 0.7451 0.9333];
h.Children(4).Color = 'r';
saveas(['dpri2_5dist_real1' studyid])

disp('Primary diameter distributions, D_KL swa 1, 4, 8, 9')
disp(klest(s_swa1(:,13),s_swa4(:,13),s_swa1(:,9),s_swa4(:,9)))
disp(klest(s_swa1(:,13),s_swa8(:,13),s_swa1(:,9),s_swa8(:,9)))
disp(klest(s_swa1(:,13),s_swa9(:,13),s_swa1(:,9),s_swa9(:,9)))
