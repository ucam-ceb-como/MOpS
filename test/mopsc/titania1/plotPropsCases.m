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
basedir  = 'vienna/real/comparisons/15ms/csvs/';
filedir  = '';
filecmp  = '';
filebase = 'Network(stage1)';

% for plotting psl data
psltime  = '0.015';
npslbins = 10;

% for saving images
studyid  = '';
projpath = 'C:\Users\Astrid\Documents\Projects\';
projdir  = 'Network temperature dependence\Recent_work\';
imagedir = 'figures\AIWSWA\real_comp\3ms_h\';
savefigs = 0;
leg_vals = '';

% order = [1 3 5 2 4 10 12 17 13 6 11 16 7 8 9 14 15]; % reverse way
% order = [1 4 2 5 3 10 13 14 15 6 11 7 9 16 17 12 8]; % this way

leg = {'D1','D2','D3','D4','D5',...
       'S1','S2','S3','S4','S5','S6',...
       'S7','S8','S9','S10','S11','S12'};
tickval = {'D1','D2','D3','D4','D5',...
           'S1','S2','S3','S4','S5','S6',...
           'S7','S8','S9','S10','S11','S12'};
texts = {'DSA','PSI$_1$','PSI$_2$','HCI$_1$','HCI$_2$',...
         'SWA','LAIW$_1$','QAIW$_1$','EAIW$_1$','EAIW$_2$',...
         'EAIW$_1+$PSI$_1$','EAIW$_1+$PSI$_2$','EAIW$_1+$EPSI$_2$',...
         'EAIW$_1+$WPSI$_2$','EAIW$_1+$BPSI$_2$','EAIW$_1+$HCI$_1$',...
         'EAIW$_1+$HCI$_2$'};

j=0;       
for i = [1:5,6,8:9,11:17]
j=i;
if i < 6
    filedir = ['dsa' num2str(i) '/']; 
    filecmp = 'swa1/';   
else
    filedir = ['swa' num2str(i-5) '/']; 
    filecmp = 'swa1/';
end

% load all data
m0dat = csvread([basedir filecmp filebase '-part.csv'],1);
part = csvread([basedir filedir filebase '-part.csv'],1);
rate = csvread([basedir filedir filebase '-part-rates.csv'],1);
pslf = csvread([basedir filedir filebase '-psl(' psltime 's).csv'],1);
chem = csvread([basedir filedir filebase '-chem.csv'],1);
cput = csvread([basedir filedir filebase '-cput.csv'],1);

% initialise saveas function if saving figures, otherwise do nothing here
if (savefigs && (j==17))
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
% 
% col1 = lines(7);
% col2 = [0 0 1;1 1 0;1 0 1;0 1 0;0 1 1;1 0 0];
% col3 = abs(ones(7,3)-col1-0.75*ones(7,3));
% cols = [col1;col2;col3];

col1 = lines(7);
col2 = [0 0 1;1 1 0;1 0 1;0 1 0;0 1 1;1 0 0];
col3 = abs(ones(7,3)-col1-0.75*ones(7,3));
cols = [col1(1,:);col1(2,:);col1(2,:);col1(3,:);col1(3,:);
        col1(4,:);col1(5,:);col1(5,:);col1(5,:);
        col1(6,:);col1(7,:);col1(7,:);col2(1,:);col2(1,:);
        col2(1,:);col2(2,:);col2(2,:)];
lins = {'-','--','-.',':','-','--','-.',':','-','--','-.',':','-',...
        '--','-.',':','-'};
    
figure(100)
set(gcf,'color','white')
bar(j,mean(part(:,3)),'FaceColor',cols(j,:))
hold on
errorbar(j,mean(part(:,3)),mean(part(:,4)),'color','k')
if j==17
    xlabel('Time (ms)','color','k')
    ylabel('Average number of particles (-)')
    addlegend(leg_vals);
    cg = [0.45 0.45 0.45];
    ax=gca;
    line([5.5 5.5],ax.YLim,'color',cg,'LineStyle',':')
    text(2.5,0.92*ax.YLim(2),'DSA','Fontsize',14,'color',cg)
    text(11,0.92*ax.YLim(2),'SWA','Fontsize',14,'color',cg)
    set(gca,'XTick',1:17,'XTickLabel',tickval,'XTickLabelRotation',90)
    saveas(['nspm' studyid])
end

figure(1000)
set(gcf,'color','white')
plot(part(:,2)*1000,part(:,3),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Number of particles (-)')
if j==17
    l = legend(leg);
    l.Location = 'EastOutside';
end
saveas(['nsp_' studyid])

figure(1)
set(gcf,'color','white')
subplot(231)
plot(part(:,2)*1000,part(:,3),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('NSP (-)')
addlegend(leg_vals);
saveas(['nsp' studyid])

subplot(232)
% figure(2)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,5),'color',cols(j,:),'LineStyle',lins{j})
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
plot(part(:,2)*1000,part(:,9)*1e9,'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Collision diameter (nm)')
addlegend(leg_vals);
saveas(['dcol_ave' studyid])

subplot(234)
% figure(4)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,21),'color',cols(j,:),'LineStyle',lins{j})
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
plot(part(:,2)*1000,part(:,37),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Number primaries per particle (-)')
addlegend(leg_vals);
saveas(['npri_ave' studyid])

subplot(236)
% figure(6)
% set(gcf,'color','white')
plot(part(:,2)*1000,part(:,39)*1e9,'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Primary average diameter (nm)')
addlegend(leg_vals);
saveas(['dpri_ave' studyid])

figure(60)
set(gcf,'color','white')
semilogy(part(:,2)*1000,part(:,33)*1e6,'color',cols(j,:),'LineStyle',lins{j})
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
plot(part(:,2)*1000,abs(part(:,5)-m0dat(:,5))./m0dat(:,5),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Relative difference in m0 (-)')
saveas(['m0_err' studyid])

figure(91)
set(gcf,'color','white')
plot(part(:,2)*1000,abs(part(:,39)-m0dat(:,39))./m0dat(:,39),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Relative difference in primary diameter (-)')
saveas(['dpri_err' studyid])

figure(92)
set(gcf,'color','white')
plot(part(:,2)*1000,abs(part(:,21)-m0dat(:,21))./m0dat(:,21),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Relative difference in m1 (-)')
saveas(['m1_err' studyid])

figure(900)
set(gcf,'color','white')
ybar = mean(abs(part(end-50:end,5)-m0dat(end-50:end,5))./m0dat(end-50:end,5));
bar(j,ybar,'FaceColor',cols(j,:))
hold on
% text(j,ybar,texts{j},'color',cols(j,:),'rotation',90,'Fontsize',14)
if j==17
    xlabel('Time (ms)')
    ylabel('Final relative difference in m0 (-)')
    set(gca,'XTick',1:17,'XTickLabel',tickval,'XTickLabelRotation',90)
    cg = [0.45 0.45 0.45];
    ax=gca;
    text(2.5,0.92*ax.YLim(2),'DSA','Fontsize',14,'color',cg)
    text(11,0.92*ax.YLim(2),'SWA','Fontsize',14,'color',cg)
    ax=gca;
    line([5.5 5.5],ax.YLim,'color',cg,'LineStyle',':')
    saveas(['m0_errf' studyid])
end

figure(910)
set(gcf,'color','white')
bar(j,mean(abs(part(end-50:end,39)-m0dat(end-50:end,39))./m0dat(end-50:end,39)),'FaceColor',cols(j,:))
hold on
if j==17
    xlabel('Time (ms)')
    ylabel('Final relative difference in ${d}_{\textrm{p}}$ (-)')
    set(gca,'XTick',1:17,'XTickLabel',tickval,'XTickLabelRotation',90)
    cg = [0.45 0.45 0.45];
    ax=gca;
    text(2.5,0.92*ax.YLim(2),'DSA','Fontsize',14,'color',cg)
    text(11,0.92*ax.YLim(2),'SWA','Fontsize',14,'color',cg)
    ax=gca;
    line([5.5 5.5],ax.YLim,'color',cg,'LineStyle',':')
    saveas(['dpri_errf' studyid])
end

figure(920)
set(gcf,'color','white')
bar(j,mean(abs(part(end-50:end,21)-m0dat(end-50:end,21))./m0dat(end-50:end,21)),'FaceColor',cols(j,:))
hold on
if j==17
    xlabel('Time (ms)')
    ylabel('Final relative difference in m1 (-)')
    set(gca,'XTick',1:17,'XTickLabel',tickval,'XTickLabelRotation',90)
    cg = [0.45 0.45 0.45];
    ax=gca;
    text(2.5,0.92*ax.YLim(2),'DSA','Fontsize',14,'color',cg)
    text(11,0.92*ax.YLim(2),'SWA','Fontsize',14,'color',cg)
    ax=gca;
    line([5.5 5.5],ax.YLim,'color',cg,'LineStyle',':')
    saveas(['m1_errf' studyid])
end

figure(930)
set(gcf,'color','white')
bar(j,mean(abs(part(end-50:end,9)-m0dat(end-50:end,9))./m0dat(end-50:end,9)),'FaceColor',cols(j,:))
hold on
if j==17
    xlabel('Time (ms)')
    ylabel('Final relative difference in ${d}_{\textrm{c}}$ (-)')
    set(gca,'XTick',1:17,'XTickLabel',tickval,'XTickLabelRotation',90)
    cg = [0.45 0.45 0.45];
    ax=gca;
    text(2.5,0.92*ax.YLim(2),'DSA','Fontsize',14,'color',cg)
    text(11,0.92*ax.YLim(2),'SWA','Fontsize',14,'color',cg)
    ax=gca;
    line([5.5 5.5],ax.YLim,'color',cg,'LineStyle',':')
    saveas(['dcol_errf' studyid])
end

%% PSL

if j<6
    figure(2001)
    set(gcf,'color','white')
    [n_d,x_d] = histwc(pslf(:,13),pslf(:,9),npslbins);
    plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
    hold on
end
if j==17
    figure(2001)
    xlabel('Primary diameter (nm)')
    ylabel('Count divided by max. count (-)')
    % addlegend(leg_vals);
    l = legend(leg{1},leg{2},leg{3},leg{4},leg{5});
    l.Location ='EastOutside';
    ax = gca;
    ax.XLim(1) = 0.0;
    saveas(['dpri_psd_large1' studyid])
end
if j==1 || j==6 || j==7 || j==8 || j==9 || j==10
    figure(2002)
    set(gcf,'color','white')
    [n_d,x_d] = histwc(pslf(:,13),pslf(:,9),npslbins);
    plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
    hold on
end
if j==17
    figure(2002)
    xlabel('Primary diameter (nm)')
    ylabel('Count divided by max. count (-)')
    % addlegend(leg_vals);
    l = legend(leg{1},leg{6},leg{7},leg{8},leg{9},leg{10});
    l.Location ='EastOutside';
    ax = gca;
    ax.XLim(1) = 0.0;
    saveas(['dpri_psd_large2' studyid])
end
if j==1 || j==6 || j==11 || j==12 || j==16 || j==17
    figure(2003)
    set(gcf,'color','white')
    [n_d,x_d] = histwc(pslf(:,13),pslf(:,9),npslbins);
    plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
    hold on
end
if j==17
    figure(2003)
    xlabel('Primary diameter (nm)')
    ylabel('Count divided by max. count (-)')
    % addlegend(leg_vals);
    l = legend(leg{1},leg{6},leg{11},leg{12},leg{16},leg{17});
    l.Location ='EastOutside';
    ax = gca;
    ax.XLim(1) = 0.0;
    saveas(['dpri_psd_large3' studyid])
end
if j==1 || j==6 || j==13 || j==14 || j==15
    figure(2004)
    set(gcf,'color','white')
    [n_d,x_d] = histwc(pslf(:,13),pslf(:,9),npslbins);
    plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
    hold on
end
if j==17
    figure(2004)
    xlabel('Primary diameter (nm)')
    ylabel('Count divided by max. count (-)')
    % addlegend(leg_vals);
    l = legend(leg{1},leg{6},leg{13},leg{14},leg{15});
    l.Location ='EastOutside';
    ax = gca;
    ax.XLim(1) = 0.0;
    saveas(['dpri_psd_large4' studyid])
end

% set(gcf,'color','white')
% [n_d,x_d] = histwc(pslf(:,13),pslf(:,9),npslbins);
% plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
% hold on
% if j==17
%     xlabel('Primary diameter (nm)')
%     ylabel('Count divided by max. count (-)')
%     % addlegend(leg_vals);
%     l = legend(leg);
%     l.Location ='EastOutside';
%     ax = gca;
%     ax.XLim(1) = 0.0;
%     saveas(['dpri_psd_large' studyid])
% end

figure(2)
set(gcf,'color','white')
[n_d,x_d] = histwc(pslf(:,3),pslf(:,9),npslbins);
subplot(231)
plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
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
plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
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
plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
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
plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
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
plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
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
plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Particle age (s)')
ylabel('Count divided by max. count (-)')
addlegend(leg_vals);
ax = gca;
ax.XLim(1) = 0.0;
saveas(['page_psd' studyid])

if i > 5
figure(20)
set(gcf,'color','white')
[n_d,x_d] = hist(pslf(:,9));
% subplot(231)
plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
hold on
    if j==17
        xlabel('Weight (-)')
        ylabel('Count divided by max. count (-)')
        % addlegend(leg_vals);
        l = legend(leg{6:end});
        l.Location ='EastOutside';
        ax = gca;
        ax.XLim(1) = 0.0;
        set(gca,'XScale','log')
        saveas(['weight_psd' studyid])
    end
end

figure(30)
set(gcf,'color','white')
% subplot(131)
[n_d,x_d] = histwc(pslf(:,11),pslf(:,9),npslbins);
plot(x_d,n_d/max(n_d),'color',cols(j,:),'LineStyle',lins{j})
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
plot(chem(:,2)*1000,chem(:,3)*1e6,'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('TiCl$_4$ conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['ticl4' studyid])

subplot(232)
% figure(2)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,39)*1e6,'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('O$_2$ conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['o2' studyid])

subplot(233)
% figure(3)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,45)*1e6,'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Cl$_2$ conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['cl2' studyid])

subplot(234)
% figure(4)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,19)*1e6,'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('TiO$_2$Cl$_3$ conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['tio2cl3' studyid])

subplot(235)
% figure(5)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,57)*1e6,'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Ar conc. (mol$\cdot$m$^{-3}$)')
addlegend(leg_vals);
saveas(['ar' studyid])

subplot(236)
% figure(6)
% set(gcf,'color','white')
plot(chem(:,2)*1000,chem(:,61),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Temperature (K)')
addlegend(leg_vals);
saveas(['temp' studyid])

%% CPUT

cadj = [0;(cput(2:end,3)-cput(2,3))/60];
figure(4)
set(gcf,'color','white')
plot(cput(:,2)*1000,cadj,'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Solver time (min)')
addlegend(leg_vals);
saveas(['cput' studyid])

figure(400)
set(gcf,'color','white')
bar(j,cadj(end),'FaceColor',cols(j,:))
hold on
% errorbar(j,cadj(end),cput(end,4)/60,'color','k')
if j==17
    xlabel('Time (ms)','color','k')
    ylabel('Final solver time (min)')
    addlegend(leg_vals);
    cg = [0.45 0.45 0.45];
    ax=gca;
    line([5.5 5.5],ax.YLim,'color',cg,'LineStyle',':')
    text(2.5,0.92*ax.YLim(2),'DSA','Fontsize',14,'color',cg)
    text(11,0.92*ax.YLim(2),'SWA','Fontsize',14,'color',cg)
%     annotation('arrow',[0.25 0.3],[0.9*ax.YLim(2) 0.9*ax.YLim(2)],'color',cg)
%     annotation('arrow',[0.3 0.35],[0.9*ax.YLim(2) 0.9*ax.YLim(2)],'color',cg)
    set(gca,'XTick',1:17,'XTickLabel',tickval,'XTickLabelRotation',90)
    saveas(['cputf' studyid])
end

%% Rates

figure(5)
set(gcf,'color','white')
subplot(221)
plot(rate(:,2)*1000,sum(rate(:,3:2:211),2),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Inc. rate (m$^{-3}\cdot$s$^{-1}$)')
addlegend(leg_vals);
saveas(['inc_rate' studyid])

subplot(222)
% figure(3)
% set(gcf,'color','white')
plot(rate(:,2)*1000,rate(:,215),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Surf. growth rate (m$^{-3}\cdot$s$^{-1}$)')
addlegend(leg_vals);
saveas(['sg_rate' studyid])

subplot(223)
% figure(3)
% set(gcf,'color','white')
plot(rate(:,2)*1000,sum(rate(:,217:2:end),2),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('Coag. rate (m$^{-3}\cdot$s$^{-1}$)')
addlegend(leg_vals);
saveas(['coag_rate' studyid])

subplot(224)
% figure(2)
% set(gcf,'color','white')
plot(rate(:,2)*1000,rate(:,213),'color',cols(j,:),'LineStyle',lins{j})
hold on
xlabel('Time (ms)')
ylabel('No inc. rate (m$^{-3}\cdot$s$^{-1}$)')
addlegend(leg_vals); 
saveas(['no_inc_rate' studyid])
end