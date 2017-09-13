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
basedir  = ''; 
filedir  = '';
filebase = '';

% for saving images
studyid  = '';
projpath = 'C:\Users\Astrid\Documents\Projects\';
projdir  = 'Network temperature dependence\IdeasFromHMMeeting\';
imagedir = 'figures\';
savefigs = 0;

% load data
pdiags = csvread([basedir filedir 'Part-split-diagnostics(stage1).csv'],1);
cdiags = csvread([basedir filedir 'Chem-split-diagnostics(stage1).csv'],1);

% update these to plot correct inception weights
wmax = 1000;
wmin = 1;
nmin = 1;
nmax = 2048;
wtfn = 'off';

% set weight function according to wtfn
if strcmp(wtfn,'L')
    % Linear scaling
    al = 0;
    bl = (wmax-wmin)/(nmax-nmin);
    cl = wmin-(bl*nmin);
    wnew = @(n)(wmin.*(n<=nmin)+(bl*n+cl).*(n>nmin));
elseif strcmp(wtfn,'Q')
    % Quadratic scaling
    aq = (wmax-wmin)/(nmax^2-2*nmax*nmin+nmin^2);
    bq = -2*aq*nmin;
    cq = wmin-(aq*nmin^2)-(bq*nmin);
    wnew = @(n)(wmin.*(n<=nmin)+(aq*n.^2+bq*n+cq).*(n>nmin));
elseif strcmp(wtfn,'E')    
    % Exponential scaling
    be = log(wmax/wmin)/(nmax-nmin);
    ae = wmin * exp(-be*nmin);
    ce = 0;
    wnew = @(n)(wmin.*(n<=nmin)+(ae*exp(be*n)).*(n>nmin));
else
    % No scaling
    wnew = @(n)(1.0);
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

%% SV, NSP, WT

figure(100)
set(gcf,'color','white')
subplot(231)
plot(pdiags(:,1)*1000,pdiags(:,4))
hold on
plot(pdiags(:,1)*1000,pdiags(:,5),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('SV (-)')

subplot(232)
plot(pdiags(:,1)*1000,pdiags(:,6))
hold on
plot(pdiags(:,1)*1000,pdiags(:,7),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('NSP (-)')

subplot(233)
plot(pdiags(:,1)*1000,wnew(pdiags(:,6)))
hold on
plot(pdiags(:,1)*1000,wnew(pdiags(:,7)),':')
plot([pdiags(1,1)*1000 pdiags(end,1)*1000],[wmax wmax],'r--')
plot([pdiags(1,1)*1000 pdiags(end,1)*1000],[wmin wmin],'g--')
set(gca,'XLim',[0 pdiags(end,1)*1000])
legend('Pre','Post','Wmax','Wmin','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Incepting weight (-)')
saveas(['sv_nsp' studyid])

subplot(234)
plot(pdiags(:,1)*1000,pdiags(:,8))
hold on
plot(pdiags(:,1)*1000,pdiags(:,9),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Total statistical weight (-)')

subplot(235)
plot(pdiags(:,1)*1000,pdiags(:,10))
hold on
plot(pdiags(:,1)*1000,pdiags(:,11),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Total particle mass (kg)')

subplot(236)
plot(pdiags(:,1)*1000,pdiags(:,14))
hold on
plot(pdiags(:,1)*1000,pdiags(:,15),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Incepting factor (-)')

figure(1000)
set(gcf,'color','white')
[vals,bins] = hist(pdiags(:,15));
plot(bins,vals/max(vals))
hold on
xlabel('Incepting factor (-)')
ylabel('Count divide by max count (-)')

%% Rates

figure(200)
set(gcf,'color','white')
subplot(331)
plot(pdiags(:,1)*1000,sum(pdiags(:,16:16+105),2))
hold on
xlabel('Time (ms)')
ylabel('Inc. events (-)')

subplot(332)
plot(pdiags(:,1)*1000,pdiags(:,16+106))
hold on
xlabel('Time (ms)')
ylabel('Surf. growth events (-)')

subplot(333)
plot(pdiags(:,1)*1000,sum(pdiags(:,16+107:end-2),2))
hold on
xlabel('Time (ms)')
ylabel('Coag. events (-)')
% saveas(['rates' studyid])

%% Cumulative rates

figure(300)
set(gcf,'color','white')
subplot(131)
loglog(pdiags(:,1)*1000,cumsum(sum(pdiags(:,16:16+105),2)))
hold on
xlabel('Time (ms)')
ylabel('Inc. events (-)')

subplot(132)
loglog(pdiags(:,1)*1000,cumsum(pdiags(:,16+106)))
hold on
xlabel('Time (ms)')
ylabel('Surf. growth events (-)')

subplot(133)
loglog(pdiags(:,1)*1000,cumsum(sum(pdiags(:,16+107:end-2),2)))
hold on
xlabel('Time (ms)')
ylabel('Coag. events (-)')
saveas(['cum_rates' studyid])

%% Rate fractions

totevs = sum(pdiags(:,8:end-2),2);

figure(400)
set(gcf,'color','white')
loglog(pdiags(:,1)*1000,sum(pdiags(:,16:16+105)./totevs,2))
hold on
loglog(pdiags(:,1)*1000,pdiags(:,16+106)./totevs,'--')
loglog(pdiags(:,1)*1000,sum(pdiags(:,16+107:end-2),2)./totevs,':')
xlabel('Time (ms)')
ylabel('Fraction of events (-)')
legend('Inc.','Surf. growth','Coag.','location','southwest')
saveas(['frac_rates' studyid])

%% Chem

figure(200)
set(gcf,'color','white')
subplot(334)
semilogy(cdiags(:,1)*1000,cdiags(:,4)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,5)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('TiCl$_4$ conc. (mol$\cdot$m$^{-3}$)')

subplot(335)
semilogy(cdiags(:,1)*1000,cdiags(:,40)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,41)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('O$_2$ conc. (mol$\cdot$m$^{-3}$)')

subplot(336)
semilogy(cdiags(:,1)*1000,cdiags(:,58)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,59)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Ar conc. (mol$\cdot$m$^{-3}$)')

subplot(337)
semilogy(cdiags(:,1)*1000,cdiags(:,42)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,43)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Cl$_2$ conc. (mol$\cdot$m$^{-3}$)')

subplot(338)
semilogy(cdiags(:,1)*1000,cdiags(:,20)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,21)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('TiO$_2$Cl$_3$ conc. (mol$\cdot$m$^{-3}$)')

subplot(339)
semilogy(cdiags(:,1)*1000,cdiags(:,44)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,45)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('O conc. (mol$\cdot$m$^{-3}$)')
saveas(['chem' studyid])
