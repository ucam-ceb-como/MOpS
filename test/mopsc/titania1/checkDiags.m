clc,
clear, close all

%% Setup

% set defaults
% set(0,'defaulttextinterpreter','latex')
% set(0,'defaultaxesfontname','Times')
% set(0,'defaulttextfontname','Times')
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
projdir  = 'Network temperature dependence\Recent_work\';
imagedir = 'figures\AIWSWA\real_comp\3ms_c\';
savefigs = 0;

%% Regular diagnostics

% initialise saveas function if saving figures,
% otherwise do nothing here
if (savefigs)
    disp('Saving figures - continuing may overwrite existing!')
    pause
    disp('Continuing...')
    saveas = @(str)(print([projpath projdir imagedir str],'-depsc'));
else
    saveas = @(str)(disp('Not saving'));
end

% load data
folder = [basedir filedir];
pdiags = csvread([folder 'Part-split-diagnostics(stage2).csv'],1);
cdiags = csvread([folder 'Chem-split-diagnostics(stage2).csv'],1);

figure(1)
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
plot(pdiags(:,1)*1000,pdiags(:,8))
hold on
plot(pdiags(:,1)*1000,pdiags(:,9),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Total weight (-)')

subplot(234)
plot(pdiags(:,1)*1000,pdiags(:,14))
hold on
plot(pdiags(:,1)*1000,pdiags(:,15),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('NPN (-)')

subplot(235)
plot(pdiags(:,1)*1000,pdiags(:,10))
hold on
plot(pdiags(:,1)*1000,pdiags(:,11),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Average diameter (-)')

%% Rates

figure(3)
set(gcf,'color','white')
subplot(331)
semilogy(pdiags(:,1)*1000,sum(pdiags(:,16:16+105),2))
hold on
xlabel('Time (ms)')
ylabel('Inc. events (-)')

subplot(332)
semilogy(pdiags(:,1)*1000,pdiags(:,16+106))
hold on
xlabel('Time (ms)')
ylabel('Surf. growth events (-)')

subplot(333)
semilogy(pdiags(:,1)*1000,sum(pdiags(:,16+107:end-2),2))
hold on
xlabel('Time (ms)')
ylabel('Coag. events (-)')

%% Cumulative rates

figure(4)
set(gcf,'color','white')
subplot(141)
semilogy(pdiags(:,1)*1000,cumsum(sum(pdiags(:,16:16+105),2)))
hold on
xlabel('Time (ms)')
ylabel('Inc. events (-)')

subplot(142)
semilogy(pdiags(:,1)*1000,cumsum(pdiags(:,16+106)))
hold on
xlabel('Time (ms)')
ylabel('Surf. growth events (-)')

subplot(143)
semilogy(pdiags(:,1)*1000,cumsum(sum(pdiags(:,16+107:end-2),2)))
hold on
xlabel('Time (ms)')
ylabel('Coag. events (-)')

%% Total rates

nterms = size(pdiags(:,16+107:end-2),2);
nterms = nterms/2;
nreals = nterms/2;

if nreals == 6
    titl_d = {'SF1: (U,U)',...
        'SF2: ($d_c$,$\frac{1}{d_c}$)',...
        'SF3: (U,$\frac{1}{d_c}$)',...
        'SF4: ($d_c$,$\frac{1}{d_c^{2}}$)',...
        'FM1: (U,$\frac{d_c^2}{\sqrt{m}}$)',...
        'FM2: ($d_c^2$,$\frac{1}{\sqrt{m}}$)'};
    ncols  = 3;
else
    titl_d = {'FM1: (U,$\frac{d_c^2\cdot w}{\sqrt{m}}$)',...
        'FM2: ($d_c^2$,$\frac{w}{\sqrt{m}}$)',...
        'FM3: ($\frac{1}{\sqrt{m}}$,$d_c^2\cdot w$)',...
        'FM4: ($\frac{d_c^2}{\sqrt{m}}$,$w$)',...
        'SF1: (U,$w$)',...
        'SF2: ($d_c$,$\frac{w}{d_c}$)',...
        'SF3: ($\frac{1}{d_c}$,$d_c\cdot w$)',...
        'SF4: (U,$\frac{w}{d_c}$)',...
        'SF5: ($d_c$,$\frac{w}{d_c^{2}}$)',...
        'SF6: ($\frac{1}{d_c^2}$,$w\cdot d_c$)',...
        'SF7: ($\frac{1}{d_c}$,$w$)'};
    ncols  = 4;
end

[maxval,ihigh] = max(sum(pdiags(:,16+107:16+107+nreals-1),1));
[minval,ilow]  = min(sum(pdiags(:,16+107:16+107+nreals-1),1));

figure(5)
set(gcf,'color','white')
for i=1:nreals
    subplot(ceil(nreals/ncols),ncols,i)
    plot(pdiags(:,1)*1000,pdiags(:,16+107+(i-1)))
    hold on
    plot(pdiags(:,1)*1000,pdiags(:,16+107+(i-1)+nterms))
    xlabel('Time (ms)')
    ylabel('Number of events (-)')
    if i==ihigh
        c = 'r';
    elseif i==ilow
        c = 'b';
    else
        c = 'k';
    end
    title(titl_d{i},'color',c,'Interpreter','latex')
    legend('Real','Fictitious','Location','NorthWest')
end

%% Rate fractions

totevs = sum(pdiags(:,8:end-2),2);

figure(4)
subplot(144)
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

figure(3)
set(gcf,'color','white')
subplot(334)
semilogy(cdiags(:,1)*1000,cdiags(:,4)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,5)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('TiCl_4 conc. (mol/m^3)')

subplot(335)
semilogy(cdiags(:,1)*1000,cdiags(:,40)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,41)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('O_2 conc. (mol/m^3)')

subplot(336)
semilogy(cdiags(:,1)*1000,cdiags(:,58)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,59)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Ar conc. (mol/m^3)')

subplot(337)
semilogy(cdiags(:,1)*1000,cdiags(:,42)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,43)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Cl_2 conc. (mol/m^3)')

subplot(338)
semilogy(cdiags(:,1)*1000,cdiags(:,20)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,21)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('TiO_2Cl_3 conc. (mol/m^3)')

subplot(339)
semilogy(cdiags(:,1)*1000,cdiags(:,44)*1e6)
hold on
semilogy(cdiags(:,1)*1000,cdiags(:,45)*1e6,':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('O conc. (mol/m^3)')

figure(7)
plot(cdiags(:,1)*1000,cdiags(:,62))
hold on
plot(cdiags(:,1)*1000,cdiags(:,63),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('Temp. (K)')
