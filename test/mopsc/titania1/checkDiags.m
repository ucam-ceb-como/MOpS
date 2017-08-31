clc, clear, close all

set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

pdiags = csvread('Part-split-diagnostics(stage1).csv',1);
cdiags = csvread('Chem-split-diagnostics(stage1).csv',1);

%% Update these to plot weights
wmax=2000;
wmin=1;
nmax=2048;

%% Plot output
% Linear scaling
al = 0;
bl = (wmax-wmin)/(nmax-1);
cl = 1-(wmax-wmin)/(nmax-1);

% Quadratic scaling
aq = (wmax-wmin)/(nmax^2-2*nmax+1);
bq = -2*aq;
cq = wmin-aq-bq;

% Exponential scaling
be = log(wmax/wmin)/(nmax-1);
ae = wmin * exp(-be);
ce = 0;

figure(1)
set(gcf,'color','white')
subplot(131)
plot(pdiags(:,1)*1000,pdiags(:,4))
hold on
plot(pdiags(:,1)*1000,pdiags(:,5),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('SV (-)')
subplot(132)
plot(pdiags(:,1)*1000,pdiags(:,6))
hold on
plot(pdiags(:,1)*1000,pdiags(:,7),':')
legend('Pre','Post','location','North','orientation','horizontal')
xlabel('Time (ms)')
ylabel('NSP (-)')
subplot(133)
plot(pdiags(:,1)*1000,pdiags(:,6).*pdiags(:,6)*aq+pdiags(:,6)*bq+cq)
hold on
plot(pdiags(:,1)*1000,pdiags(:,6).*pdiags(:,6)*aq+pdiags(:,7)*bq+cq,':')
plot(pdiags(:,1)*1000,pdiags(:,6).*pdiags(:,6)*al+pdiags(:,6)*bl+cl,'--')
plot(pdiags(:,1)*1000,pdiags(:,6).*pdiags(:,6)*al+pdiags(:,7)*bl+cl,'-.')
plot(pdiags(:,1)*1000,ae*exp(pdiags(:,6)*be),'-')
plot(pdiags(:,1)*1000,ae*exp(pdiags(:,6)*be),':')
plot([pdiags(1,1)*1000 pdiags(end,1)*1000],[wmax wmax],'r--')
plot([pdiags(1,1)*1000 pdiags(end,1)*1000],[wmin wmin],'g--')
set(gca,'XLim',[0 pdiags(end,1)*1000])
legend('Pre, q','Post, q','Pre, l','Post, l','Pre, e','Post, e','Wmax','Wmin','location','North','orientation','vertical')
xlabel('Time (ms)')
ylabel('Incepting weight (-)')

figure(2)
set(gcf,'color','white')
subplot(331)
plot(pdiags(:,1)*1000,sum(pdiags(:,8:8+105),2))
hold on
xlabel('Time (ms)')
ylabel('Inc. events (-)')
subplot(332)
plot(pdiags(:,1)*1000,pdiags(:,8+106))
hold on
xlabel('Time (ms)')
ylabel('Surf. growth events (-)')
subplot(333)
plot(pdiags(:,1)*1000,sum(pdiags(:,8+107:end-2),2))
hold on
xlabel('Time (ms)')
ylabel('Coag. events (-)')

figure(3)
set(gcf,'color','white')
subplot(131)
loglog(pdiags(:,1)*1000,cumsum(sum(pdiags(:,8:8+105),2)))
xlabel('Time (ms)')
ylabel('Inc. events (-)')
subplot(132)
loglog(pdiags(:,1)*1000,cumsum(pdiags(:,8+106)))
xlabel('Time (ms)')
ylabel('Surf. growth events (-)')
subplot(133)
loglog(pdiags(:,1)*1000,cumsum(sum(pdiags(:,8+107:end-2),2)))
xlabel('Time (ms)')
ylabel('Coag. events (-)')

totevs = sum(pdiags(:,8:end-2),2);

figure(4)
set(gcf,'color','white')
loglog(pdiags(:,1)*1000,sum(pdiags(:,8:8+105)./totevs,2))
hold on
loglog(pdiags(:,1)*1000,pdiags(:,8+106)./totevs,'--')
loglog(pdiags(:,1)*1000,sum(pdiags(:,8+107:end-2),2)./totevs,':')
xlabel('Time (ms)')
ylabel('Fraction of events (-)')
legend('Inc.','Surf. growth','Coag.','location','southwest')

figure(2)
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
