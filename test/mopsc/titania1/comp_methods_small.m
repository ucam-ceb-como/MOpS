% clc, clear, close all

%% Defaults defaults

set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

%% Simulation params

t0      = 0;
tf      = 0.01;
npts    = 1000;
tau     = 0.0013;

%% Inception only

I0      = 1e16;
K0      = 0;

m1inc   = (2*79.87)/(6.022*10^26);
m0f     = @(t,m0)(I0-(1/tau)*m0);
m1f     = @(t,m1)(I0*m1inc-(1/tau)*m1);

OPTIONS = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t1,m0] = ode23s(m0f,[t0 tf],0,OPTIONS);
[t2,m1] = ode23s(m1f,[t0 tf],0,OPTIONS);

data1   = csvread('outputs/incep_10ms_cd/Network(stage1)-part.csv',1);

figure(1)
set(gcf,'color','white')
plot(data1(:,2),data1(:,5),'--')
hold on
tdat = [data1(:,2)' fliplr(data1(:,2)')];
c1 = (data1(:,5)-data1(:,6))';
c2 = (data1(:,5)+data1(:,6))';
inbtw = [c1,fliplr(c2)];
fill(tdat,inbtw,'c','FaceAlpha',0.3,'EdgeColor','none')
plot(t1,m0,'k')
legend('Sim','Sim uncertainty','Numerical','Location','SouthEast')
xlabel('Time (s)')
ylabel('$M_0$ (m$^{-3}$)')
title('Number density: inception')

figure(3)
set(gcf,'color','white')
plot(data1(:,2),data1(:,21),'--')
hold on
tdat = [data1(:,2)' fliplr(data1(:,2)')];
c1 = (data1(:,21)-data1(:,22))';
c2 = (data1(:,21)+data1(:,22))';
inbtw = [c1,fliplr(c2)];
fill(tdat,inbtw,'c','FaceAlpha',0.3,'EdgeColor','none')
plot(t2,m1,'k')
legend('Sim','Sim uncertainty','Numerical','Location','SouthEast')
xlabel('Time (s)')
ylabel('$M_1$ (kg$\cdot$m$^{-3}$)')
title('Mass density: inception')

%% Inception and coagulation

I0      = 1e16;
K0      = 0.01;
tf      = 0.002;

m0f     = @(t,m0)(I0-(1/tau)*m0-(K0/2)*m0^2);
m1f     = @(t,m1)(I0*m1inc-(1/tau)*m1);

[t1,m0] = ode23s(m0f,[t0 tf],0,OPTIONS);
[t2,m1] = ode23s(m1f,[t0 tf],0,OPTIONS);

data2   = csvread('Network(stage1)-part.csv',1);

figure(2)
set(gcf,'color','white')
plot(data2(:,2),data2(:,5),'--')
hold on
tdat = [data2(:,2)' fliplr(data2(:,2)')];
c1 = (data2(:,5)-data2(:,6))';
c2 = (data2(:,5)+data2(:,6))';
inbtw = [c1,fliplr(c2)];
fill(tdat,inbtw,'c','FaceAlpha',0.3,'EdgeColor','none')
plot(t1,m0,'k')
legend('Sim','Sim uncertainty','Numerical','Location','SouthEast')
xlabel('Time (s)')
ylabel('$M_0$ (m$^{-3}$)')
title('Number density: incep and coag')

figure(4)
set(gcf,'color','white')
plot(data2(:,2),data2(:,21),'--')
hold on
tdat = [data2(:,2)' fliplr(data2(:,2)')];
c1 = (data2(:,21)-data2(:,22))';
c2 = (data2(:,21)+data2(:,22))';
inbtw = [c1,fliplr(c2)];
fill(tdat,inbtw,'c','FaceAlpha',0.3,'EdgeColor','none')
plot(t2,m1,'k')
legend('Sim','Sim uncertainty','Numerical','Location','SouthEast')
xlabel('Time (s)')
ylabel('$M_1$ (kg$\cdot$m$^{-3}$)')
title('Mass density: incep and coag')

figure(5)
set(gcf,'color','white')
plot(data2(:,2),data2(:,3),'--')
legend('Sim','Location','SouthEast')
xlabel('Time (s)')
ylabel('$N$ (-)')
title('Number of sim particles: incep and coag')

[~,m0] = ode23s(m0f,linspace(t0,tf,1000),0,OPTIONS);
[t1,m1] = ode23s(m1f,linspace(t0,tf,1000),0,OPTIONS);

figure(6)
set(gcf,'color','white')
plot(data2(:,2),data2(:,5)./data2(:,21),'--')
hold on
plot(t1,m0./m1,'k')
legend('Sim','Numerical','Location','SouthEast')
xlabel('Time (s)')
ylabel('$M_0/M_1$ (kg$^{-1}$)')
title('Ratio $M_0/M_1$: incep and coag')
set(gca,'YScale','log')