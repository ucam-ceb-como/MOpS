clc, clear, close all

set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

Fo2 = [1015.6,1002.4,979.51,1002.9,1026.1,1045.9,1057.1,1057.6,...
       1084.1,1117.4,1156.2,1198.7,1246.1,1317.7,1394.9,1475.6,1562.2,...
       1641,1625.3,... % pfr
       1613.5,1604.2,1596.8,1590.4,1584.5,1578.9,1573.3,...
       1567.7,1562,1556.2,1550.1,1543.7,1537,1530,1522.5,1514.6,1506.1,...
       1496.8,1486.9,1476.1,1464.3,1451.3,1436.5,1420,1401.4,1379.8,...
       1354.5,1324.4,1287.4,1243.5,1191.7,1133.5,1072.5,1011.4,951.67,...
       892.05,834.75,778.7,724.21,671.69,621.67,574.93,538.28,517.95,...
       502.86,491.19,481.76,473.73,466.41,459.77,453.46,447.38,441.49,...
       435.78,430.15,424.82,419.7,414.81,410.16,405.68,401.54,397.69,...
       394.12,390.87,387.95,385.29,382.85,380.63,378.63,376.8,375.1,...
       373.56,372.14,370.84,369.64,368.53,367.49,366.5,365.56,364.65,];
Fticl4 = [33.722,153.8,248.77,291.1,333.16,371.73,401.71,...
          423.08,482.6,549.06,620.48,696.01,777.71,891.31,1011.4,...
          1134.3,1263,1378.7,1363.1,...% pfr
          1351.2,1341.9,1334.6,1328.1,1322.3,...
          1316.7,1311,1305.4,1299.8,1293.9,1287.9,1281.4,1274.8,...
          1267.7,1260.3,1252.4,1243.9,1234.6,1224.7,1213.9,1202.1,...
          1189,1174.2,1157.8,1139.1,1117.5,1092.2,1062.2,1025.1,...
          981.28,929.46,871.3,810.21,749.21,689.43,629.81,572.51,...
          516.47,461.98,409.45,359.44,312.7,276.05,255.72,240.63,...
          228.97,219.54,211.5,204.19,197.55,191.23,185.16,179.27,...
          173.56,167.93,162.6,157.49,152.6,147.95,143.47,139.33,...
          135.48,131.91,128.66,125.74,123.08,120.65,118.42,116.43,...
          114.6,112.89,111.35,109.94,108.64,107.44,106.33,105.29,...
          104.3,103.35,102.45];
Fcl2 = [2.78E-03,26.101,84.183,138.85,193.54,254.46,331.79,...
        439.52,558.76,667.17,763.46,855.69,945.2,1028.9,1106.7,1173.5,...
        1228.8,1271.9,1304.9,... % pfr
        1329.7,1349,1364.1,1377.2,1389,1400.3,...
        1411.6,1422.8,1434.1,1445.7,1457.8,1470.7,1483.9,1497.9,...
        1512.7,1528.4,1545.3,1563.8,1583.4,1604.8,1628.2,1654.1,...
        1683.4,1715.9,1752.8,1795.3,1845,1903.9,1975.9,2060.7,...
        2160,2270.4,2385.3,2499.7,2611.4,2722.7,2829.3,2933.2,3033.9,...
        3130.5,3222.1,3307.1,3966.9,4005.4,4034.7,4057.7,4076.5,4092.7,...
        4107.4,4120.7,4133.3,4145.3,4156.8,4167.9,4178.8,4189,4198.8,...
        4208,4216.9,4225.3,4233.2,4240.5,4247.4,4253.7,4259.4,4264.6,...
        4269.4,4273.9,4277.9,4281.6,4285.1,4288.3,4291.3,4294.1,4296.7,...
        4299.1,4301.5,4303.8,4306,4308.2];
    
eqreac1 = linspace(0,1,4);  % 1
eqreac2 = linspace(1,4,7);  % 2 3 4 
eqreac3 = linspace(4,7,6);  % 5 6 7
eqreac4 = linspace(7,10,6); % 8 9 10
eqreac5 = linspace(10,11,81); % 11
eqreac  = [eqreac1(2:end),eqreac2(2:end),eqreac3(2:end),eqreac4(2:end),...
           eqreac5(2:end)];

kldivs  = 0;
savfigs = 1;
fdir    = 'C:\Users\Astrid\Documents\Projects\Network temperature dependence\Recent_work\figures\networks\half_EA\';

sid={'1','2','3','4','5','6'};%,'7','8','9','10','11','12'};
n = 6;
x = (1:n);
stage = sid{n};
sigma = 0.07;

c_g = [0.5 0.5 0.5];

if kldivs
    disp('Computing KL divergences...')
    npts = 20000;
else
    disp(['Plotting profiles and distributions (for stage ' stage ')...'])
    npts = 1000;
end

if savfigs
    disp('Saving figures - continuing will overwrite existing!')
    pause 
    disp('Continuing...')
    saveas = @(str)(print([fdir str],'-depsc'));
else
    saveas = @(str)(disp('Not saving'));
end

%% start
s1=csvread(['../hm_revsd_halfEa/Network(stage12)-psl(1.5s).csv'],1);
p1=csvread(['../hm_revsd_halfEa/Network(stage12)-part.csv'],1);

if ~kldivs
    c=zeros(n,1);
    ca=c;
    cb=c;
    cc=c;
    cd=c;
    ce=c;
    cf=c;
    cg=c;
    ch=c;
    ci=c; % TiCl4
    cj=c; % O2
    ck=c; % Cl2
    ct=0;
    for i=1:n
        if i < n-1
            fn = sid{i};
        elseif i < n
            fn = '11';
        else
            fn = '12';
        end
        c1=csvread(['../hm_revsd_halfEa/Network(stage' fn ')-chem.csv'],1);
        c(i)=c1(end,61);
        ci(i)=c1(end,3);
        cj(i)=c1(end,39);
        ck(i)=c1(end,45);
        c1=csvread(['../hm_revsd_halfEa/Network(stage' fn ')-part.csv'],1);
        cd(i)=c1(end,43); % sr
        ce(i)=c1(end,9);  % dc
        cf(i)=c1(end,39); % dp
        cg(i)=c1(end,37); % np
        ch(i)=c1(end,41); % sl
        c1=csvread(['../hm_revsd_halfEa/Network(stage' fn ')-part-rates.csv'],1);
        ca(i)=sum(c1(end,3:2:212));
        cb(i)=c1(end,215);
        cc(i)=c1(end,217);
        c1=csvread(['../hm_revsd_halfEa/Network(stage' fn ')-cput.csv'],1);
        ct=ct+c1(end,3);
    end
    
    
    %% Conc ratio
    figure(1000)
    set(gcf,'color','white')
    plot(x,cj./ci,'-b+')
    hold on
    plot(eqreac,Fo2./Fticl4,':r')
    ylabel('O$_{2}$:TiCl$_{4}$ ratio (mol$/$mol)')
    xlabel('Reactor index (distance)')
    legend('Sim','HM data')
%     set(gca,'YLim',[600 1100],'XLim',[0 n+1],'XTick',x)
%     saveas(['ticl4_o2' caseno])

    figure(1001)
    set(gcf,'color','white')
    plot(x,cj./ck,'-b+')
    hold on
    plot(eqreac,Fo2./Fcl2,':r')
    ylabel('O$_{2}$:Cl$_{2}$ ratio (mol$/$mol)')
    xlabel('Reactor index (distance)')
    legend('Sim','HM data')
%     set(gca,'YLim',[600 1100],'XLim',[0 n+1],'XTick',x)
%     saveas(['cl2_o2' caseno])

    %% Temperature
    figure(100)
    set(gcf,'color','white')
    plot(x,c,'--bo')
    hold on
    ylabel('Temperature (K)')
    xlabel('Reactor index (distance)')
    for i = 1:n
        if i==1 || i==2 || i==3 || i==4
            text(x(i)*0.95,c(i)*1.05,sid{i},'Interpreter','latex','Background','y','EdgeColor','y')
        elseif i==5
            text(x(i)*0.95,c(i)*1.05,sid{i},'Interpreter','latex','Background','r','EdgeColor','r')
        elseif i==6
            text(x(i)*0.95,c(i)*1.05,sid{i},'Interpreter','latex','Background','c','EdgeColor','c')
        else
            text(x(i)*0.95,c(i)*1.05,sid{i},'Interpreter','latex','Background',c_g,'EdgeColor',c_g,'Color','white')
        end
    end
%     set(gca,'YLim',[600 1100],'XLim',[0 n+1],'XTick',x)
    saveas(['temp'])
    
    %% Process rates 1-3
    figure(1)
    set(gcf,'color','white')
    semilogy(x,ca,'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Inception rate (mol$\cdot$m$^{-3}\cdot$s$^{-1}$)')
%     set(gca,'YScale','log','XLim',[0 x(end)+1],'XTick',x)
%     saveas(['inc' caseno])
    
    figure(2)
    set(gcf,'color','white')
    semilogy(x,cb,'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Surface growth rate (mol$\cdot$m$^{-3}\cdot$s$^{-1}$)')
%     set(gca,'YScale','log','YLim',[10^24 10^28],'XLim',[0 x(end)+1],'XTick',x)
%     saveas(['add' caseno])
    
    figure(3)
    set(gcf,'color','white')
    semilogy(x,cc,'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Coagulation rate (mol$\cdot$m$^{-3}\cdot$s$^{-1}$)')
%     set(gca,'YScale','log','XLim',[0 x(end)+1],'XTick',x)
%     saveas(['coa' caseno])
    
    figure(4)
    set(gcf,'color','white')
    semilogy(x,cd,'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Sintering rate (m$^{2}\cdot$m$^{-3}\cdot$s$^{-1}$)')
%     set(gca,'YScale','log','YLim',[10^-25 10^-15],'XLim',[0 x(end)+1],'XTick',x)
%     saveas(['sin' caseno])
    
    %% Relative process rates 1-3
    figure(5)
    set(gcf,'color','white')
    semilogy(x,ca./(cc+ca+cb),'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Inception rate fraction')
%     set(gca,'YScale','log','XLim',[0 x(end)+1],'XTick',x)
%     saveas(['pinc' caseno])
    
    figure(6)
    set(gcf,'color','white')
    semilogy(x,cb./(cc+ca+cb),'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Surface growth rate fraction')
%     set(gca,'YScale','log','XLim',[0 x(end)+1],'XTick',x)
%     saveas(['padd' caseno])
    
    figure(7)
    set(gcf,'color','white')
    semilogy(x,cc./(cc+ca+cb),'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Coagulation rate fraction')
%     set(gca,'YScale','log','XLim',[0 x(end)+1],'XTick',x)
%     saveas(['pcoa' caseno])
    
    %% Particle properties
    figure(8)
    set(gcf,'color','white')
    plot(x,ce*1e9,'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Collision diameter (nm)')
    set(gca,'XLim',[0 x(end)+1],'XTick',x)
    saveas(['dcol'])
    
    figure(9)
    set(gcf,'color','white')
    plot(x,cf*1e9,'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Primary diameter (nm)')
    set(gca,'XLim',[0 x(end)+1],'XTick',x)
    saveas(['dpri'])
    
    figure(10)
    set(gcf,'color','white')
    plot(x,cg,'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Number of primaries per particle')
    set(gca,'XLim',[0 x(end)+1],'XTick',x)
    saveas(['npri'])
    
    figure(11)
    set(gcf,'color','white')
    plot(x,(ce./cf).^3,'-b+')
    hold on
    xlabel('Reactor index (distance)')
    ylabel('Degree of aggregation')
    set(gca,'XLim',[0 x(end)+1],'XTick',x)
    saveas(['dagg'])
    
    %% CPU time
%     figure(12)
%     set(gcf,'color','white')
%     plot(1,ct/3600,'b+')
%     hold on
%     xlabel('Case')
%     ylabel('Total CPU time (hr)')
%     set(gca,'XLim',[0 8],'XTick',1:7,'XTickLabel',leg,...
%         'TickLabelInterpreter','latex','XTickLabelRotation',25)
%     saveas(['cput' caseno])
end

%% Normal distributions
% figure(2)
% set(gcf,'color','white')
% [f,x,u] = ksdensity(s1(:,3));
% plot(x,f,'-b')
% hold on

% figure(3)
% set(gcf,'color','white')
% [f,x,w] = ksdensity(s1(:,14));
% plot(x,f,'-b')
% hold on

% figure(4)
% set(gcf,'color','white')
% [f,x,v] = ksdensity(s1(:,13));
% plot(x,f,'-b')
% hold on

% figure(5)
% set(gcf,'color','white')
% [f,x,z] = ksdensity(s1(:,15));
% plot(x,f,'-b')
% hold on

% figure(6)
% set(gcf,'color','white')
% [f,x,q] = ksdensity(s1(:,6));
% plot(x,f,'-b')
% hold on

%% Lognormal distributions of particle properties
dvec1=[];
[D,DgofD11] = kernelest(p1(end,5),s1(:,3),sigma,npts,dvec1);
if kldivs
    dvec1 = D;
    fofD11 = (DgofD11./D)/p1(end,5);
else
    figure(20)
    set(gcf,'color','white')
    plot(D,DgofD11,'-b')
    hold on
    xlabel('Collision diameter, $d_{\textrm{c}}$ (nm)')
    ylabel('$\textrm{d}n/\textrm{dln}(d_{\textrm{c}})$ (m$^{-3}$)')
    set(gca,'YScale','log','XScale','log')
    set(gca,'YLim',[1 10^17])
    saveas(['lkdcol_s' stage ''])
end

dvec2=[];
[D,DgofD12] = kernelest(p1(end,5),s1(:,14),sigma,npts,dvec2);
if kldivs
    dvec2 = D;
    fofD12 = (DgofD12./D)/p1(end,5);
else
    figure(30)
    set(gcf,'color','white')
    plot(D,DgofD12,'-b')
    hold on
    xlabel('Sintering level, $s$ (-)')
    ylabel('$\textrm{d}n/\textrm{dln}s$')
    set(gca,'YScale','log')
    set(gca,'XLim',[0 1],'YLim',[1 10^17])
    saveas(['lksl_s' stage ''])
end

dvec3=[];
[D,DgofD13] = kernelest(p1(end,5),s1(:,13),sigma,npts,dvec3);
if kldivs
    dvec3 = D;
    fofD13 = (DgofD13./D)/p1(end,5);
else
    figure(40)
    set(gcf,'color','white')
    plot(D,DgofD13,'-b')
    hold on
    xlabel('Primary particle diameter, $d_{\textrm{p}}$ (nm)')
    ylabel('$\textrm{d}n/\textrm{dln}(d_{\textrm{p}})$ (m$^{-3}$)')
    set(gca,'YScale','log','XScale','log')
    set(gca,'YLim',[1 10^17])
    saveas(['lkdpri_s' stage ''])
end

dvec4=[];
[D,DgofD14] = kernelest(p1(end,5),s1(:,15),sigma,npts,dvec4);
if kldivs
    dvec4 = D;
    fofD14 = (DgofD14./D)/p1(end,5);
else
    figure(50)
    set(gcf,'color','white')
    plot(D,DgofD14,'-b')
    hold on
    xlabel('Total sintering time, $t_s$ (s)')
    ylabel('$\textrm{d}n/\textrm{dln}t_s$')
    set(gca,'YScale','log','XScale','log')
    set(gca,'YLim',[1 10^17])
    saveas(['lktsin_s' stage ''])
end

dvec5=[];
[D,DgofD15] = kernelest(p1(end,5),s1(:,6),sigma,npts,dvec5);
if kldivs
    dvec5 = D;
    fofD15 = (DgofD15./D)/p1(end,5);
else
    figure(60)
    set(gcf,'color','white')
    loglog(D,DgofD15,'-b')
    hold on
    xlabel('Particle age, $t_p$ (s)')
    ylabel('$\textrm{d}n/\textrm{dln}t_p$')
    set(gca,'YLim',[1,10^17])
    saveas(['lkage_s' stage ''])
end

dvec6=[];
[D,DgofD16] = kernelest(p1(end,5),s1(:,12),sigma,npts,dvec6);
if kldivs
    dvec6 = D;
    fofD16 = (DgofD16./D)/p1(end,5);
else
    figure(70)
    set(gcf,'color','white')
    plot(D,DgofD16,'-b')
    hold on
    xlabel('Number of primaries, $n_{\textrm{pri}}$ (s)')
    ylabel('$\textrm{d}n/\textrm{dln}(n_{\textrm{pri}})$ (m$^{-3}$)')
    set(gca,'YScale','log')
    set(gca,'YLim',[1,10^17])
    saveas(['lknpri_s' stage ''])
end

dvec7=[];
[D,DgofD17] = kernelest(p1(end,5),s1(:,11),sigma,npts,dvec7);
if kldivs
    dvec7 = D;
    fofD16 = (DgofD17./D)/p1(end,5);
else
    figure(80)
    set(gcf,'color','white')
    plot(D,DgofD17,'-b')
    hold on
    xlabel('Number of TiO$_2$ units, $n_{\textrm{TiO2}}$ (s)')
    ylabel('$\textrm{d}n/\textrm{dln}(n_{\textrm{TiO2}})$ (m$^{-3}$)')
    set(gca,'YScale','log')
    set(gca,'YLim',[1,10^17])
    saveas(['lkntio2_s' stage ''])
end