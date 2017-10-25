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

runmulti = 0;

% for finding files
basedir  = ''; %vienna/real/comparisons/3ms/csvs_c/
filedir  = '';
filebase = '';

% for saving images
studyid  = '';
projpath = 'C:\Users\Astrid\Documents\Projects\';
projdir  = 'Network temperature dependence\Recent_work\';
imagedir = 'figures\AIWSWA\real_comp\3ms_c\';
savefigs = 1;

% update these to plot correct inception weights
wmax = 500;
wmin = 10;
nmin = 1;
nmax = 4096;
wtfn = 'E';

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

%% Multi-diagnostics

if runmulti
    col1  = lines(7);
    col2  = [0 0 1;1 1 0;1 0 1;0 1 0;0 1 1;1 0 0];
    col3  = abs(ones(7,3)-col1-0.75*ones(7,3));
    %     cols  = [col1;col2;col3];
    cols = [col1(1,:);col1(2,:);col1(2,:);col1(3,:);col1(3,:);
        col1(4,:);col1(5,:);col1(5,:);col1(5,:);
        col1(6,:);col1(7,:);col1(7,:);col2(1,:);col2(1,:);
        col2(1,:);col2(2,:);col2(2,:)];
    lins = {'-','--','-.',':','-','--','-.',':','-','--','-.',':','-',...
        '--','-.',':','-'};
    
    for i=1:17
        % initialise saveas function if saving figures,
        % otherwise do nothing here
        if (savefigs && i == 17)
            disp('Saving figures - continuing may overwrite existing!')
            pause
            disp('Continuing...')
            saveas = @(str)(print([projpath projdir imagedir str],'-depsc'));
        else
            saveas = @(str)(disp('Not saving'));
        end
        % load data
        if i < 6
            filedir = ['dsa' num2str(i) '/'];
        else
            filedir = ['swa' num2str(i-5) '/'];
        end
        fdr = [basedir filedir];
        pdiags = csvread([fdr 'Part-split-diagnostics(stage1).csv'],1);
        cdiags = csvread([fdr 'Chem-split-diagnostics(stage1).csv'],1);
        
        figure(10000)
        set(gcf,'color','white')
        semilogy(pdiags(:,1)*1000,pdiags(:,5),'color',cols(i,:),'LineStyle',lins{i})
        hold on
        if i == 17
            xlabel('Time (ms)')
            ylabel('Sample volume (-)')
            l = legend('D1','D2','D3','D4','D5',...
                'S1','S2','S3','S4','S5',...
                'S6','S7','S8','S9','S10','S11','S12');
            l.Location = 'EastOutside';
            saveas(['vsmp' studyid])
        end
        
        figure(10001)
        set(gcf,'color','white')
        plot(pdiags(:,1)*1000,pdiags(:,13),'color',cols(i,:),'LineStyle',lins{i})
        hold on
        if i == 17
            xlabel('Time (ms)')
            ylabel('Incepting weight (-)')
            l = legend('D1','D2','D3','D4','D5',...
                'S1','S2','S3','S4','S5',...
                'S6','S7','S8','S9','S10','S11','S12');
            l.Location = 'EastOutside';
            saveas(['inc_weight' studyid])
        end
    end
else
    %% Regular diagnostics
    
    % initialise saveas function if saving figures,
    % otherwise do nothing here
    if (savefigs && i == 17)
        disp('Saving figures - continuing may overwrite existing!')
        pause
        disp('Continuing...')
        saveas = @(str)(print([projpath projdir imagedir str],'-depsc'));
    else
        saveas = @(str)(disp('Not saving'));
    end
    
    % load data
    folder = [basedir filedir];
    pdiags = csvread([folder 'Part-split-diagnostics(stage1).csv'],1);
    cdiags = csvread([folder 'Chem-split-diagnostics(stage1).csv'],1);
    
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
    plot(pdiags(:,1)*1000,wnew(pdiags(:,6)))
    hold on
    plot(pdiags(:,1)*1000,wnew(pdiags(:,7)),':')
    plot([pdiags(1,1)*1000 pdiags(end,1)*1000],[wmax wmax],'r--')
    plot([pdiags(1,1)*1000 pdiags(end,1)*1000],[wmin wmin],'g--')
    plot(pdiags(:,1)*1000,pdiags(:,12),'--')
    plot(pdiags(:,1)*1000,pdiags(:,13),'-.')
    set(gca,'XLim',[0 pdiags(end,1)*1000])
    legend('Pre','Post','Wmax','Wmin','Pre sim','Post sim',...
        'location','North','orientation','vertical')
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
    
    figure(2)
    set(gcf,'color','white')
    [vals,bins] = hist(pdiags(:,15));
    plot(bins,vals/max(vals))
    hold on
    xlabel('Incepting factor (-)')
    ylabel('Count divide by max count (-)')
    
    %% Rates
    
    figure(3)
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
    
    figure(4)
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
        title(titl_d{i},'color',c)
        legend('Real','Fictitious','Location','NorthWest')
    end
    
    %% Rate fractions
    
    totevs = sum(pdiags(:,8:end-2),2);
    
    figure(6)
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
end