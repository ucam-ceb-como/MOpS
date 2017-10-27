clc, 
clear, close all

% Setup
% set defaults
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','Times')
set(0,'defaulttextfontname','Times')
set(0,'DefaultLineLineWidth',2.5);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10)

basedir = ''; %coag events - delete this to make space/

% for saving images
studyid  = '';
projpath = 'C:\Users\Astrid\Documents\Projects\';
projdir  = 'Network temperature dependence\Recent_work\';
imagedir = 'figures\coagEvents\'; 
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

% coagdat = csvread([basedir 'Coag-event-diagnostics-1micros.csv'],1);
% flagdat(:,1)=(coagdat(:,4)*1e9>100);
% flagdat(:,2)=(coagdat(:,5)*1e9<=49);
% 
% figure(100)
% set(gcf,'color','white')
% y1 = length(find(flagdat(:,1)==1))/length(flagdat(:,1));
% y2 = length(find(flagdat(:,2)==1))/length(flagdat(:,2));
% bar(1,y1,'FaceColor','y');
% hold on
% bar(2,y2,'FaceColor','g');
% y3 = length(find((flagdat(:,1)+flagdat(:,2))==2))/length(flagdat(:,1));
% bar(3,y3,'FaceColor','b');
% set(gca,'XTick',[1 2 3],'XTickLabel',{'Large','Small','Both'},'YLim',[0 1.2])
% ylabel('Fraction meeting size condition')
% xlabel('Particle')
% text(0.4,y1+0.05,'$d_{\textrm{c}}>100$ nm','FontSize',14)
% text(1.5,y2+0.05,'$d_{\textrm{c}}\leq49$ nm','FontSize',14)
% title('$t\in[0,100]$ $\mu$s')
% saveas(['fl_1micros' studyid])
% 
% figure(1)
% set(gcf,'color','white')
% h = histogram(coagdat(:,6),'NumBins',100);
% h.FaceColor = 'y';
% hold on;
% h = histogram(coagdat(:,7),'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Statistical weight (-)')
% title('$t\in[0,100]$ $\mu$s')
% saveas(['sw_1micros' studyid])
% 
% figure(2)
% set(gcf,'color','white')
% h = histogram(coagdat(:,4)*1e9,'NumBins',100);
% h.FaceColor = 'y';
% hold on
% h = histogram(coagdat(:,5)*1e9,'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Collision diameter (nm)')
% title('$t\in[0,100]$ $\mu$s')
% saveas(['dc_1micros' studyid])
% 
% figure(200)
% set(gcf,'color','white')
% dr = coagdat(:,4)./coagdat(:,5);
% mdr = mean(dr);
% h = histogram(dr,'NumBins',50);
% h.FaceColor = 'b';
% hold on
% ax = gca;
% set(ax,'YScale','log')
% plot([mdr mdr],[ax.YLim(1) ax.YLim(2)],':r')
% text(mdr,ax.YLim(2)*0.3,'mean','color','r','FontSize',14)
% ylabel('Count')
% xlabel('Collision diameter ratio (-)')
% title('$t\in[0,100]$ $\mu$s')
% saveas(['dr_1micros' studyid])
% 
% coagdat = csvread([basedir 'Coag-event-diagnostics-2micros.csv'],1);
% flagdat = [];
% flagdat(:,1)=(coagdat(:,4)*1e9>100);
% flagdat(:,2)=(coagdat(:,5)*1e9<=49);
% figure(300)
% set(gcf,'color','white')
% y1 = length(find(flagdat(:,1)==1))/length(flagdat(:,1));
% y2 = length(find(flagdat(:,2)==1))/length(flagdat(:,2));
% bar(1,y1,'FaceColor','y');
% hold on
% bar(2,y2,'FaceColor','g');
% y3 = length(find((flagdat(:,1)+flagdat(:,2))==2))/length(flagdat(:,1));
% bar(3,y3,'FaceColor','b');
% set(gca,'XTick',[1 2 3],'XTickLabel',{'Large','Small','Both'},'YLim',[0 1.2])
% ylabel('Fraction meeting size condition')
% xlabel('Particle')
% text(0.4,y1+0.05,'$d_{\textrm{c}}>100$ nm','FontSize',14)
% text(1.5,y2+0.05,'$d_{\textrm{c}}\leq49$ nm','FontSize',14)
% title('$t\in(100,200]$ $\mu$s')
% saveas(['fl_2micros' studyid])
% 
% figure(3)
% set(gcf,'color','white')
% h = histogram(coagdat(:,6),'NumBins',100);
% h.FaceColor = 'y';
% hold on;
% h = histogram(coagdat(:,7),'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Statistical weight (-)')
% title('$t\in(100,200]$ $\mu$s')
% saveas(['sw_2micros' studyid])
% 
% figure(4)
% set(gcf,'color','white')
% h = histogram(coagdat(:,4)*1e9,'NumBins',100);
% h.FaceColor = 'y';
% hold on
% h = histogram(coagdat(:,5)*1e9,'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Collision diameter (nm)')
% title('$t\in(100,200]$ $\mu$s')
% saveas(['dc_2micros' studyid])
% 
% figure(400)
% set(gcf,'color','white')
% dr = coagdat(:,4)./coagdat(:,5);
% mdr = mean(dr);
% h = histogram(dr,'NumBins',50);
% h.FaceColor = 'b';
% hold on
% ax = gca;
% set(ax,'YScale','log')
% plot([mdr mdr],[ax.YLim(1) ax.YLim(2)],':r')
% text(mdr,ax.YLim(2)*0.3,'mean','color','r','FontSize',14)
% ylabel('Count')
% xlabel('Collision diameter ratio (-)')
% title('$t\in(100,200]$ $\mu$s')
% saveas(['dr_2micros' studyid])
% 
% coagdat = csvread([basedir 'Coag-event-diagnostics-3micros.csv'],1);
% flagdat = [];
% flagdat(:,1)=(coagdat(:,4)*1e9>100);
% flagdat(:,2)=(coagdat(:,5)*1e9<=49);
% 
% figure(500)
% set(gcf,'color','white')
% y1 = length(find(flagdat(:,1)==1))/length(flagdat(:,1));
% y2 = length(find(flagdat(:,2)==1))/length(flagdat(:,2));
% bar(1,y1,'FaceColor','y');
% hold on
% bar(2,y2,'FaceColor','g');
% y3 = length(find((flagdat(:,1)+flagdat(:,2))==2))/length(flagdat(:,1));
% bar(3,y3,'FaceColor','b');
% set(gca,'XTick',[1 2 3],'XTickLabel',{'Large','Small','Both'},'YLim',[0 1.2])
% ylabel('Fraction meeting size condition')
% xlabel('Particle')
% text(0.4,y1+0.05,'$d_{\textrm{c}}>100$ nm','FontSize',14)
% text(1.5,y2+0.05,'$d_{\textrm{c}}\leq49$ nm','FontSize',14)
% title('$t\in(200,300]$ $\mu$s')
% saveas(['fl_3micros' studyid])
% 
% figure(5)
% set(gcf,'color','white')
% h = histogram(coagdat(:,6),'NumBins',100);
% h.FaceColor = 'y';
% hold on;
% h = histogram(coagdat(:,7),'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Statistical weight (-)')
% title('$t\in(200,300]$ $\mu$s')
% saveas(['sw_3micros' studyid])
% 
% figure(6)
% set(gcf,'color','white')
% h = histogram(coagdat(:,4)*1e9,'NumBins',100);
% h.FaceColor = 'y';
% hold on
% h = histogram(coagdat(:,5)*1e9,'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Collision diameter (nm)')
% title('$t\in(200,300]$ $\mu$s')
% saveas(['dc_3micros' studyid])
% 
% figure(600)
% set(gcf,'color','white')
% dr = coagdat(:,4)./coagdat(:,5);
% mdr = mean(dr);
% h = histogram(dr,'NumBins',50);
% h.FaceColor = 'b';
% hold on
% ax = gca;
% set(ax,'YScale','log')
% plot([mdr mdr],[ax.YLim(1) ax.YLim(2)],':r')
% text(mdr,ax.YLim(2)*0.3,'mean','color','r','FontSize',14)
% ylabel('Count')
% xlabel('Collision diameter ratio (-)')
% title('$t\in(200,300]$ $\mu$s')
% saveas(['dr_3micros' studyid])
% 
% coagdat = csvread([basedir 'Coag-event-diagnostics-4micros.csv'],1);
% flagdat = [];
% flagdat(:,1)=(coagdat(:,4)*1e9>100);
% flagdat(:,2)=(coagdat(:,5)*1e9<=49);
% 
% figure(700)
% set(gcf,'color','white')
% y1 = length(find(flagdat(:,1)==1))/length(flagdat(:,1));
% y2 = length(find(flagdat(:,2)==1))/length(flagdat(:,2));
% bar(1,y1,'FaceColor','y');
% hold on
% bar(2,y2,'FaceColor','g');
% y3 = length(find((flagdat(:,1)+flagdat(:,2))==2))/length(flagdat(:,1));
% bar(3,y3,'FaceColor','b');
% set(gca,'XTick',[1 2 3],'XTickLabel',{'Large','Small','Both'},'YLim',[0 1.2])
% ylabel('Fraction meeting size condition')
% xlabel('Particle')
% text(0.4,y1+0.05,'$d_{\textrm{c}}>100$ nm','FontSize',14)
% text(1.5,y2+0.05,'$d_{\textrm{c}}\leq49$ nm','FontSize',14)
% title('$t\in(300,400]$ $\mu$s')
% saveas(['fl_4micros' studyid])
% 
% figure(7)
% set(gcf,'color','white')
% h = histogram(coagdat(:,6),'NumBins',100);
% h.FaceColor = 'y';
% hold on;
% h = histogram(coagdat(:,7),'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Statistical weight (-)')
% title('$t\in(300,400]$ $\mu$s')
% saveas(['sw_4micros' studyid])
% 
% figure(8)
% set(gcf,'color','white')
% h = histogram(coagdat(:,4)*1e9,'NumBins',100);
% h.FaceColor = 'y';
% hold on
% h = histogram(coagdat(:,5)*1e9,'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Collision diameter (nm)')
% title('$t\in(300,400]$ $\mu$s')
% saveas(['dc_4micros' studyid])
% 
% figure(800)
% set(gcf,'color','white')
% dr = coagdat(:,4)./coagdat(:,5);
% mdr = mean(dr);
% h = histogram(dr,'NumBins',50);
% h.FaceColor = 'b';
% hold on
% ax = gca;
% set(ax,'YScale','log')
% plot([mdr mdr],[ax.YLim(1) ax.YLim(2)],':r')
% text(mdr,ax.YLim(2)*0.3,'mean','color','r','FontSize',14)
% ylabel('Count')
% xlabel('Collision diameter ratio (-)')
% title('$t\in(300,400]$ $\mu$s')
% saveas(['dr_4micros' studyid])
% 
% coagdat = csvread([basedir 'Coag-event-diagnostics-5-51micros.csv'],1);
% flagdat = [];
% flagdat(:,1)=(coagdat(:,4)*1e9>100);
% flagdat(:,2)=(coagdat(:,5)*1e9<=49);
% 
% figure(900)
% set(gcf,'color','white')
% y1 = length(find(flagdat(:,1)==1))/length(flagdat(:,1));
% y2 = length(find(flagdat(:,2)==1))/length(flagdat(:,2));
% bar(1,y1,'FaceColor','y');
% hold on
% bar(2,y2,'FaceColor','g');
% y3 = length(find((flagdat(:,1)+flagdat(:,2))==2))/length(flagdat(:,1));
% bar(3,y3,'FaceColor','b');
% set(gca,'XTick',[1 2 3],'XTickLabel',{'Large','Small','Both'},'YLim',[0 1.2])
% ylabel('Fraction meeting size condition')
% xlabel('Particle')
% text(0.4,y1+0.05,'$d_{\textrm{c}}>100$ nm','FontSize',14)
% text(1.5,y2+0.05,'$d_{\textrm{c}}\leq49$ nm','FontSize',14)
% title('$t\in(500,510]$ $\mu$s')
% saveas(['fl_5_1micros' studyid])
% 
% figure(9)
% set(gcf,'color','white')
% h = histogram(coagdat(:,6),'NumBins',100);
% h.FaceColor = 'y';
% hold on;
% h = histogram(coagdat(:,7),'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Statistical weight (-)')
% title('$t\in(500,510]$ $\mu$s')
% saveas(['sw_5_1micros' studyid])
% 
% figure(10)
% set(gcf,'color','white')
% h = histogram(coagdat(:,4)*1e9,'NumBins',100);
% h.FaceColor = 'y';
% hold on
% h = histogram(coagdat(:,5)*1e9,'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Collision diameter (nm)')
% title('$t\in(500,510]$ $\mu$s')
% saveas(['dc_5_1micros' studyid])
% 
% figure(1000)
% set(gcf,'color','white')
% dr = coagdat(:,4)./coagdat(:,5);
% mdr = mean(dr);
% h = histogram(dr,'NumBins',50);
% h.FaceColor = 'b';
% hold on
% ax = gca;
% set(ax,'YScale','log')
% plot([mdr mdr],[ax.YLim(1) ax.YLim(2)],':r')
% text(mdr,ax.YLim(2)*0.3,'mean','color','r','FontSize',14)
% ylabel('Count')
% xlabel('Collision diameter ratio (-)')
% title('$t\in(500,510]$ $\mu$s')
% saveas(['dr_5_1micros' studyid])
% 
% coagdat = csvread([basedir 'Coag-event-diagnostics-1-1ms.csv'],1);
% flagdat = [];
% flagdat(:,1)=(coagdat(:,4)*1e9>100);
% flagdat(:,2)=(coagdat(:,5)*1e9<=49);
% 
% figure(1100)
% set(gcf,'color','white')
% y1 = length(find(flagdat(:,1)==1))/length(flagdat(:,1));
% y2 = length(find(flagdat(:,2)==1))/length(flagdat(:,2));
% bar(1,y1,'FaceColor','y');
% hold on
% bar(2,y2,'FaceColor','g');
% y3 = length(find((flagdat(:,1)+flagdat(:,2))==2))/length(flagdat(:,1));
% bar(3,y3,'FaceColor','b');
% set(gca,'XTick',[1 2 3],'XTickLabel',{'Large','Small','Both'},'YLim',[0 1.2])
% ylabel('Fraction meeting size condition')
% xlabel('Particle')
% text(0.4,y1+0.05,'$d_{\textrm{c}}>100$ nm','FontSize',14)
% text(1.5,y2+0.05,'$d_{\textrm{c}}\leq49$ nm','FontSize',14)
% title('$t\in(1,1.1]$ ms')
% saveas(['fl_1_1ms' studyid])
% 
% figure(11)
% set(gcf,'color','white')
% h = histogram(coagdat(:,6),'NumBins',100);
% h.FaceColor = 'y';
% hold on;
% h = histogram(coagdat(:,7),'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Statistical weight (-)')
% title('$t\in(1,1.1]$ ms')
% saveas(['sw_1_1ms' studyid])
% 
% figure(12)
% set(gcf,'color','white')
% h = histogram(coagdat(:,4)*1e9,'NumBins',100);
% h.FaceColor = 'y';
% hold on
% h = histogram(coagdat(:,5)*1e9,'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Collision diameter (nm)')
% title('$t\in(1,1.1]$ ms')
% saveas(['dc_1_1ms' studyid])
% 
% figure(1200)
% set(gcf,'color','white')
% dr = coagdat(:,4)./coagdat(:,5);
% mdr = mean(dr);
% h = histogram(dr,'NumBins',50);
% h.FaceColor = 'b';
% hold on
% ax = gca;
% set(ax,'YScale','log')
% plot([mdr mdr],[ax.YLim(1) ax.YLim(2)],':r')
% text(mdr,ax.YLim(2)*0.3,'mean','color','r','FontSize',14)
% ylabel('Count')
% xlabel('Collision diameter ratio (-)')
% title('$t\in(1,1.1]$ $\mu$s')
% saveas(['dr_1_1ms' studyid])
% 
% coagdat = csvread([basedir 'Coag-event-diagnostics-2-1ms.csv'],1);
% flagdat = [];
% flagdat(:,1)=(coagdat(:,4)*1e9>100);
% flagdat(:,2)=(coagdat(:,5)*1e9<=49);
% 
% figure(1300)
% set(gcf,'color','white')
% y1 = length(find(flagdat(:,1)==1))/length(flagdat(:,1));
% y2 = length(find(flagdat(:,2)==1))/length(flagdat(:,2));
% bar(1,y1,'FaceColor','y');
% hold on
% bar(2,y2,'FaceColor','g');
% y3 = length(find((flagdat(:,1)+flagdat(:,2))==2))/length(flagdat(:,1));
% bar(3,y3,'FaceColor','b');
% set(gca,'XTick',[1 2 3],'XTickLabel',{'Large','Small','Both'},'YLim',[0 1.2])
% ylabel('Fraction meeting size condition')
% xlabel('Particle')
% text(0.4,y1+0.05,'$d_{\textrm{c}}>100$ nm','FontSize',14)
% text(1.5,y2+0.05,'$d_{\textrm{c}}\leq49$ nm','FontSize',14)
% title('$t\in(2,2.1]$ ms')
% saveas(['fl_2_1ms' studyid])
% 
% figure(13)
% set(gcf,'color','white')
% h = histogram(coagdat(:,6),'NumBins',100);
% h.FaceColor = 'y';
% hold on;
% h = histogram(coagdat(:,7),'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Statistical weight (-)')
% title('$t\in(2,2.1]$ ms')
% saveas(['sw_2_1ms' studyid])
% 
% figure(14)
% set(gcf,'color','white')
% h = histogram(coagdat(:,4)*1e9,'NumBins',100);
% h.FaceColor = 'y';
% hold on
% h = histogram(coagdat(:,5)*1e9,'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Collision diameter (nm)')
% title('$t\in(2,2.1]$ ms')
% saveas(['dc_2_1ms' studyid])
% 
% figure(1400)
% set(gcf,'color','white')
% dr = coagdat(:,4)./coagdat(:,5);
% mdr = mean(dr);
% h = histogram(dr,'NumBins',50);
% h.FaceColor = 'b';
% hold on
% ax = gca;
% set(ax,'YScale','log')
% plot([mdr mdr],[ax.YLim(1) ax.YLim(2)],':r')
% text(mdr,ax.YLim(2)*0.3,'mean','color','r','FontSize',14)
% ylabel('Count')
% xlabel('Collision diameter ratio (-)')
% title('$t\in(2,2.1]$ $\mu$s')
% saveas(['dr_2_1ms' studyid])
% 
% coagdat = csvread([basedir 'Coag-event-diagnostics-3-1ms.csv'],1);
% flagdat = [];
% flagdat(:,1)=(coagdat(:,4)*1e9>100);
% flagdat(:,2)=(coagdat(:,5)*1e9<=49);
% 
% figure(1500)
% set(gcf,'color','white')
% y1 = length(find(flagdat(:,1)==1))/length(flagdat(:,1));
% y2 = length(find(flagdat(:,2)==1))/length(flagdat(:,2));
% bar(1,y1,'FaceColor','y');
% hold on
% bar(2,y2,'FaceColor','g');
% y3 = length(find((flagdat(:,1)+flagdat(:,2))==2))/length(flagdat(:,1));
% bar(3,y3,'FaceColor','b');
% set(gca,'XTick',[1 2 3],'XTickLabel',{'Large','Small','Both'},'YLim',[0 1.2])
% ylabel('Fraction meeting size condition')
% xlabel('Particle')
% text(0.4,y1+0.05,'$d_{\textrm{c}}>100$ nm','FontSize',14)
% text(1.5,y2+0.05,'$d_{\textrm{c}}\leq49$ nm','FontSize',14)
% title('$t\in(3,3.1]$ ms')
% saveas(['fl_3_1ms' studyid])
% 
% figure(15)
% set(gcf,'color','white')
% h = histogram(coagdat(:,6),'NumBins',100);
% h.FaceColor = 'y';
% hold on;
% h = histogram(coagdat(:,7),'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Statistical weight (-)')
% title('$t\in(3,3.1]$ ms')
% saveas(['sw_3_1ms' studyid])
% 
% figure(16)
% set(gcf,'color','white')
% h = histogram(coagdat(:,4)*1e9,'NumBins',100);
% h.FaceColor = 'y';
% hold on
% h = histogram(coagdat(:,5)*1e9,'NumBins',100);
% h.FaceColor = 'g';
% set(gca,'YScale','log')
% legend('Large particle','Small particle','Location','North')
% ylabel('Count')
% xlabel('Collision diameter (nm)')
% title('$t\in(3,3.1]$ ms')
% saveas(['dc_3_1ms' studyid])
% 
% figure(1600)
% set(gcf,'color','white')
% dr = coagdat(:,4)./coagdat(:,5);
% mdr = mean(dr);
% h = histogram(dr,'NumBins',50);
% h.FaceColor = 'b';
% hold on
% ax = gca;
% set(ax,'YScale','log')
% plot([mdr mdr],[ax.YLim(1) ax.YLim(2)],':r')
% text(mdr,ax.YLim(2)*0.7,'mean','color','r','FontSize',14)
% ylabel('Count')
% xlabel('Collision diameter ratio (-)')
% title('$t\in(3,3.1]$ $\mu$s')
% saveas(['dr_3_1ms' studyid])

coagdat = csvread([basedir 'Coag-event-diagnostics-3micros.csv'],1);
pscdat = csvread([basedir 'PSC-event-diagnostics-3micros.csv'],1);
flagdat(:,1)=(coagdat(:,4)*1e9>100);
flagdat(:,2)=(coagdat(:,5)*1e9<=49);

figure(1)
set(gcf,'color','white')
y1 = length(find(flagdat(:,1)==1))/length(flagdat(:,1));
y2 = length(find(flagdat(:,2)==1))/length(flagdat(:,2));
bar(1,y1,'FaceColor','y');
hold on
bar(2,y2,'FaceColor','g');
y3 = length(find((flagdat(:,1)+flagdat(:,2))==2))/length(flagdat(:,1));
bar(3,y3,'FaceColor','b');
set(gca,'XTick',[1 2 3],'XTickLabel',{'Large','Small','Both'},'YLim',[0 1.2])
ylabel('Fraction meeting size condition')
xlabel('Particle')
text(0.4,y1+0.05,'$d_{\textrm{c}}>100$ nm','FontSize',14)
text(1.5,y2+0.05,'$d_{\textrm{c}}\leq49$ nm','FontSize',14)
title('$t\in(300,310]$ $\mu$s')
saveas(['fl_3micros_pcs' studyid])

figure(2)
set(gcf,'color','white')
h = histogram(coagdat(:,6),'NumBins',100);
h.FaceColor = 'y';
hold on;
h = histogram(coagdat(:,7),'NumBins',100);
h.FaceColor = 'g';
set(gca,'YScale','log')
legend('Large particle','Small particle','Location','North')
ylabel('Count')
xlabel('Statistical weight (-)')
title('$t\in(300,310]$ $\mu$s')
saveas(['sw_3micros_pcs' studyid])

figure(3)
set(gcf,'color','white')
h = histogram(coagdat(:,4)*1e9,'NumBins',100);
h.FaceColor = 'y';
hold on
h = histogram(coagdat(:,5)*1e9,'NumBins',100);
h.FaceColor = 'g';
set(gca,'YScale','log')
legend('Large particle','Small particle','Location','North')
ylabel('Count')
xlabel('Collision diameter (nm)')
title('$t\in(300,310]$ $\mu$s')
saveas(['dc_3micros_pcs' studyid])

figure(4)
set(gcf,'color','white')
dr = coagdat(:,4)./coagdat(:,5);
mdr = mean(dr);
h = histogram(dr,'NumBins',50);
h.FaceColor = 'b';
hold on
ax = gca;
set(ax,'YScale','log')
plot([mdr mdr],[ax.YLim(1) ax.YLim(2)],':r')
text(mdr,ax.YLim(2)*0.3,'mean','color','r','FontSize',14)
ylabel('Count')
xlabel('Collision diameter ratio (-)')
title('$t\in(300,310]$ $\mu$s')
saveas(['dr_3micros_pcs' studyid])

figure(5)
set(gcf,'color','white')
deltac1 = (coagdat(:,8)-coagdat(:,4))./coagdat(:,4);
deltac2 = (coagdat(:,9)-coagdat(:,5))./coagdat(:,5);
deltaw1 = (coagdat(:,6)-coagdat(:,10))./coagdat(:,6);
deltaw2 = (coagdat(:,7)-coagdat(:,11))./coagdat(:,7);
deltaca = (coagdat(:,13)-coagdat(:,12))./coagdat(:,12);
% dca_0 = coagdat(:,12);
% dca_1 = coagdat(:,13);
% t = coagdat(:,1);
subplot(221)
h1 = histogram(deltac1,'NumBins',10);
hold on
ax = gca;
set(ax,'YScale','log')
subplot(222)
h2 = histogram(deltac2,'NumBins',10);
hold on
ax = gca;
set(ax,'YScale','log')
subplot(223)
h3 = histogram(deltaw1,'NumBins',10);
hold on
ax = gca;
set(ax,'YScale','log')
subplot(224)
h4 = histogram(deltaw2,'NumBins',10);
hold on
ax = gca;
set(ax,'YScale','log')

figure(6) 
set(gcf,'color','white')
h = histogram(pscdat(:,2)*1e9,'NumBins',100);
h.FaceColor = 'y';
hold on
set(gca,'YScale','log')
ylabel('Count')
xlabel('Collision diameter (nm)')
title('$t\in(300,310]$ $\mu$s')
saveas(['pdc_3micros_pcs' studyid])

figure(7) 
set(gcf,'color','white')
h = histogram(pscdat(:,4),'NumBins',100);
h.FaceColor = 'y';
hold on
set(gca,'YScale','log')
ylabel('Count')
xlabel('Statistical weight (-)')
title('$t\in(300,310]$ $\mu$s')
saveas(['psw_3micros_pcs' studyid])

figure(8)
set(gcf,'color','white')
plot(pscdat(:,1),(pscdat(:,5)-pscdat(:,4))./pscdat(:,4))