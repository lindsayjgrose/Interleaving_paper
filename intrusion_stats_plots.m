clear 
close all
%% Calculate n^2 for Argo
addpath(genpath("/Users/lindsay.grose/Documents/URI/git.nosync/"))
addpath(genpath("/Users/lindsay.grose/Documents/MATLAB/"))


% em = load('all_CBstatistics_griddeddata.mat');
em = load('all_CBstatistics_griddeddata2.mat');
glider = load('intrusion_statistics_new.mat');
all_WF = load('all_WFstatistics_new.mat');

%% colors

colorwf = rgb('Red');
colorem = rgb('Black');


%% WF only
axisfontsize = 40;
labelfontsize = 50;

%%horizontal scales
figure
set(gcf,'Position',[-190 1 2700 1330],'color','w')
ax1=subplot(2,6,1:2);
histogram(all_WF.statistics.all.horizontal_scales,'BinWidth',1.8,'FaceColor',colorwf)
xlim([0,75])
ylim([0,110])
ax1.XAxis.FontSize = axisfontsize;
ax1.YAxis.FontSize = axisfontsize;
% ylim([0,75])
% set(gca,'YScale','Log')
xlabel('Length (km)','FontWeight','Bold','FontSize',labelfontsize)
ylabel('No. of Features','FontWeight','Bold','FontSize',labelfontsize)
% legend({'WF','EM-APEX','Glider'},'FontWeight','Bold','FontSize',30)
xLimits = ax1.XLim;
yLimits = ax1.YLim;
text(xLimits(1)-1, yLimits(2)+2, 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)

%%vertical scales
% figure
% set(gcf,'Position',[806 205 1460 1132],'color','w')
ax2=subplot(2,6,3:4);
histogram(all_WF.statistics.all.vertical_scales,'BinWidth',1.8,'FaceColor',colorwf)
xlim([0,80])
ylim([0,55])
ax2.XAxis.FontSize = axisfontsize;
ax2.YAxis.FontSize = axisfontsize;
% set(gca,'YScale','Log')
xlabel('Height (m)','FontWeight','Bold','FontSize',labelfontsize)
% ylabel('Count','FontWeight','Bold','FontSize',24)
% legend({'WF','EM-APEX','Glider'},'FontWeight','Bold','FontSize',20)
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-1, yLimits(2)+1, 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)

%%mean density
% figure
% set(gcf,'Position',[806 205 1460 1132],'color','w')
ax3=subplot(2,6,5:6);
histogram(all_WF.statistics.all.mean_density,'BinWidth',.05,'FaceColor',colorwf)
ax3.XAxis.FontSize = axisfontsize;
ax3.YAxis.FontSize = axisfontsize;
xlim([25.4,27.7])
ylim([0,55])
% set(gca,'YScale','Log')
xlabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',labelfontsize)
% xlabel('\textbf{Measured} \(\bar{\mathbf{x}}\)','interpreter','latex')
% ylabel('Count','FontWeight','Bold','FontSize',24)
% legend({'WF','EM-APEX','Glider'},'FontWeight','Bold','FontSize',20)
xLimits = ax3.XLim;
yLimits = ax3.YLim;
text(xLimits(1)-.05, yLimits(2)+1, 'c)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)
% Define shading limits (example: shade between 5 and 10 days)
x1 = 26;
x2 = 26.5;
yLimits = [0, 130]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray'), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Define shading limits (example: shade between 5 and 10 days)
x1 = 26.5;
x2 = 27;
yLimits = [0, 130]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray')-.2, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Define shading limits (example: shade between 5 and 10 days)
x1 = 27;
x2 = 27.5;
yLimits = [0, 130]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray')-.4, 'FaceAlpha', 0.3, 'EdgeColor', 'none');


meanN = 0.0034;
f = 2*7.2921e-5*sind(-37);

char_slope = f/meanN;

%%slope depth
% figure
% set(gcf,'Position',[806 205 1460 1132],'color','w')
ax5=subplot(2,6,8:9);
histogram(all_WF.statistics.all.slope_depth,'BinWidth',.0013,'FaceColor',colorwf)
ax5.XAxis.FontSize = axisfontsize;
ax5.YAxis.FontSize = axisfontsize;
xlim([-0.03,0.03])
ylim([0,115])
xline(char_slope,'color',rgb('RoyalBlue'),'LineWidth',5)
xline(-char_slope,'color',rgb('RoyalBlue'),'LineWidth',5)
% set(gca,'YScale','Log')
xlabel({'Slope (dz/dx)'},'FontWeight','Bold','FontSize',labelfontsize)
ylabel('No. of Features','FontWeight','Bold','FontSize',labelfontsize)
% ylabel('Count','FontWeight','Bold','FontSize',24)
% legend({'WF','EM-APEX','Glider'},'FontWeight','Bold','FontSize',20)
xLimits = ax5.XLim;
yLimits = ax5.YLim;
text(xLimits(1)-.001, yLimits(2)-.5, 'd)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)


%%slope density
% figure
% set(gcf,'Position',[806 205 1460 1132],'color','w')
ax6=subplot(2,6,10:11);
histogram(all_WF.statistics.all.slope_density*100000,'BinWidth',.0000015*100000,'FaceColor',colorwf)
xlim([-3e-5*100000,3e-5*100000])
ylim([0,55])
ax6.XAxis.FontSize = axisfontsize;
ax6.YAxis.FontSize = axisfontsize;
% set(gca,'YScale','Log')
xlabel({'Slope (d\rho/dx)'},'FontWeight','Bold','FontSize',labelfontsize)
% ylabel('Count','FontWeight','Bold','FontSize',24)
% legend({'WF','EM-APEX','Glider'},'FontWeight','Bold','FontSize',20)
% xtickformat('%.1e');  % Use scientific notation with 1 digit precision
xLimits = ax6.XLim;
yLimits = ax6.YLim;
text(xLimits(1)-.1e-5*100000, yLimits(2)+1, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)

saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/WF_statistics.svg');


%% time, length vs density - EM

cmap = brewermap([],'YlGn');
cmap(1:50,:) = [];

axisfontsize = 40;
labelfontsize = 50;


figure
set(gcf,'Position',[-9 -106 2504 1336],'color','w')


ax4=subplot(2,6,1:2);
% histogram(em.statistics.all.time,'BinWidth',.4,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',1)
histogram(em.statistics.all.mean_density-1000,'BinWidth',.05,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',1)
ax4.XAxis.FontSize = axisfontsize;
ax4.YAxis.FontSize = axisfontsize;
ylim([0,99])
% set(gca,'YScale','Log')
xlabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',labelfontsize)
ylabel('No. of Features','FontWeight','Bold','FontSize',labelfontsize)
% legend({'WF','EM-APEX','Glider'},'FontWeight','Bold','FontSize',20)
xlim([25,27.6])
xLimits = ax4.XLim;
yLimits = ax4.YLim;
text(xLimits(1)-.05, yLimits(2)+2, 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)
% Define shading limits (example: shade between 5 and 10 days)
x1 = 26;
x2 = 26.5;
yLimits = [0, 99]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray'), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Define shading limits (example: shade between 5 and 10 days)
x1 = 26.5;
x2 = 27;
yLimits = [0, 99]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray')-.2, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Define shading limits (example: shade between 5 and 10 days)
x1 = 27;
x2 = 27.5;
yLimits = [0, 99]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray')-.4, 'FaceAlpha', 0.3, 'EdgeColor', 'none');


%%vertical scales
ax1=subplot(2,6,3:4);
histogram(em.statistics.all.vertical_scales,'BinWidth',1.7,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',1)
xlim([0,90])
ylim([0,79])
ax1.XAxis.FontSize = axisfontsize;
ax1.YAxis.FontSize = axisfontsize;
% ylim([0,75])
% set(gca,'YScale','Log')
xlabel('Height (m)','FontWeight','Bold','FontSize',labelfontsize)
% ylabel('No. of Features','FontWeight','Bold','FontSize',labelfontsize)
xLimits = ax1.XLim;
yLimits = ax1.YLim;
text(xLimits(1)-1, yLimits(2)+1, 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)

%%time
ax2=subplot(2,6,5:6);
histogram(em.statistics.all.time,'BinWidth',0.4,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',1)
xlim([0,21])
ylim([0,129])
ax2.XAxis.FontSize = axisfontsize;
ax2.YAxis.FontSize = axisfontsize;
% set(gca,'YScale','Log')
xlabel('Lifetime (d)','FontWeight','Bold','FontSize',labelfontsize)
% ylabel('Count','FontWeight','Bold','FontSize',24)
% legend({'WF','EM-APEX','Glider'},'FontWeight','Bold','FontSize',20)
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-1, yLimits(2)+1, 'c)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)




ax=subplot(2,6,8:9);
scatter(em.statistics.all.mean_density-1000,em.statistics.all.time,130,rgb('DodgerBlue'),'filled')
% scatter(em.statistics.all.mean_density-1000,em.statistics.all.time,130,em.statistics.all.vertical_scales,'filled')
% colormap(ax,cmap)
% clim(ax,[0,45])
% hold on
% plot(glider.statistics.glider.mean_density-1000,glider.statistics.glider.time,'.','MarkerSize',50,'Color',colorglid)
ax.XAxis.FontSize = 40;
ax.YAxis.FontSize = 40;
% xlabel('Dist. from Retroflection','FontWeight','Bold','FontSize',50)
ylabel('Lifetime (d)','FontWeight','Bold','FontSize',50)
xlim([25.4,27.7])
ylim([0,24])
xlabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',50)
% title('EM only','FontWeight','Bold','FontSize',50)
xLimits = ax.XLim;
yLimits = ax.YLim;
text(xLimits(1)-.25, yLimits(2)-1, 'd)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)
pos = get(ax, 'Position');
set(ax, 'Position', pos+[-0.01,0,0,0]); 
% Define shading limits (example: shade between 5 and 10 days)
x1 = 26;
x2 = 26.5;
yLimits = [0, 99]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray'), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Define shading limits (example: shade between 5 and 10 days)
x1 = 26.5;
x2 = 27;
yLimits = [0, 99]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray')-.2, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Define shading limits (example: shade between 5 and 10 days)
x1 = 27;
x2 = 27.5;
yLimits = [0, 99]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray')-.4, 'FaceAlpha', 0.3, 'EdgeColor', 'none');


ax=subplot(2,6,10:11);
% scatter(em.statistics.all.mean_density-1000,em.statistics.all.horizontal_scales,130,rgb('DodgerBlue'),'filled')
% % scatter(em.statistics.all.mean_density-1000,em.statistics.all.horizontal_scales,130,em.statistics.all.vertical_scales,'filled')
% % colormap(ax,cmap)
% % clim(ax,[0,45])
% % hold on
% % plot(glider.statistics.glider.mean_density-1000,glider.statistics.glider.time,'.','MarkerSize',50,'Color',colorglid)
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% ylabel('Length (km)','FontWeight','Bold','FontSize',50)
% xlim([25.4,27.7])
% ylim([0,490])
% xlabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',50)
% % title('EM only','FontWeight','Bold','FontSize',50)
% xLimits = ax.XLim;
% yLimits = ax.YLim;
% text(xLimits(1)-.05, yLimits(2)+1, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)

% ax=subplot(2,3,6);
scatter(em.statistics.all.mean_density-1000,em.statistics.all.vertical_scales,130,rgb('DodgerBlue'),'filled');
% scatter(em.statistics.all.mean_density-1000,em.statistics.all.vertical_scales,130,em.statistics.all.vertical_scales,'filled')
% colormap(ax,cmap)
% clim(ax,[0,45])
% cb1 = colorbar;
% ylabel(cb1,'Vertical height','FontWeight','Bold','FontSize',40);
ax.XAxis.FontSize = 40;
ax.YAxis.FontSize = 40;
% xlabel('Dist. from Retroflection','FontWeight','Bold','FontSize',50)
ylabel('Height (m)','FontWeight','Bold','FontSize',50)
xlabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',50)
xlim([25.4,27.7])
ylim([0,124])
xLimits = ax.XLim;
yLimits = ax.YLim;
text(xLimits(1)-.35, yLimits(2)-5, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)
pos = get(ax, 'Position');
set(ax, 'Position', pos+[+0.03,0,0,0]); 
% Define shading limits (example: shade between 5 and 10 days)
x1 = 26;
x2 = 26.5;
yLimits = [0, 130]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray'), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Define shading limits (example: shade between 5 and 10 days)
x1 = 26.5;
x2 = 27;
yLimits = [0, 130]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray')-.2, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Define shading limits (example: shade between 5 and 10 days)
x1 = 27;
x2 = 27.5;
yLimits = [0, 130]; % or use ylim if you want dynamic adjustment
hold on
fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     rgb('LightGray')-.4, 'FaceAlpha', 0.3, 'EdgeColor', 'none');


saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/EM_statistics_length_height_density.svg');


