%% float 09965 
clear
close all

%% load ssh data
path = '/Users/lindsay.grose/Documents/MATLAB/Datasets/SSH.nosync/';
file = 'Cape_Basin_2023_03_04_2023_11_13.nc';
ssh_data = ncload([path,file]);

ogdate = datetime([2023,03,04,0,0,0],'timezone','UTC','format','yyyy-MM-dd HH:mm:ss');

trial =zeros(length(ssh_data.time));
ssh_data.datetime = trial(:,1);

for i = 1:length(ssh_data.time)
    ssh_data.datetime(i,1) = datenum(ogdate + days(i-1));
end

ssh_data.datetime = dateshift(datetime(ssh_data.datetime,'ConvertFrom','datenum'),'start','day');

mymap = brewermap([],'Dark2');

colorwf = mymap(3,:);
% colorwf = [0,0,0];
colorglid = mymap(5,:);
colorem = mymap(4,:);

%% load DSC data and depth data
load('/Users/lindsay.grose/Documents/MATLAB/Datasets/EM_Apex_Floats/Processes_data/Density/09965_str.mat')
load('/Users/lindsay.grose/Documents/MATLAB/Datasets/EM_Apex_Floats/Processes_data/Depth/09965_str.mat')

float_density_data_grid.datetime = dateshift(datetime(float_density_data_grid.time,'ConvertFrom','datenum'),'Start','day');
float_density_data_grid.potrho_anom = float_density_data_grid.potential_density - 1000;

%% calculate N2

[N2,pmid] = gsw_Nsquared(float_depth_data_grid.absolute_sal, float_depth_data_grid.conservative_temp, float_depth_data_grid.depth);

%% calc dynamic height

p_ref = 975;
dyn_height = gsw_geo_strf_dyn_height(float_depth_data_grid.absolute_sal, float_depth_data_grid.conservative_temp, float_depth_data_grid.depth, p_ref);


%integrated value
int_dyn_height_100_float = dyn_height(51,:)./9.81;

%% FSLE data

path = '/Users/lindsay.grose/Documents/MATLAB/Datasets/FSLE/FSLE5/';
file = 'combined_fsle.nc';
fsle_data_back = ncload([path,file]);

fsle_data_back.dt = datetime([2023,3,15]) + days(fsle_data_back.time);

file = 'forward/combined_fsle.nc';
fsle_data_for = ncload([path,file]);

fsle_data_for.dt = datetime([2023,3,15]) + days(fsle_data_for.time);

%get the time series
fsle_back = [];
fsle_for = [];

for i=1:length(float_density_data_grid.datetime(1,:))
    mylat = float_density_data_grid.lat(1,i);
    mylon = float_density_data_grid.lon(1,i);
    mytime = float_density_data_grid.datetime(1,i);

    latindnow = findClosest(mylat,fsle_data_back.lat);
    lonindnow = findClosest(mylon,fsle_data_back.lon);
    timeindnow = find(mytime == fsle_data_back.dt);

    fsle_back = [fsle_back,fsle_data_back.lambda1(latindnow,lonindnow,timeindnow)];
    fsle_for = [fsle_for,fsle_data_for.lambda1(latindnow,lonindnow,timeindnow)];

end

%% plot in depth space


fig = figure;
set(gcf,'Position',[72 1 1212 1315],'color','w')



ax1 = subplot(12,4,9:16);
ttt = pcolor(float_depth_data_grid.dist,float_depth_data_grid.depth,float_depth_data_grid.absolute_sal);
set(ttt,'EdgeColor','none')
set(ax1,'YDir','reverse');
colormap(ax1,cmocean('haline'))
cb = colorbar(ax1,'Position',[0.909486510008703,0.634980988593155,0.008703220191471,0.120912547528515]);
cb.FontSize = 20;
ax1.XAxis.FontSize = 24;
ax1.YAxis.FontSize = 24;
% xlabel('Distance (km)','FontWeight','Bold','FontSize',60)
ylabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',30)
ylabel(cb,'g/kg','FontWeight','Bold','FontSize',30)
title('Salinity','FontWeight','Bold','FontSize',30)
clim([34.3,35.8])
ylim(ax1,[0,1300])
hold on
line([mydists_pts(1),mydists_pts(1)],ylim,'Color',mycolors(1,:),'LineWidth',6)
line([mydists_pts(2),mydists_pts(2)],ylim,'Color',mycolors(2,:),'LineWidth',6)
line([mydists_pts(3),mydists_pts(3)],ylim,'Color',mycolors(3,:),'LineWidth',6)
line([mydists_pts(4),mydists_pts(4)],ylim,'Color',mycolors(4,:),'LineWidth',6)
pos = get(ax1, 'Position');
set(ax1, 'Position', pos -[0,0.03,0,0]);
xLimits = ax1.XLim;
yLimits = ax1.YLim;
text(xLimits(1)-150, yLimits(2)-1.6, 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
xlim(ax1,[0,maxdist])
xticklabels([]);

findClosest(1027,float_depth_data_grid.potrho(:,my_inds(2)))


float_depth_data_grid.depth(262,my_inds(2))


%% plot -- long

mydists_pts = [300,525,800,1100];
my_inds = findClosest(mydists_pts,float_density_data_grid.dist(1,:));

mycolors = [rgb('LightGreen');rgb('MediumSeaGreen');rgb('ForestGreen');rgb('DarkGreen')];
mycolors = [[247,247,247]./256;[204,204,204]./256;[150,150,150]./256;[82,82,82]./256];

subplot_labels = {'a)','b)','c)','d)'};

minlon = 12;
maxlon = 20;
minlat = -41;
maxlat = -34;

fig = figure;
set(gcf,'Position',[72 1 1212 1315],'color','w')

for i = 1:length(my_inds)
    ind = find(float_density_data_grid.datetime(1,my_inds(i)) == ssh_data.datetime);
    ax = subplot(12,4,[i,i+4]);
    ttt = pcolor(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,ind)');
    set(ttt,'EdgeColor','none')
    colormap(cmocean('rain'))
    hold on 
    contour(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,ind)',[1.5 1.25 1 0.75 0.5 0.1],'k',"ShowText",true,'LineWidth',5);
    plot(float_density_data_grid.lon(1,my_inds(i)-30:my_inds(i)+30),float_density_data_grid.lat(1,my_inds(i)-30:my_inds(i)+30),'color',colorem,'LineWidth',16)
    scatter(float_density_data_grid.lon(1,my_inds(i)),float_density_data_grid.lat(1,my_inds(i)),400,mycolors(i,:),'filled')
    set(ax,'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
    axis([minlon,maxlon,minlat,maxlat]);
    clim([-0.25,1.5])
    title(string(ssh_data.datetime(ind)),'FontWeight','Bold','FontSize',30)
    set(ax,'XTick',[14.0 18.0],'XTickLabel',{['14',char(176),'E'],['18',char(176),'E']})
    % set(ax,'XTick',[12.0 14.0 16.0 18.0 20.0 24.0],'XTickLabel',{['12',char(176),'E'],['16',char(176),'E'],['20',char(176),'E'],['24',char(176),'E']})
    set(ax,'YTick',[-40.0 -36.0 -32.0],'YTickLabel',{['40',char(176),'S'],['36',char(176),'S'],['32',char(176),'S']})
    ax.XAxis.FontSize = 26;
    ax.YAxis.FontSize = 26;
    ax.FontWeight = 'Bold'; 
    if i==length(my_inds)
        cb2 = colorbar(ax,'Position',[0.909486510008703,0.834980988593155,0.008703220191471,0.120152091254754]);
        cb2.FontSize = 20;
        ylabel(cb2, 'SSH (m)','FontWeight','Bold','FontSize',30) 
    end
    pos = get(ax, 'Position');
    % pos(3) = 0.15;  % Set the width of all subplots to be equal
    % pos(4) = 0.18;  % Set the height of all subplots to be equal
    set(ax, 'Position', pos +[-0.08+(i/40),+.03,0,0]);
    if i==1
        xLimits = ax.XLim;
        yLimits = ax.YLim;
        text(xLimits(1)-4, yLimits(2)+.5, string(subplot_labels{i}), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
    end 
end



maxdist = max(float_density_data_grid.dist(1,:));
minrho = 26.0;
maxrho= 27.5;


%calc N2 and DSC for density range
minrho = 1027.0;
maxrho = 1027.5;

ind1026 = findClosest(minrho,float_density_data_grid.potential_density);
ind10275 = findClosest(maxrho,float_density_data_grid.potential_density);

meanDSC_float = mean(abs(float_density_data_grid.DSC(ind1026:ind10275,:)),'omitmissing');


meanN2_float = nan(size(meanDSC_float));

for i=1:length(meanDSC_float)
    if length(find(~isnan(float_depth_data_grid.potrho(:,i)))) > 1
        ind1026 = findClosest(minrho,float_depth_data_grid.potrho(:,i));
        ind10275 = findClosest(maxrho,float_depth_data_grid.potrho(:,i));
        meanN2_float(i) = mean(real(sqrt(N2(ind1026:ind10275,i))),'omitmissing');
    end
end

ax8 = subplot(12,4,45:48);
plot(float_density_data_grid.dist(1,:),mirror_blackman_filter(meanN2_float',15)*1000,'LineWidth',3,'Color','Black')
pos = get(ax8, 'Position');
set(ax8, 'Position', pos+[0,-.05,0,0]);
ax8.FontSize = 24;
xlabel('Distance (km)','FontWeight','Bold','FontSize',30)
ylabel('N','FontWeight','Bold','FontSize',30)
xlim(ax8,[0,maxdist])
xLimits = ax8.XLim;
yLimits = ax8.YLim;
text(xLimits(1)-160, yLimits(2)+.5, 'i)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
ylim(ax8,[2,5])

ax7 = subplot(12,4,41:44);
plot(float_density_data_grid.dist(1,:),mirror_blackman_filter(meanDSC_float',15),'LineWidth',3,'Color','Black')
pos = get(ax7, 'Position');
set(ax7, 'Position', pos+[0,-.05,0,0]);
xticklabels([]);
ax7.FontSize = 24;
ylabel('|DSC|','FontWeight','Bold','FontSize',30)
xlim(ax7,[0,maxdist])
xLimits = ax7.XLim;
yLimits = ax7.YLim;
text(xLimits(1)-150, yLimits(2)-13, 'h)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
ylim(ax7,[0,80])

%calc N2 and DSC for density range
minrho = 1026.5;
maxrho = 1027.0;

ind1026 = findClosest(minrho,float_density_data_grid.potential_density);
ind10275 = findClosest(maxrho,float_density_data_grid.potential_density);

meanDSC_float = mean(abs(float_density_data_grid.DSC(ind1026:ind10275,:)),'omitmissing');


meanN2_float = nan(size(meanDSC_float));

for i=1:length(meanDSC_float)
    if length(find(~isnan(float_depth_data_grid.potrho(:,i)))) > 1
        ind1026 = findClosest(minrho,float_depth_data_grid.potrho(:,i));
        ind10275 = findClosest(maxrho,float_depth_data_grid.potrho(:,i));
        meanN2_float(i) = mean(real(sqrt(N2(ind1026:ind10275,i))),'omitmissing');
    end
end

ax6 = subplot(12,4,37:40);
plot(float_density_data_grid.dist(1,:),mirror_blackman_filter(meanN2_float',15)*1000,'LineWidth',3,'Color','Black')
pos = get(ax6, 'Position');
set(ax6, 'Position', pos+[0,-.05,0,0]);
ax6.FontSize = 24;
xticklabels([]);
% xlabel('Distance (km)','FontWeight','Bold','FontSize',30)
ylabel('N','FontWeight','Bold','FontSize',30)
xlim(ax6,[0,maxdist])
xLimits = ax6.XLim;
yLimits = ax6.YLim;
text(xLimits(1)-150, yLimits(2)-.3, 'g)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
ylim(ax6,[3.1,5.1])

ax5 = subplot(12,4,33:36);
plot(float_density_data_grid.dist(1,:),mirror_blackman_filter(meanDSC_float',15),'LineWidth',3,'Color','Black')
pos = get(ax5, 'Position');
set(ax5, 'Position', pos+[0,-.05,0,0]);
xticklabels([]);
ax5.FontSize = 24;
ylabel('|DSC|','FontWeight','Bold','FontSize',30)
xlim(ax5,[0,maxdist])
xLimits = ax5.XLim;
yLimits = ax5.YLim;
text(xLimits(1)-155, yLimits(2)+7, 'f)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
ylim(ax5,[0,80])



ax4 = subplot(12,4,29:32);
plot(float_density_data_grid.dist(1,:),fsle_back*-1,'LineWidth',3,'Color','Blue')
hold on
plot(float_density_data_grid.dist(1,:),fsle_for,'LineWidth',3,'Color','Red')
pos = get(ax4, 'Position');
set(ax4, 'Position', pos+[0,-.05,0,0]);
ax4.FontSize = 24;
xticklabels([]);
% xlabel('Distance (km)','FontWeight','Bold','FontSize',30)
ylabel('FSLE','FontWeight','Bold','FontSize',30)
xlim(ax4,[0,maxdist])
ylim(ax4,[0,2.5])
xLimits = ax4.XLim;
yLimits = ax4.YLim;
text(xLimits(1)-150, yLimits(2)-.6, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)


ax3 = subplot(12,4,25:28);
plot(float_density_data_grid.dist(1,:),int_dyn_height_100_float,'LineWidth',3,'Color','Black')
pos = get(ax3, 'Position');
set(ax3, 'Position', pos+[0,-0.05,0,0]);
ax3.FontSize = 24;
ylabel('\phi','FontWeight','Bold','FontSize',30)
xlim(ax3,[0,maxdist])
xticklabels([]);
xLimits = ax3.XLim;
yLimits = ax3.YLim;
text(xLimits(1)-150, yLimits(2)-.1, 'd)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)


minrho = 26.0;
maxrho= 27.5;

ax2 = subplot(12,4,17:24);
ttt = pcolor(float_density_data_grid.dist,float_density_data_grid.potrho_anom,float_density_data_grid.DSC);
set(ttt,'EdgeColor','none')
set(ax2,'YDir','reverse');
colormap(ax2,cmocean('balance'))
cb = colorbar(ax2,'Position',[0.909486510008703,0.476806083650189,0.008703220191471,0.120912547528515]);
cb.FontSize = 20;
ax2.XAxis.FontSize = 24;
ax2.YAxis.FontSize = 24;
% xlabel('Distance (km)','FontWeight','Bold','FontSize',60)
ylabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',30)
ylabel(cb,'m^3/kg','FontWeight','Bold','FontSize',30)
title('Spiciness Curvature','FontWeight','Bold','FontSize',30)
clim([-150,150])
ylim(ax2,[minrho,maxrho])
hold on
line([mydists_pts(1),mydists_pts(1)],ylim,'Color',mycolors(1,:),'LineWidth',6)
line([mydists_pts(2),mydists_pts(2)],ylim,'Color',mycolors(2,:),'LineWidth',6)
line([mydists_pts(3),mydists_pts(3)],ylim,'Color',mycolors(3,:),'LineWidth',6)
line([mydists_pts(4),mydists_pts(4)],ylim,'Color',mycolors(4,:),'LineWidth',6)
pos = get(ax2, 'Position');
set(ax2, 'Position', pos +[0,-0.05,0,0]);
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-150, yLimits(2)-1.6, 'c)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
xlim(ax2,[0,maxdist])
xticklabels([]);

ax1 = subplot(12,4,9:16);
ttt = pcolor(float_density_data_grid.dist,float_density_data_grid.potrho_anom,float_density_data_grid.absolute_salinity);
set(ttt,'EdgeColor','none')
set(ax1,'YDir','reverse');
colormap(ax1,cmocean('haline'))
cb = colorbar(ax1,'Position',[0.909486510008703,0.634980988593155,0.008703220191471,0.120912547528515]);
cb.FontSize = 20;
ax1.XAxis.FontSize = 24;
ax1.YAxis.FontSize = 24;
% xlabel('Distance (km)','FontWeight','Bold','FontSize',60)
ylabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',30)
ylabel(cb,'g/kg','FontWeight','Bold','FontSize',30)
title('Salinity','FontWeight','Bold','FontSize',30)
clim([34.3,35.8])
ylim(ax1,[minrho,maxrho])
hold on
line([mydists_pts(1),mydists_pts(1)],ylim,'Color',mycolors(1,:),'LineWidth',6)
line([mydists_pts(2),mydists_pts(2)],ylim,'Color',mycolors(2,:),'LineWidth',6)
line([mydists_pts(3),mydists_pts(3)],ylim,'Color',mycolors(3,:),'LineWidth',6)
line([mydists_pts(4),mydists_pts(4)],ylim,'Color',mycolors(4,:),'LineWidth',6)
pos = get(ax1, 'Position');
set(ax1, 'Position', pos -[0,0.03,0,0]);
xLimits = ax1.XLim;
yLimits = ax1.YLim;
text(xLimits(1)-150, yLimits(2)-1.6, 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
xlim(ax1,[0,maxdist])
xticklabels([]);

% saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/float09965_ssh_dsc.svg');

%% plot short -- changed - mutliaxis figure

mydists_pts = [300,525,800,1100];
my_inds = findClosest(mydists_pts,float_density_data_grid.dist(1,:));

mycolors = [rgb('LightGreen');rgb('MediumSeaGreen');rgb('ForestGreen');rgb('DarkGreen')];
mycolors = [[247,247,247]./256;[204,204,204]./256;[150,150,150]./256;[82,82,82]./256];

subplot_labels = {'a)','b)','c)','d)'};

minlon = 12;
maxlon = 20;
minlat = -41;
maxlat = -34;

%calc N2 and DSC for density range
minrho = 1026.5;
maxrho = 1027.0;

ind1026 = findClosest(minrho,float_density_data_grid.potential_density);
ind10275 = findClosest(maxrho,float_density_data_grid.potential_density);

meanDSC_float_t = mean(abs(float_density_data_grid.DSC(ind1026:ind10275,:)),'omitmissing');


meanN2_float_t = nan(size(meanDSC_float_t));

for i=1:length(meanDSC_float_t)
    if length(find(~isnan(float_depth_data_grid.potrho(:,i)))) > 1
        ind1026 = findClosest(minrho,float_depth_data_grid.potrho(:,i));
        ind10275 = findClosest(maxrho,float_depth_data_grid.potrho(:,i));
        meanN2_float_t(i) = mean(real(sqrt(N2(ind1026:ind10275,i))),'omitmissing');
        % meanN2_float_t(i) = mean(N2(ind1026:ind10275,i),'omitmissing');
    end
end

%calc N2 and DSC for density range
minrho = 1027.0;
maxrho = 1027.5;

ind1026 = findClosest(minrho,float_density_data_grid.potential_density);
ind10275 = findClosest(maxrho,float_density_data_grid.potential_density);

meanDSC_float_i = mean(abs(float_density_data_grid.DSC(ind1026:ind10275,:)),'omitmissing');


meanN2_float_i = nan(size(meanDSC_float_t));

for i=1:length(meanDSC_float_t)
    if length(find(~isnan(float_depth_data_grid.potrho(:,i)))) > 1
        ind1026 = findClosest(minrho,float_depth_data_grid.potrho(:,i));
        ind10275 = findClosest(maxrho,float_depth_data_grid.potrho(:,i));
        meanN2_float_i(i) = mean(real(sqrt(N2(ind1026:ind10275,i))),'omitmissing');
        % meanN2_float_i(i) = mean(N2(ind1026:ind10275,i),'omitmissing');
    end
end

fig = figure;
set(gcf,'Position',[72 1 1212 1315],'color','w')

for i = 1:length(my_inds)
    ind = find(float_density_data_grid.datetime(1,my_inds(i)) == ssh_data.datetime);
    ax = subplot(12,4,[i,i+4]);
    ttt = pcolor(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,ind)');
    set(ttt,'EdgeColor','none')
    colormap(cmocean('rain'))
    hold on 
    contour(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,ind)',[1.5 1.25 1 0.75 0.5 0.1],'k',"ShowText",true,'LineWidth',5);
    plot(float_density_data_grid.lon(1,my_inds(i)-30:my_inds(i)+30),float_density_data_grid.lat(1,my_inds(i)-30:my_inds(i)+30),'color',colorem,'LineWidth',16)
    scatter(float_density_data_grid.lon(1,my_inds(i)),float_density_data_grid.lat(1,my_inds(i)),400,mycolors(i,:),'filled')
    set(ax,'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
    axis([minlon,maxlon,minlat,maxlat]);
    clim([-0.25,1.5])
    title(string(ssh_data.datetime(ind)),'FontWeight','Bold','FontSize',30)
    set(ax,'XTick',[14.0 18.0],'XTickLabel',{['14',char(176),'E'],['18',char(176),'E']})
    % set(ax,'XTick',[12.0 14.0 16.0 18.0 20.0 24.0],'XTickLabel',{['12',char(176),'E'],['16',char(176),'E'],['20',char(176),'E'],['24',char(176),'E']})
    set(ax,'YTick',[-40.0 -36.0 -32.0],'YTickLabel',{['40',char(176),'S'],['36',char(176),'S'],['32',char(176),'S']})
    ax.XAxis.FontSize = 26;
    ax.YAxis.FontSize = 26;
    ax.FontWeight = 'Bold'; 
    if i==length(my_inds)
        cb2 = colorbar(ax,'Position',[0.909486510008703,0.834980988593155,0.008703220191471,0.120152091254754]);
        cb2.FontSize = 20;
        ylabel(cb2, 'SSH (m)','FontWeight','Bold','FontSize',30) 
    end
    pos = get(ax, 'Position');
    % pos(3) = 0.15;  % Set the width of all subplots to be equal
    % pos(4) = 0.18;  % Set the height of all subplots to be equal
    set(ax, 'Position', pos +[-0.08+(i/40),+.03,0,0]);
    xLimits = ax.XLim;
    yLimits = ax.YLim;
    if i==1
        text(xLimits(1)-4, yLimits(2)+.5, string(subplot_labels{i}), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
    end
    if i>=2
        text(xLimits(1)-2, yLimits(2)+.5, string(subplot_labels{i}), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
    end
end



maxdist = max(float_density_data_grid.dist(1,:));


ax5 = subplot(12,4,33:40);
hold on
% Left y-axis
yyaxis left
plot(float_density_data_grid.dist(1,:),mirror_blackman_filter(meanN2_float_t',15),'LineWidth',6,'Color',[230,171,2]./255)
hold on
plot(float_density_data_grid.dist(1,:),mirror_blackman_filter(meanN2_float_i',15),'--','LineWidth',3,'Color',[230,171,2]./255);
ylabel('N','FontWeight','Bold','FontSize',30)
ylim(ax5,[2.5e-3,6e-3])
ax5.YColor = [230,171,2]./255;
% ax5.YColor = [127,201,127]./255;

% Right y-axis
yyaxis right
plot(float_density_data_grid.dist(1,:),mirror_blackman_filter(meanDSC_float_t',15),'LineWidth',6,'Color',rgb('DarkSlateBlue'))
plot(float_density_data_grid.dist(1,:),mirror_blackman_filter(meanDSC_float_i',15), '--','LineWidth',3,'Color',rgb('DarkSlateBlue'));
ylabel('|DSC|','FontWeight','Bold','FontSize',30)
ylim(ax5,[0,74])
% ax5.YColor = [56,108,176]./255;   
ax5.YColor = rgb('DarkSlateBlue');   

pos = get(ax5, 'Position');
set(ax5, 'Position', pos+[0,-.1,0,0]);
xlabel('Distance (km)','FontWeight','Bold','FontSize',30)
ax5.FontSize = 24;
xlim(ax5,[0,maxdist])
xLimits = ax5.XLim;
yLimits = ax5.YLim;
text(xLimits(1)-155, yLimits(2)+4, 'i)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
% xlabel('Distance (km)','FontWeight','Bold','FontSize',30)


ax4 = subplot(12,4,29:32);
plot(float_density_data_grid.dist(1,:),fsle_back*-1,'LineWidth',3,'Color',[31,120,180]./255)
hold on
plot(float_density_data_grid.dist(1,:),fsle_for,'LineWidth',3,'Color',[227,26,28]./255)
pos = get(ax4, 'Position');
set(ax4, 'Position', pos+[0,-.07,0,0]);
ax4.FontSize = 24;
xticklabels([]);
% xlabel('Distance (km)','FontWeight','Bold','FontSize',30)
ylabel('FSLE','FontWeight','Bold','FontSize',30)
xlim(ax4,[0,maxdist])
ylim(ax4,[0,2.5])
xLimits = ax4.XLim;
yLimits = ax4.YLim;
text(xLimits(1)-150, yLimits(2)-.6, 'h)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)


ax3 = subplot(12,4,25:28);
plot(float_density_data_grid.dist(1,:),int_dyn_height_100_float,'LineWidth',3,'Color','Black')
pos = get(ax3, 'Position');
set(ax3, 'Position', pos+[0,-0.06,0,0]);
ax3.FontSize = 24;
ylabel('\phi','FontWeight','Bold','FontSize',30)
xlim(ax3,[0,maxdist])
xticklabels([]);
xLimits = ax3.XLim;
yLimits = ax3.YLim;
text(xLimits(1)-150, yLimits(2)-.1, 'g)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)


minrho = 26.0;
maxrho= 27.5;

ax2 = subplot(12,4,17:24);
ttt = pcolor(float_density_data_grid.dist,float_density_data_grid.potrho_anom,float_density_data_grid.DSC);
set(ttt,'EdgeColor','none')
set(ax2,'YDir','reverse');
colormap(ax2,cmocean('balance'))
cb = colorbar(ax2,'Position',[0.909486510008703,0.476806083650189,0.008703220191471,0.120912547528515]);
cb.FontSize = 20;
ax2.XAxis.FontSize = 24;
ax2.YAxis.FontSize = 24;
% xlabel('Distance (km)','FontWeight','Bold','FontSize',60)
ylabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',30)
ylabel(cb,'m^3/kg','FontWeight','Bold','FontSize',30)
title('Spiciness Curvature','FontWeight','Bold','FontSize',30)
clim([-150,150])
ylim(ax2,[minrho,maxrho])
hold on
line([mydists_pts(1),mydists_pts(1)],ylim,'Color',mycolors(1,:),'LineWidth',6)
line([mydists_pts(2),mydists_pts(2)],ylim,'Color',mycolors(2,:),'LineWidth',6)
line([mydists_pts(3),mydists_pts(3)],ylim,'Color',mycolors(3,:),'LineWidth',6)
line([mydists_pts(4),mydists_pts(4)],ylim,'Color',mycolors(4,:),'LineWidth',6)
pos = get(ax2, 'Position');
set(ax2, 'Position', pos +[0,-0.05,0,0]);
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-150, yLimits(2)-1.6, 'f)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
xlim(ax2,[0,maxdist])
xticklabels([]);

ax1 = subplot(12,4,9:16);
ttt = pcolor(float_density_data_grid.dist,float_density_data_grid.potrho_anom,float_density_data_grid.absolute_salinity);
set(ttt,'EdgeColor','none')
set(ax1,'YDir','reverse');
colormap(ax1,cmocean('haline'))
cb = colorbar(ax1,'Position',[0.909486510008703,0.634980988593155,0.008703220191471,0.120912547528515]);
cb.FontSize = 20;
ax1.XAxis.FontSize = 24;
ax1.YAxis.FontSize = 24;
% xlabel('Distance (km)','FontWeight','Bold','FontSize',60)
ylabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',30)
ylabel(cb,'g/kg','FontWeight','Bold','FontSize',30)
title('Salinity','FontWeight','Bold','FontSize',30)
clim([34.3,35.8])
ylim(ax1,[minrho,maxrho])
hold on
line([mydists_pts(1),mydists_pts(1)],ylim,'Color',mycolors(1,:),'LineWidth',6)
line([mydists_pts(2),mydists_pts(2)],ylim,'Color',mycolors(2,:),'LineWidth',6)
line([mydists_pts(3),mydists_pts(3)],ylim,'Color',mycolors(3,:),'LineWidth',6)
line([mydists_pts(4),mydists_pts(4)],ylim,'Color',mycolors(4,:),'LineWidth',6)
pos = get(ax1, 'Position');
set(ax1, 'Position', pos -[0,0.03,0,0]);
xLimits = ax1.XLim;
yLimits = ax1.YLim;
text(xLimits(1)-150, yLimits(2)-1.6, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',40)
xlim(ax1,[0,maxdist])
xticklabels([]);

saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/float09965_ssh_dsc2.svg');


%% plot of the T/S over rapid transition

mydists_pts = [480,530];
% mydists_pts = [420,460];
% mydists_pts = [515,550];
my_inds = findClosest(mydists_pts,float_density_data_grid.dist(1,:));

figure
set(gcf,'Position',[72 1 2395 1315],'color','w')

minrho = 26.0;
maxrho= 27.5;

ind = find(float_density_data_grid.datetime(1,floor(median(my_inds))) == ssh_data.datetime);
ax = subplot(2,4,2:3);
ttt = pcolor(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,ind)');
set(ttt,'EdgeColor','none')
colormap(cmocean('rain'))
hold on 
contour(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,ind)',[1.5 1.25 1 0.75 0.5 0.1],'k',"ShowText",true,'LineWidth',5);
plot(float_density_data_grid.lon(1,:),float_density_data_grid.lat(1,:),'color',colorem,'LineWidth',16)
scatter(float_density_data_grid.lon(1,my_inds(1):my_inds(2)),float_density_data_grid.lat(1,my_inds(1):my_inds(2)),400,'White','filled')
set(ax,'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
axis([minlon,maxlon,minlat,maxlat]);
clim([-0.25,1.5])
set(ax,'XTick',[14.0 18.0],'XTickLabel',{['14',char(176),'E'],['18',char(176),'E']})
% set(ax,'XTick',[12.0 14.0 16.0 18.0 20.0 24.0],'XTickLabel',{['12',char(176),'E'],['16',char(176),'E'],['20',char(176),'E'],['24',char(176),'E']})
set(ax,'YTick',[-40.0 -36.0 -32.0],'YTickLabel',{['40',char(176),'S'],['36',char(176),'S'],['32',char(176),'S']})
ax.XAxis.FontSize = 26;
ax.YAxis.FontSize = 26;
ax.FontWeight = 'Bold'; 



ax1 = subplot(2,4,5:8);
ttt = pcolor(float_density_data_grid.dist,float_density_data_grid.potrho_anom,float_density_data_grid.DSC);
set(ttt,'EdgeColor','none')
set(ax1,'YDir','reverse');
colormap(ax1,cmocean('balance'))
cb = colorbar;
cb.FontSize = 50;
ax1.XAxis.FontSize = 54;
ax1.YAxis.FontSize = 54;
xlabel('Distance (km)','FontWeight','Bold','FontSize',60)
ylabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',60)
ylabel(cb,'m^3/kg','FontWeight','Bold','FontSize',60)
title('Spiciness Curvature','FontWeight','Bold','FontSize',60)
clim([-150,150])
ylim(ax1,[minrho,maxrho])
hold on
line([mydists_pts(1),mydists_pts(1)],ylim,'Color',mycolors(1,:),'LineWidth',12)
line([mydists_pts(2),mydists_pts(2)],ylim,'Color',mycolors(2,:),'LineWidth',12)
pos = get(ax1, 'Position');
set(ax1, 'Position', pos +[0,0.02,0,0]);
xLimits = ax1.XLim;
yLimits = ax1.YLim;
text(xLimits(1)-90, yLimits(2)-1.6, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)

timeeee = datetime(float_density_data_grid.time(1,:),"ConvertFrom",'datenum');

cmap = cmocean('amp');              % Choose a colormap (e.g., 'jet', 'parula', etc.)
numColors = size(cmap, 1);
% Normalize z to map to colormap indices
zNormalized = (timeeee(my_inds(1):my_inds(2)) - min(timeeee(my_inds(1):my_inds(2)))) / (max(timeeee(my_inds(1):my_inds(2))) - min(timeeee(my_inds(1):my_inds(2))));  % Normalize to range [0, 1]
colorIndices = round(zNormalized * (numColors - 1)) + 1;

%%
load('pot_density_grid')
% Plot each line with corresponding color
figure
set(gcf,'Position',[684 204 1582 1133],'color','w')
ax4=subplot(1,1,1);
hold on;
for i = my_inds(2):-1:my_inds(1)
    color = cmap(colorIndices(i-(my_inds(1)-1)), :);  % Get the color for this line
    plot(float_density_data_grid.absolute_salinity(:,i), float_density_data_grid.conservative_temperature(:,i), 'Color', color, 'LineWidth', 5);
end
scatter(float_density_data_grid.absolute_salinity(:,i), float_density_data_grid.conservative_temperature(:,i)-100,50,float_density_data_grid.conservative_temperature(:,i)-100)
colormap(cmap)
cb1 = colorbar;
cb1.FontSize=30;
ylabel(cb1,'Hours since first profile')
clim([0,hours(max(timeeee(my_inds(1):my_inds(2))) - min(timeeee(my_inds(1):my_inds(2))))])
ax4.XAxis.FontSize = 44;
ax4.YAxis.FontSize = 44;
hold on
[contour_handle,~] = contour(pot_density_grid.sal,pot_density_grid.temp,pot_density_grid.pot_rho-1000,[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5],'k','HandleVisibility','off','color',rgb('DarkGray'));
labelhand = clabel(contour_handle,'FontSize',30);
for k = 1:length(labelhand)
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(27))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) + 0.5; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) + 2.5; % Shift left
    end
end
ylabel(['Temperature (',char(176),'C)'],'FontWeight','Bold','FontSize',80);
xlabel("Salinity (g/kg)",'FontWeight','Bold','FontSize',80);
% plot(meanprof_raw.SA_smooth,meanprof_raw.CT_smooth,'color',rgb('Cyan'),'linewidth',4)
ylim([2,22])
xlim([34.5,35.6])
% saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_gpn_vary2.svg');


%% plot with depth
figure
set(gcf,'Position',[684 204 1582 1133],'color','w')
ax4=subplot(1,1,1);
hold on;
for i = my_inds(2):-1:my_inds(1)
    color = cmap(colorIndices(i-(my_inds(1)-1)), :);  % Get the color for this line
    plot(float_density_data_grid.absolute_salinity(:,i), float_density_data_grid.depth(:,i), 'Color', color, 'LineWidth', 5);
end
colormap(cmap)
cb1 = colorbar;
cb1.FontSize=30;
ylabel(cb1,'Hours since first profile')
clim([0,hours(max(timeeee(my_inds(1):my_inds(2))) - min(timeeee(my_inds(1):my_inds(2))))])
ax4.XAxis.FontSize = 44;
ax4.YAxis.FontSize = 44;
hold on
ylabel(['Pressure'],'FontWeight','Bold','FontSize',80);
xlabel("Salinity (g/kg)",'FontWeight','Bold','FontSize',80);
% plot(meanprof_raw.SA_smooth,meanprof_raw.CT_smooth,'color',rgb('Cyan'),'linewidth',4)
ylim([0,1600])
xlim([34.5,35.6])
set(ax4,'YDir','Reverse')
% saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_gpn_vary2.svg');

figure
set(gcf,'Position',[684 204 1582 1133],'color','w')
ax4=subplot(1,1,1);
hold on;
for i = my_inds(2):-1:my_inds(1)
    color = cmap(colorIndices(i-(my_inds(1)-1)), :);  % Get the color for this line
    plot(float_density_data_grid.conservative_temperature(:,i), float_density_data_grid.depth(:,i), 'Color', color, 'LineWidth', 5);
end
colormap(cmap)
cb1 = colorbar;
cb1.FontSize=30;
ylabel(cb1,'Hours since first profile')
clim([0,hours(max(timeeee(my_inds(1):my_inds(2))) - min(timeeee(my_inds(1):my_inds(2))))])
ax4.XAxis.FontSize = 44;
ax4.YAxis.FontSize = 44;
hold on
ylabel(['Pressure'],'FontWeight','Bold','FontSize',80);
xlabel(['Temperature (',char(176),'C)'],'FontWeight','Bold','FontSize',80);
% plot(meanprof_raw.SA_smooth,meanprof_raw.CT_smooth,'color',rgb('Cyan'),'linewidth',4)
ylim([0,1600])
xlim([2,22])
set(ax4,'YDir','Reverse')
% saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_gpn_vary2.svg');

[N2,pmid] = gsw_Nsquared(float_depth_data_grid.absolute_sal,float_depth_data_grid.conservative_temp,float_depth_data_grid.depth);

figure
set(gcf,'Position',[684 204 1582 1133],'color','w')
ax4=subplot(1,1,1);
hold on;
for i = my_inds(2):-1:my_inds(1)
    color = cmap(colorIndices(i-(my_inds(1)-1)), :);  % Get the color for this line
    plot(movmean(N2(:,i),[5,5]), float_depth_data_grid.depth(1:end-1,i), 'Color', color, 'LineWidth', 5);
end
colormap(cmap)
cb1 = colorbar;
cb1.FontSize=30;
ylabel(cb1,'Hours since first profile')
clim([0,hours(max(timeeee(my_inds(1):my_inds(2))) - min(timeeee(my_inds(1):my_inds(2))))])
ax4.XAxis.FontSize = 44;
ax4.YAxis.FontSize = 44;
hold on
ylabel(['Pressure'],'FontWeight','Bold','FontSize',80);
xlabel(['N2'],'FontWeight','Bold','FontSize',80);
% plot(meanprof_raw.SA_smooth,meanprof_raw.CT_smooth,'color',rgb('Cyan'),'linewidth',4)
ylim([0,1000])
xlim([-1e-5,3e-4])
set(ax4,'YDir','Reverse')
% saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_gpn_vary2.svg');
