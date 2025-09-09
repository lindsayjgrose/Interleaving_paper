%% make a move of the cruise stuff
clear
close all

addpath(genpath("/Users/lindsay.grose/Documents/URI/git.nosync/"))
addpath(genpath("/Users/lindsay.grose/Documents/MATLAB/"))

D = dir;
cruise_path = '/Users/lindsay.grose/Documents/MATLAB/Datasets/EM_Apex_Floats/Processes_data/Depth/';
cd(cruise_path);

%%

%WF transects, floats, and glider tracks on top of SSH 

path = '/Users/lindsay.grose/Documents/MATLAB/Datasets/SSH.nosync/';
%/Volumes/MGL2312/raw/adcp/proc/os75bb/contour/os75bb.nc
file = 'Cape_Basin_2023_03_04_2023_11_13.nc';
ssh_data = ncload([path,file]);

ogdate = datetime([2023,03,04,0,0,0],'timezone','UTC','format','yyyy-MM-dd HH:mm:ss');

trial =zeros(length(ssh_data.time));
ssh_data.datetime = trial(:,1)

for i = 1:length(ssh_data.time)
    ssh_data.datetime(i,1) = datenum(ogdate + days(i-1))
end


%load WF timestamps and lat lon
load('wf_all.mat')

%%

datenow = dateshift(datetime(ssh_data.datetime,'ConvertFrom','datenum'),'start','day');
two_month_ind = find(datenow < dateshift(datetime([2023,04,30,0,0,0],'format','yyyy-MM-dd HH:mm:ss'),'start','day'));

alldata.WF_lat = [];
alldata.WF_lon = [];

for wf=1:10
    alldata.WF_lat = [alldata.WF_lat;wf_all.lat{wf}];
    alldata.WF_lon = [alldata.WF_lon;wf_all.lon{wf}];
end

alldata_deep.WF_lat = [];
alldata_deep.WF_lon = [];

deepnums = [5,6,7,8];
for wf=1:length(deepnums)
    alldata_deep.WF_lat = [alldata_deep.WF_lat;wf_all.lat{deepnums(wf)}];
    alldata_deep.WF_lon = [alldata_deep.WF_lon;wf_all.lon{deepnums(wf)}];
end

    
alldata.ssh_time = datenow;
alldata.ssh_lat = ssh_data.latitude;
alldata.ssh_lon = ssh_data.longitude;
alldata.ssh_adt = ssh_data.adt;


minlon = 10;
maxlon = 25;
minlat = -42;
maxlat = -30;

current_date_ssh = datenow((19));
sshind = find(current_date_ssh == alldata.ssh_time);

load('Agul_bath.mat')
bath_data.lat = bath_data.lat(1:3:end);
bath_data.lon = bath_data.lon(1:3:end);
bath_data.elevation = bath_data.elevation(1:3:end,1:3:end);
load('BrownToBlue.mat')

ind = find((bath_data.elevation*-1) < 0);
bath_data.elevation(ind) = NaN;

BrownToBlue(112:128,:) = [];

%%% need to be in float depth directory
pat = '(\d\d\d\d\d.*)';
dir_names = [];
for ii = 1:length(D) 
  a = regexp(D(ii).name,pat);
  if(a==1)
     dir_names = [dir_names;D(ii).name];
  end
end

alldata.float_lat = [];
alldata.float_lon = [];
for i = 1:length(dir_names)
    load([pwd,'/',dir_names(i,:)])

    section_name = dir_names(i,:);
    section_name = section_name(1:length(section_name)-4);
    if ~strcmp(section_name, '09969_all') && ~strcmp(section_name, '09970_all') && ~strcmp(section_name, '09965_end')
        alldata.float_lat = [alldata.float_lat,float_depth_data_grid.lat(1,:)];
        alldata.float_lon = [alldata.float_lon,float_depth_data_grid.lon(1,:)];
    end
end

%%

mymap = brewermap([],'Dark2');

% colorwf = mymap(3,:);
colorwf = [0,0,0];
colorglid = mymap(5,:);
colorem = mymap(4,:);

original_color = colorglid;  % Given color
darkening_factor = 0.77;  % Reduce brightness (0.8 = 20% darker)
darker_color = original_color * darkening_factor;

smoothed = smoothdata((bath_data.elevation')*-1, 1, 'gaussian', 25);  % smooth down columns
smoothed = smoothdata(smoothed, 2, 'gaussian', 25);  % then across rows


figure
ax1 = subplot(1,2,1);
set(gcf,'Position',[75 385 2322 952],'color','w')
ttt = pcolor(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,sshind)');
set(ttt,'EdgeColor','none')
colormap(cmocean('rain'))
auxh = colorbar();
hold on 
contour(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,sshind)',[1.5 1 0.75 0.25 0 -0.25],'k',"ShowText",true,'LineWidth',2);
contour(bath_data.lon,bath_data.lat,smoothed,[500,1000,2000,3000],'LineColor',rgb('Beige'),"ShowText",true,'LineWidth',3)
% plot(alldata.float_lon(floatind),alldata.float_lat(floatind),'.','color',mymap(4,:),'MarkerSize',40,'HandleVisibility','off')
plot(alldata.WF_lon,alldata.WF_lat,'.','color',rgb('RoyalBlue'),'MarkerSize',40,'HandleVisibility','off')
plot(alldata_deep.WF_lon,alldata_deep.WF_lat,'.','color',rgb('Red'),'MarkerSize',40,'HandleVisibility','off')
h_text = plot(nan,nan,'LineStyle','none');
b = plot(12,-50,'.','color',rgb('RoyalBlue'),'MarkerSize',70);
b = plot(12,-50,'.','color',rgb('Red'),'MarkerSize',70);
legend({'','','','Wire Flyer:','Shallow Dives','Deep Dives'},'FontWeight','Bold','FontSize',50)
set(ax1,'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
axis([minlon,maxlon,minlat,maxlat]);
clim([-.25,1.5])
ax1.FontWeight = 'Bold'; 
cb1 = colorbar(ax1,'Position',[0.413867355727821,0.112394957983193,0.02,0.810924369747898]);
cb1.FontSize = 32;
ylabel(cb1, 'SSH (m)','FontWeight','Bold','FontSize',60) 
set(ax1,'XTick',[12.0 16.0 20.0 24.0],'XTickLabel',{['12',char(176),'E'],['16',char(176),'E'],['20',char(176),'E'],['24',char(176),'E']})
set(ax1,'YTick',[-40.0 -36.0 -32.0],'YTickLabel',{['40',char(176),'S'],['36',char(176),'S'],['32',char(176),'S']})
ax1.XAxis.FontSize = 60;
ax1.YAxis.FontSize = 60;
ax1.XAxis.FontWeight = 'bold';
ax1.YAxis.FontWeight = 'bold';
pos = get(ax1, 'Position');
set(ax1, 'Position', pos +[-0.06,0,0,0]);
% lon1 = 20.1;
% lat1 = -38.6;
% plot(lon1,lat1,'p','MarkerSize',40,'MarkerFaceColor','Red','MarkerEdgeColor','red','HandleVisibility','off')
xLimits = ax1.XLim;
yLimits = ax1.YLim;
text(xLimits(1)-1.7, yLimits(2), 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
text(11.4, -37.2, 'C-A','color',rgb('Khaki'),'FontWeight','Bold','FontSize',60)
text(14.4, -38.8, 'SC','color',rgb('Khaki'),'FontWeight','Bold','FontSize',60)
text(16.4, -36, 'A','color',rgb('Khaki'),'FontWeight','Bold','FontSize',60)


ax2 = subplot(1,2,2);
ttt = pcolor(bath_data.lon,bath_data.lat,(bath_data.elevation')*-1);
set(ttt,'EdgeColor','none')
colormap(ax2,flipud(BrownToBlue))
axis(ax2,[minlon,maxlon,minlat,maxlat])
clim([0,6000])
hold on
plot(alldata.float_lon,alldata.float_lat,'.','color',rgb('Black'),'MarkerSize',20,'HandleVisibility','off')
t = plot(12,-50,'.','color',rgb('Black'),'MarkerSize',70);
legend({'', sprintf('Profiling\n  Float')}, 'FontWeight', 'Bold', 'FontSize', 50)
set(ax2,'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
ax2.FontWeight = 'Bold'; 
set(ax2,'XTick',[12.0 16.0 20.0 24.0],'XTickLabel',{['12',char(176),'E'],['16',char(176),'E'],['20',char(176),'E'],['24',char(176),'E']})
set(ax2,'YTick',[-40.0 -36.0 -32.0],'YTickLabel',{['40',char(176),'S'],['36',char(176),'S'],['32',char(176),'S']})
ax2.XAxis.FontSize = 60;
ax2.YAxis.FontSize = 60;
cb2 = colorbar(ax2,'Position',[0.912575366063738,0.112394957983193,0.02,0.810924369747898]);
cb2.FontSize = 32;
set(cb2, 'YDir', 'reverse' );
ylabel(cb2, 'Depth (m)','FontWeight','Bold','FontSize',60) 
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-1.7, yLimits(2), 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
yticklabels([])
text(12.4, -39, {'Schmidt-','Ott'},'color',rgb('White'),'FontWeight','Bold','FontSize',20,'HorizontalAlignment', 'center')
text(14.4, -38.5, 'Erica','color',rgb('White'),'FontWeight','Bold','FontSize',20)

saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/data_map.svg');

