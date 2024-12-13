%% make a move of the cruise stuff
clear
close all

addpath(genpath("/Users/lindsay.grose/Documents/URI/git.nosync/"))
addpath(genpath("/Users/lindsay.grose/Documents/MATLAB/"))
%%

%WF transects, floats, and glider tracks on top of SSH 


path = '/Users/lindsay.grose/Documents/MATLAB/Datasets/SSH.nosync/';
%/Volumes/MGL2312/raw/adcp/proc/os75bb/contour/os75bb.nc
file = 'Cape_Basin_2023-03-04_2023-11-13.nc';
ssh_data = ncload([path,file]);

ogdate = datetime([2023,03,04,0,0,0],'timezone','UTC','format','yyyy-MM-dd HH:mm:ss');

trial =zeros(length(ssh_data.time));
ssh_data.datetime = trial(:,1)

for i = 1:length(ssh_data.time)
    ssh_data.datetime(i,1) = datenum(ogdate + days(i-1))
end

%load adcp
load('ADCP_data_JH_edited.mat')

%load WF timestamps and lat lon
load('wf_all.mat')

%load EM-APEX timestamps and lat lon
load('FloatsCPMarch.mat')

datenow = dateshift(datetime(ssh_data.datetime,'ConvertFrom','datenum'),'start','day');
two_month_ind = find(datenow < dateshift(datetime([2023,04,30,0,0,0],'format','yyyy-MM-dd HH:mm:ss'),'start','day'));

alldata.WF_time = [];
alldata.WF_lat = [];
alldata.WF_lon = [];

for wf=1:10
    wftimes = wf_all.time{wf};
    wftimes = dateshift(unixtime2datetime(double(wftimes)/1e6),'start','day');
    alldata.WF_time = [alldata.WF_time;wftimes];
    alldata.WF_lat = [alldata.WF_lat;wf_all.lat{wf}];
    alldata.WF_lon = [alldata.WF_lon;wf_all.lon{wf}];
end

% get glider data in here
glider_data = ncload('/Users/lindsay.grose/Documents/MATLAB/Datasets/Glider/ds_ATD.nc');
glider_data.time = glider_data.time + 80.35674769; %% this is in julian day

glider_data.datetime = datetime(convertTo(datetime(2023,1,1,0,0,0),'juliandate') + glider_data.time,'convertfrom','Juliandate');

alldata.glider_time = dateshift(glider_data.datetime,'start','day');
alldata.glider_lat = glider_data.lat;
alldata.glider_lon = glider_data.long;

alldata.float_time = dateshift(datetime(gr.mt,'ConvertFrom','datenum'),'start','day');
alldata.float_lat = gr.lat;
alldata.float_lon = gr.lon;

alldata.adcp_time = dateshift(unixtime2datetime(double(ADCP_data.timestamp)/1e6),'start','day');
alldata.adcp_lat = ADCP_data.lat;
alldata.adcp_lon = ADCP_data.lon;
    
alldata.ssh_time = datenow;
alldata.ssh_lat = ssh_data.latitude;
alldata.ssh_lon = ssh_data.longitude;
alldata.ssh_adt = ssh_data.adt;


minlon = 10;
maxlon = 25;
minlat = -42;
maxlat = -30;

current_date = datenow((85));
current_date_ssh = datenow((19));
sshind = find(current_date_ssh == alldata.ssh_time);
floatind = find(current_date >= alldata.float_time);
wfind = find(current_date >= alldata.WF_time);
adcpind = find(current_date >= alldata.adcp_time);
gliderind = find(current_date >= alldata.glider_time);

%%
figure
ax1 = axes;
set(gcf,'Position',[88 52 1640 1123],'color','w')
ttt = pcolor(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,sshind)');
set(ttt,'EdgeColor','none')
colormap(cmocean('rain'))
auxh = colorbar();
hold on 
contour(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,sshind)',[1.5 1 0.75 0.25 0 -0.25],'k',"ShowText",true,'LineWidth',2);
% plot(alldata.adcp_lon(adcpind),alldata.adcp_lat(adcpind),'color',rgb('Black'),'linewidth',3,'HandleVisibility','off')
plot(alldata.float_lon(floatind),alldata.float_lat(floatind),'.','color',rgb('Red'),'MarkerSize',40,'HandleVisibility','off')
plot(alldata.WF_lon(wfind),alldata.WF_lat(wfind),'.','color',rgb('Blue'),'MarkerSize',40,'HandleVisibility','off')
plot(alldata.glider_lon(gliderind),alldata.glider_lat(gliderind),'.','color',rgb('Lime'),'MarkerSize',40,'HandleVisibility','off')

% plot(alldata.adcp_lon(1),alldata.adcp_lat(1),'color',rgb('Black'),'linewidth',3)
t = plot(alldata.float_lon(1),alldata.float_lat(1),'.','color',rgb('Red'),'MarkerSize',40);
b = plot(alldata.WF_lon(1),alldata.WF_lat(1),'.','color',rgb('Blue'),'MarkerSize',40);
g = plot(alldata.glider_lon(1),alldata.glider_lat(1),'.','color',rgb('Lime'),'MarkerSize',40);

legend({'','','EM-APEX','Wire Flyer','Glider'},'FontWeight','Bold','FontSize',40)
set(ax1,'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
axis([minlon,maxlon,minlat,maxlat]);
clim([-.5,1.5])
% title(string(current_date_ssh),'FontWeight','Bold','FontSize',30)
auxh.FontSize = 44;
auxh.FontWeight = 'bold';
ylabel(auxh,'SSH (m)','FontSize',80)
set(ax1,'XTick',[12.0 16.0 20.0 24.0],'XTickLabel',{['12',char(176),'E'],['16',char(176),'E'],['20',char(176),'E'],['24',char(176),'E']})
set(ax1,'YTick',[-40.0 -36.0 -32.0],'YTickLabel',{['40',char(176),'S'],['36',char(176),'S'],['32',char(176),'S']})
ax1.XAxis.FontSize = 80;
ax1.YAxis.FontSize = 80;
ax1.XAxis.FontWeight = 'bold';
ax1.YAxis.FontWeight = 'bold';

saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/data_map.png');


