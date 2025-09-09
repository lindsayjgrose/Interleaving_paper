clear
close all

addpath(genpath("/Users/lindsay.grose/Documents/URI/git.nosync/"))
addpath(genpath("/Users/lindsay.grose/Documents/MATLAB/"))


%%
D = dir;
cruise_path = '/Users/lindsay.grose/Documents/MATLAB/Datasets/WireFlyer/QUICCHE/Dives/Processed/';
cd(cruise_path)
pat = '(\d\d\d\d\d\d\d\d[_]\d\d\d\d\d\d)';
dir_names = []
for ii = 1:length(D) 
  a = regexp(D(ii).name,pat);
  if(a==1)
     dir_names = [dir_names;D(ii).name]  
  end
end

minlon = 12.5;
maxlon = 21;
minlat = -40;
maxlat = -34;

minlon = 12.5;
maxlon = 19;
minlat = -40;
maxlat = -35;


path = '/Users/lindsay.grose/Documents/MATLAB/Datasets/SSH.nosync/';
file = 'Cape_Basin_2023_03_04_2023_11_13.nc';
ssh_data = ncload([path,file]);

ogdate = datetime([2023,03,04,0,0,0],'timezone','UTC','format','yyyy-MM-dd HH:mm:ss');

trial =zeros(length(ssh_data.time));
ssh_data.datenum = trial(:,1)

for i = 1:length(ssh_data.time)
    ssh_data.datenum(i,1) = datenum(ogdate + days(i-1))
end

ssh_data.datetime = datetime(ssh_data.datenum,'ConvertFrom','datenum');


mymap = brewermap([],'Dark2');

colorwf = mymap(3,:);
colorwf = [0,0,0];

%%
figure
set(gcf,'Position',[4 1 2155 1336],'color','w')
mydives = [8,7,5];
for i=1:3
    cd([pwd,'/',dir_names(mydives(i),:)])	
    DD = dir([pwd,'/*_gridded_data.mat']);
    load([DD.folder,'/',DD.name]);
    DD = dir([pwd,'/*merged.mat']);
    load([DD.folder,'/',DD.name]);

    lon = gps_gprmc_t_GPS_GPRMC_DATA.lon;
    lat = gps_gprmc_t_GPS_GPRMC_DATA.lat;

    sshind = findClosest(median(unixtime2datetime(gps_processed_DATA.timestamp/1e6),'omitmissing'),ssh_data.datetime);
    fprintf(1,'%d \n',sshind);

    ax = subplot(4,3,i);
    ttt = pcolor(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,sshind)');
    set(ttt,'EdgeColor','none')
    colormap(ax,cmocean('rain'))
    axis(ax,[minlon,maxlon,minlat,maxlat])
    clim([-0.25,1.5])
    hold on 
    contour(ssh_data.longitude,ssh_data.latitude,ssh_data.adt(:,:,sshind)',[1.5 1.25 1 0.75 0.5 0.1],'k',"ShowText",true,'LineWidth',2);
    plot(lon,lat,'color',colorwf,'LineWidth',20)
    ind = find(~isnan(lon));
    if i==2
        plot(lon(ind(1)),lat(ind(1)),'.','color','Red',MarkerSize=80)
        plot(lon(ind(end)),lat(ind(end)),'.','color','Green',MarkerSize=80)
    else
        plot(lon(ind(1)),lat(ind(1)),'.','color','Green',MarkerSize=80)
        plot(lon(ind(end)),lat(ind(end)),'.','color','Red',MarkerSize=80)
    end

    minlon = min(lon(ind(1)),lon(ind(end)))-1.5;
    maxlon = max(lon(ind(1)),lon(ind(end)))+1.5;
    minlat = min(lat(ind(1)),lat(ind(end)))-1.5;
    maxlat = max(lat(ind(1)),lat(ind(end)))+1.5;

    axis(ax,[minlon,maxlon,minlat,maxlat])
    if i==1
        xLimits = ax.XLim;
        yLimits = ax.YLim;
        text(xLimits(1)-1.7, yLimits(2), 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xLimits = ax.XLim;
        yLimits = ax.YLim;
        text(xLimits(1)-2, yLimits(2), 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xLimits = ax.XLim;
        yLimits = ax.YLim;
        text(xLimits(1)-2.2, yLimits(2), 'c)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end

    % if i == 1
    %     text(lon(ind(1)),lat(ind(1))+.4,'S','color',rgb('Red'),'FontSize',28,'FontWeight','Bold');
    %     text(lon(ind(end))+.3,lat(ind(end)),'E','color',rgb('Red'),'FontSize',28,'FontWeight','Bold');
    % end
    % if i == 2
    %     text(lon(ind(1))-.5,lat(ind(1))+.5,'S','color',rgb('Red'),'FontSize',28,'FontWeight','Bold');
    %     text(lon(ind(end))-.5,lat(ind(end))+.4,'E','color',rgb('Red'),'FontSize',28,'FontWeight','Bold');
    % end
    % if i == 3
    %     text(lon(ind(1))+.2,lat(ind(1))+.1,'S','color',rgb('Red'),'FontSize',28,'FontWeight','Bold');
    %     text(lon(ind(end))+.15,lat(ind(end))+.4,'E','color',rgb('Red'),'FontSize',28,'FontWeight','Bold');
    % end
    set(ax,'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
    ax.FontWeight = 'Bold'; 
    pos = get(ax, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    pos(4) = 0.22;  % Set the height of all subplots to be equal
    set(ax, 'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  ); 
    set(ax,'XTick',[14.0 16.0 18.0 20.0],'XTickLabel',{['14',char(176),'E'],['16',char(176),'E'],['18',char(176),'E'],['20',char(176),'E']})
    set(ax,'YTick',[-40.0 -38.0 -36.0 -34.0 -32.0],'YTickLabel',{['40',char(176),'S'],['38',char(176),'S'],['36',char(176),'S'],['34',char(176),'S'],['32',char(176),'S']})
    % set(ax1,'XTick',[12.0 16.0 20.0],'XTickLabel',{['12',char(176),'E'],['16',char(176),'E'],['20',char(176),'E']})
    % set(ax1,'YTick',[-40.0 -36.0 -32.0],'YTickLabel',{['40',char(176),'S'],['36',char(176),'S'],['32',char(176),'S']})
    ax.XAxis.FontSize = 26;
    ax.YAxis.FontSize = 26;
    if i==3
        cb2 = colorbar;
        cb2.FontSize = 30;
        ylabel(cb2, 'SSH (m)','FontWeight','Bold','FontSize',44) 
    end
    % title(string(dateshift(median(unixtime2datetime(gps_processed_DATA.timestamp/1e6),'omitmissing'),'start','day')),'FontWeight','Bold','FontSize',36)
    if i ==3
        title('Cyclone-Anticyclone Edge','FontWeight','Bold','FontSize',36)
        set(ax, 'Position', pos +[-.056,-0.02,0,0]); 
    end
    if i ==2
        title('Seamount City','FontWeight','Bold','FontSize',36)
        set(ax, 'Position', pos +[-.065,-0.02,0,0]); 
    end    
    if i ==1
        title('Anticyclone','FontWeight','Bold','FontSize',36)
        set(ax, 'Position', pos +[-.072,-0.02,0,0]); 
    end

    % get rid of beginning and end profiles
    if i==3
        num = 200:1:285;
    else
        num = 285:1:380;
    end
    ind200 = findClosest(num,WF_data.conservative_temperature.depth);
    [blah,badcol] = find(~isnan(WF_data.conservative_temperature.data(ind200,:)));
    badcol = unique(badcol);

    % Find the columns NOT in badcol
    allCols = 1:size(WF_data.conservative_temperature.data, 2); % All column indices
    goodcols = setdiff(allCols, badcol);

    %TEMPERATURE PLOT
    map =cmocean('thermal');
    map([1:2:50,170:2:end],:) = [];
    ax1 = subplot(4,3,i+3);
    if i==2
        ttt = pcolor(max(WF_data.conservative_temperature.distance(goodcols))-WF_data.conservative_temperature.distance(goodcols),WF_data.conservative_temperature.depth,WF_data.conservative_temperature.data(:,goodcols));
    else
        ttt = pcolor(WF_data.conservative_temperature.distance(goodcols),WF_data.conservative_temperature.depth,WF_data.conservative_temperature.data(:,goodcols));
    end
    set(ttt,'EdgeColor','none')
    ylim([285,815]);
    set(ax1,'YDir','reverse');
    ax1.XAxis.FontSize = 26;
    ax1.YAxis.FontSize = 26;
    colormap(ax1,map); %flipud(slanCM('RdYlBu')))
    clim(ax1,[4.5,13.4])
    % xlabel('Distance (km)','FontWeight','Bold','FontSize',16)
    pos = get(ax1, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    pos(4) = 0.21;  % Set the height of all subplots to be equal
    xticklabels([]);
    if i==1
        ylabel(ax1, 'Depth (m)','FontWeight','Bold','FontSize',40)
        set(ax1, 'Position', pos +[-0.072,-0.04,0,0]);
    end
    if i==2
        yticklabels([]);
        set(ax1, 'Position', pos +[-.065,-0.04,0,0]);
    end
    if i==3
        yticklabels([]);
        auxh = colorbar(ax1);
        auxh.FontSize = 26;
        ylabel(auxh,{'Temperature',['(',char(176),'C)']},'FontWeight','Bold','FontSize',40)
        set(ax1, 'Position', pos +[-.056,-0.04,0,0]);
    end
    % xlim([0,115])
    if i==1
        xlim([0,105])
    end
    if i==1
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-430, 'd)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+10, yLimits(2)-430, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-430, 'f)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
 
 
    % %SALINITY
    ax2 = subplot(4,3,i+6);
    if i==2
        ttt = pcolor(max(WF_data.absolute_salinity.distance(goodcols))-WF_data.absolute_salinity.distance(goodcols),WF_data.absolute_salinity.depth,WF_data.absolute_salinity.data(:,goodcols));
    else
        ttt = pcolor(WF_data.absolute_salinity.distance(goodcols),WF_data.absolute_salinity.depth,WF_data.absolute_salinity.data(:,goodcols));
    end
    set(ttt,'EdgeColor','none')
    ylim([285,815]);
    set(ax2,'YDir','reverse');
    ax2.XAxis.FontSize = 26;
    ax2.YAxis.FontSize = 26;
    colormap(ax2,cmocean('haline')); %flipud(slanCM('RdYlBu')))
    clim(ax2,[34.45,35.35])
    xticklabels([]);
    % xlabel('Distance (km)','FontWeight','Bold','FontSize',16)
    pos = get(ax2, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    pos(4) = 0.21;  % Set the height of all subplots to be equal
    if i==1
        ylabel(ax2, 'Depth (m)','FontWeight','Bold','FontSize',40)
        set(ax2, 'Position', pos+[-0.074,-.04,0,0]);
    end
    if i==2
        yticklabels([]);
        set(ax2, 'Position', pos+[-0.065,-.04,0,0]);
    end
    if i==3
        yticklabels([]);
        auxh = colorbar(ax2);
        auxh.FontSize = 26;
        ylabel(auxh,{'Salinity', '(g/kg)'},'FontWeight','Bold','FontSize',40)
        set(ax2, 'Position', pos+[-0.056,-.04,0,0]);
    end
    % xlim([0,115])
    if i==1
        xlim([0,105])
    end
    if i==1
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-430, 'g)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+10, yLimits(2)-430, 'h)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-430, 'i)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end

    % oxygen
    map = slanCM('ocean');
    map(1:70,:) = [];

    ax3 = subplot(4,3,i+9);
    if i==2
        ttt = pcolor(max(WF_data.oxygen_concentration.distance(goodcols))-WF_data.oxygen_concentration.distance(goodcols),WF_data.oxygen_concentration.depth,WF_data.oxygen_concentration.data(:,goodcols));
    else
        ttt = pcolor(WF_data.oxygen_concentration.distance(goodcols),WF_data.oxygen_concentration.depth,WF_data.oxygen_concentration.data(:,goodcols));
    end
    set(ttt,'EdgeColor','none')
    ylim([285,815]);
    set(ax3,'YDir','reverse');
    ax3.XAxis.FontSize = 26;
    ax3.YAxis.FontSize = 26;
    colormap(ax3,map); %flipud(slanCM('RdYlBu')))
    clim(ax3,[190,230])
    xlabel('Distance (km)','FontWeight','Bold','FontSize',40)
    pos = get(ax3, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    pos(4) = 0.21;  % Set the height of all subplots to be equal
    if i==1
        ylabel(ax3, 'Depth (m)','FontWeight','Bold','FontSize',40)
    end
    if i==2
        yticklabels([]);
        set(ax3, 'Position', pos+[-0.065,-0.04,0,0]);
    end
    if i==3
        yticklabels([]);
        auxh = colorbar(ax3);
        auxh.FontSize = 26;
        ylabel(auxh,{'Oxygen', '(umol/kg)'},'FontWeight','Bold','FontSize',40)
        set(ax3, 'Position', pos+[-0.056,-0.04,0,0]);
    end
    % xlim([0,115])
    if i==1
        xlim([0,105])
        hold on
        plot(105,815,'.','color','Red',MarkerSize=60)
        set(ax3, 'Position', pos+[-0.074,-0.04,0,0]);
    end
    if i==1
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-430, 'j)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+10, yLimits(2)-430, 'k)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-430, 'l)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    hold on
    plot(0,815,'.','color','Green',MarkerSize=60)
    if i~=1
        plot(max(WF_data.oxygen_concentration.distance(goodcols)),815,'.','color','Red',MarkerSize=60)
    end


    cd(cruise_path)
end

saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/dive_panel_depth.svg');


%% density 
% set(gcf,'Position',[1 1 2509 1295],'color','w')

figure
set(gcf,'Position',[4 1 2155 1336],'color','w')
mydives = [8,7,5];
for i=1:3
    cd([pwd,'/',dir_names(mydives(i),:)])	
    DD = dir([pwd,'/*_derived_quantities.mat']);
    load([DD.folder,'/',DD.name]);

  
    ax = subplot(3,3,i);

    if i==2
        ttt = pcolor(max(max(data.density.dist))-data.density.dist,data.density.potential_density,data.density.SA);
    else
        ttt = pcolor(data.density.dist,data.density.potential_density,data.density.SA);
    end
    % ttt = pcolor(data.density.dist,data.density.potential_density,data.density.SA);
    set(ttt,'EdgeColor','none')
    colormap(ax,cmocean('haline')); %flipud(slanCM('RdYlBu')))
    ylim(ax,[26.35,27.35])
    clim(ax,[34.45,35.35])
    ax.XAxis.FontSize = 26;
    ax.YAxis.FontSize = 26;
    xticklabels([]);
    pos = get(ax, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    % pos(4) = 0.21;  % Set the height of all subplots to be equal
    pos(4) = 0.29;  % Set the height of all subplots to be equal
    if i==1
        ylabel(ax, '\sigma (kg/m^3)','FontWeight','Bold','FontSize',40)
        % set(ax, 'Position', pos +[-0.072,-0.04,0,0],'YDir','Reverse');
        set(ax, 'Position', pos +[-0.072,-0.03,0,0],'YDir','Reverse');
    end
    if i==2
        yticklabels([]);
        % set(ax, 'Position', pos +[-0.0685,-0.04,0,0],'YDir','Reverse');
        set(ax, 'Position', pos +[-0.0685,-0.03,0,0],'YDir','Reverse');
    end
    if i==3
        yticklabels([]);
        auxh = colorbar(ax);
        auxh.FontSize = 26;
        ylabel(auxh,{'Salinity', '(g/kg)'},'FontWeight','Bold','FontSize',40)
        % set(ax, 'Position', pos +[-0.065,-0.04,0,0],'YDir','Reverse');
        set(ax, 'Position', pos +[-0.065,-0.03,0,0],'YDir','Reverse');
    end
    if i==1
        xlim([0,105])
    end
    if i==1
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-.85, 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xlim([0,115])
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+10, yLimits(2)-.85, 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xlim([0,110])
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-.85, 'c)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i == 3
        t = title('Cyclone-Anticyclone Edge','FontWeight','Bold','FontSize',36);
        t.Position(2) = t.Position(2) + 0.01; % Move the title down (adjust the value as needed)
    end
    if i ==2
        t = title('Seamount City','FontWeight','Bold','FontSize',36);
        t.Position(2) = t.Position(2) + 0.01; % Move the title down (adjust the value as needed)
    end    
    if i ==1
        t = title('Anticyclone','FontWeight','Bold','FontSize',36);
        t.Position(2) = t.Position(2) + 0.01; % Move the title down (adjust the value as needed)
    end
  


    %DSC PLOT
    ax1 = subplot(3,3,i+3);
    
    if i==2
        ttt = pcolor(max(max(data.density.dist))-data.density.dist,data.density.potential_density,data.density.DSC);
    else
        ttt = pcolor(data.density.dist,data.density.potential_density,data.density.DSC);
    end
    % ttt = pcolor(data.density.dist,data.density.potential_density,data.density.DSC);
    set(ttt,'EdgeColor','none')
    set(ax1,'YDir','reverse');
    ylim(ax1,[26.35,27.35])
    ax1.XAxis.FontSize = 26;
    ax1.YAxis.FontSize = 26;
    colormap(ax1,cmocean('balance')); %flipud(slanCM('RdYlBu')))
    clim(ax1,[-150,150])
    xticklabels([]);
    % xlabel('Distance (km)','FontWeight','Bold','FontSize',16)
    pos = get(ax1, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    % pos(4) = 0.21;  % Set the height of all subplots to be equal
    pos(4) = 0.29;  % Set the height of all subplots to be equal
    if i==1
        ylabel(ax1, '\sigma (kg/m^3)','FontWeight','Bold','FontSize',40)
        set(ax1, 'Position', pos +[-0.072,-0.04,0,0]);
    end
    if i==2
        yticklabels([]);
        set(ax1, 'Position', pos +[-0.0685,-0.04,0,0]);
    end
    if i==3
        yticklabels([]);
        auxh = colorbar(ax1);
        auxh.FontSize = 26;
        ylabel(auxh,{'Spiciness','Curvature (m^3/kg)'},'FontWeight','Bold','FontSize',40)
        set(ax1, 'Position', pos +[-0.065,-0.04,0,0]);
    end
    if i==1
        xlim([0,105])
    end
    if i==1
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-.85, 'd)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xlim([0,115])
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+10, yLimits(2)-.85, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xlim([0,110])
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-.85, 'f)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end

    % Turner Angle
    vmin = -100;       % Minimum data value
    vmax = 100;        % Maximum data value
    flat_min = -45;   % Start of the flat region
    flat_max = 45;    % End of the flat region
    
    % Define the custom colormap
    n = 256; % Number of colors
    cmap = brewermap([],'RdYlBu'); % Use the 'jet' colormap as a base
    
    % Map the flat region to a single color (e.g., the middle of 'jet')
    flat_color = cmap(round(n/2), :);
    
    % Adjust the colormap
    x = linspace(vmin, vmax, n); % Full range of data
    cmap_adjusted = interp1(x, cmap, x, 'linear', 'extrap'); % Interpolate
    cmap_adjusted(x >= flat_min & x <= flat_max, :) = repmat(flat_color, sum(x >= flat_min & x <= flat_max), 1);

    ax2 = subplot(3,3,i+6);
         % Define shading limits (example: shade between 5 and 10 days)
     x1 = 0;
    if i==1
        x2 = 105;
    end
    if i>1
        x2 = 115;
    end
    yLimits = [26, 26.5]; % or use ylim if you want dynamic adjustment
    hold on
    % fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
    %      rgb('LightGray'), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    % yline(26.5,'Color',rgb('LightGray')-.2,'LineWidth',2)
    % Define shading limits (example: shade between 5 and 10 days)
    yLimits = [26.5, 27]; % or use ylim if you want dynamic adjustment
    hold on
    % fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
    %      rgb('LightGray')-.2, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    if i==2
        yline(27,'Color',[0 0 0],'LineWidth',2)
    end
    yLimits = [27, 27.5]; % or use ylim if you want dynamic adjustment
    % yline(27,'Color',rgb('LightGray')-.4,'LineWidth',2)
    hold on
    % fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
    %      rgb('LightGray')-.4, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on
    if i==2
        ttt = pcolor(max(max(data.density.dist))-data.density.dist,data.density.potential_density,data.density.Tu);
    else
        ttt = pcolor(data.density.dist,data.density.potential_density,data.density.Tu);
    end
    % ttt = pcolor(data.density.dist,data.density.potential_density,data.density.Tu);
    set(ttt,'EdgeColor','none')
    set(ax2,'YDir','reverse');
    ylim(ax2,[26.35,27.35])
    ax2.XAxis.FontSize = 26;
    ax2.YAxis.FontSize = 26;
    colormap(ax2,cmap_adjusted); %flipud(slanCM('RdYlBu')))
    clim(ax2,[vmin,vmax])
    % xticklabels([]);
    xlabel('Distance (km)','FontWeight','Bold','FontSize',40)
    pos = get(ax2, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    % pos(4) = 0.21;  % Set the height of all subplots to be equal
    pos(4) = 0.29;  % Set the height of all subplots to be equal
    if i==1
        ylabel(ax2, '\sigma (kg/m^3)','FontWeight','Bold','FontSize',40)
        % set(ax2, 'Position', pos+[-0.072,-0.08,0,0]);
        set(ax2, 'Position', pos+[-0.072,-0.05,0,0]);
    end
    if i==2
        yticklabels([]);
        set(ax2, 'Position', pos+[-0.0685,-0.05,0,0]);
    end
    if i==3
        yticklabels([]);
        auxh = colorbar(ax2);
        auxh.FontSize = 26;
        ylabel(auxh,{'Turner Angle', ['(',char(176),')']},'FontWeight','Bold','FontSize',40)
        auxh.Ticks = [-100,-45,45,100]; % Positions of the ticks
        set(ax2, 'Position', pos+[-0.065,-0.05,0,0]);
    end
    if i==1
        xlim([0,105])
    end
    if i==1
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-.85, 'g)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xlim([0,115])
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+10, yLimits(2)-.85, 'h)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xlim([0,110])
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-.85, 'i)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    box on
    % Strain
    % ax3 = subplot(4,3,i+9);
    % ttt = pcolor(data.density.dist,data.density.potential_density,data.density.strain);
    % set(ttt,'EdgeColor','none')
    % set(ax3,'YDir','reverse');
    % ylim(ax3,[26.35,27.35])
    % ax3.XAxis.FontSize = 26;
    % ax3.YAxis.FontSize = 26;
    % colormap(ax3,cmocean('curl')); 
    % clim(ax3,[-1,1])
    % xlabel('Distance (km)','FontWeight','Bold','FontSize',40)
    % if i==1
    %     ylabel(ax3, '\sigma (kg/m^3)','FontWeight','Bold','FontSize',40)
    % end
    % if i==2
    %     yticklabels([]);
    % end
    % if i==3
    %     yticklabels([]);
    %     auxh = colorbar(ax3);
    %     auxh.FontSize = 26;
    %     ylabel(auxh,{'Strain'},'FontWeight','Bold','FontSize',40)
    % end
    % pos = get(ax3, 'Position');
    % pos(3) = 0.25;  % Set the width of all subplots to be equal
    % pos(4) = 0.21;  % Set the height of all subplots to be equal
    % set(ax3, 'Position', pos+[-0.08+(i/80),-0.03,0,0]);
    % if i==3
    %     xlim([0,105])
    % end
    % if i==1
    %     xLimits = ax3.XLim;
    %     yLimits = ax3.YLim;
    %     text(xLimits(1)-18, yLimits(2)-.95, 'd)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    % end

    % cmap = cividis();
    % % Buoyancy frequency 
    % ax3 = subplot(4,3,i+9);
    % ttt = pcolor(data.density.dist,data.density.potential_density,real(sqrt(data.density.N2)));
    % set(ttt,'EdgeColor','none')
    % colormap(ax3,cmap); %flipud(slanCM('RdYlBu')))
    % ylim(ax3,[26.35,27.35])
    % clim(ax3,[2.9e-3,4.7e-3])
    % ax3.XAxis.FontSize = 26;
    % ax3.YAxis.FontSize = 26;
    % xticklabels([]);
    % if i==1
    %     ylabel(ax3, 'Depth (m)','FontWeight','Bold','FontSize',40)
    % end
    % if i==2
    %     yticklabels([]);
    % end
    % if i==3
    %     yticklabels([]);
    %     auxh = colorbar(ax3);
    %     auxh.FontSize = 26;
    %     ylabel(auxh,{'Buoyancy', 'Frequency (s^{-1})'},'FontWeight','Bold','FontSize',40)
    % end
    % pos = get(ax3, 'Position');
    % pos(3) = 0.25;  % Set the width of all subplots to be equal
    % pos(4) = 0.21;  % Set the height of all subplots to be equal
    % set(ax3, 'Position', pos+[-0.08+(i/80),-0.01,0,0],'YDir','Reverse'); 
    % if i==3
    %     xlim([0,105])
    % end
    % if i==1
    %     xLimits = ax3.XLim;
    %     yLimits = ax3.YLim;
    %     text(xLimits(1)-18, yLimits(2)-500, 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    % end


    cd(cruise_path)
end

saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/dive_panel_denisty.svg');
% saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/dive_panel_denisty_shade.svg');

%% shear and others in depth space


figure
set(gcf,'Position',[4 1 2155 1336],'color','w')
mydives = [8,7,5];
for i=1:3
    cd([pwd,'/',dir_names(mydives(i),:)])	
    DD = dir([pwd,'/*_derived_quantities.mat']);
    load([DD.folder,'/',DD.name]);

    cmap = cividis();
    % Buoyancy frequency 
    ax = subplot(3,3,i+3);
    if i==2
        ttt = pcolor(max(max(data.depth.dist))-data.depth.dist,data.depth.depth,real(sqrt(data.depth.N2)));
    else
        ttt = pcolor(data.depth.dist,data.depth.depth,real(sqrt(data.depth.N2)));
    end
    % ttt = pcolor(data.depth.dist,data.depth.depth,real(sqrt(data.depth.N2)));
    set(ttt,'EdgeColor','none')
    colormap(ax,cmap); %flipud(slanCM('RdYlBu')))
    ylim(ax,[285,815])
    clim(ax,[2.5e-3,4.7e-3])
    ax.XAxis.FontSize = 26;
    ax.YAxis.FontSize = 26;
    xticklabels([]);
    pos = get(ax, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    % pos(4) = 0.21;  % Set the height of all subplots to be equa;
    pos(4) = 0.28;  % Set the height of all subplots to be equa;
    if i==1
        ylabel(ax, 'Depth (m)','FontWeight','Bold','FontSize',40)
        set(ax, 'Position', pos+[-0.072,-0.035,0,0],'YDir','Reverse'); 
    end
    if i==2
        yticklabels([]);
        set(ax, 'Position', pos+[-0.0685,-0.035,0,0],'YDir','Reverse'); 
    end
    if i==3
        yticklabels([]);
        auxh = colorbar(ax);
        auxh.FontSize = 26;
        ylabel(auxh,{'Buoyancy', 'Frequency (s^{-1})'},'FontWeight','Bold','FontSize',40)
        set(ax, 'Position', pos+[-0.065,-0.035,0,0],'YDir','Reverse'); 
    end
    if i==1
        xlim([0,105])
    end
    if i==1
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-450, 'd)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+10, yLimits(2)-450, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-450, 'f)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end


    %DSC PLOT
    ax1 = subplot(3,3,i);
    if i==2
        ttt = pcolor(max(max(data.depth.dist))-data.depth.dist,data.depth.depth,data.depth.DSC);
    else
        ttt = pcolor(data.depth.dist,data.depth.depth,data.depth.DSC);
    end
    % ttt = pcolor(data.depth.dist,data.depth.depth,data.depth.DSC);
    set(ttt,'EdgeColor','none')
    set(ax1,'YDir','reverse');
    ylim(ax1,[285,815])
    ax1.XAxis.FontSize = 26;
    ax1.YAxis.FontSize = 26;
    colormap(ax1,cmocean('balance')); %flipud(slanCM('RdYlBu')))
    pos = get(ax1, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    % pos(4) = 0.21;  % Set the height of all subplots to be equal   
    pos(4) = 0.28;  % Set the height of all subplots to be equal   
    clim(ax1,[-150,150])
    hold on
    if i==1
        [C, h] = contour(data.depth.dist(80:end-50,15:end-15),data.depth.depth(80:end-50,15:end-15),data.depth.potential_density(80:end-50,15:end-15),[26.4 26.5 26.6 26.7 26.8 26.9 27],'color', rgb('Cyan'),"ShowText",false,'LineWidth',4);
    end
    if i==2
        [C, h] = contour(max(max(data.depth.dist))-data.depth.dist(100:end-65,15:end-15),data.depth.depth(100:end-65,15:end-15),data.depth.potential_density(100:end-65,15:end-15),[26.7 26.8 26.9 27 27.1],'color',rgb('Cyan'),"ShowText",false,'LineWidth',4);
    end
    if i==3
        [C, h] = contour(data.depth.dist(51:end-50,15:end-15),data.depth.depth(51:end-50,15:end-15),data.depth.potential_density(51:end-50,15:end-15),[26.6 26.7 26.8 26.9 27],'color',rgb('Cyan'),"ShowText",false,'LineWidth',4);
    end
        % clabel(C, h, 'LabelSpacing', 300); % Set distance between labels
    xticklabels([]);
    if i==1
        ylabel(ax1, 'Depth (m)','FontWeight','Bold','FontSize',40)
        set(ax1, 'Position', pos +[-0.072,-0.02,0,0]);
    end
    if i==2
        yticklabels([]);
        set(ax1, 'Position', pos +[-0.0685,-0.02,0,0]);
    end
    if i==3
        yticklabels([]);
        auxh = colorbar(ax1);
        auxh.FontSize = 26;
        ylabel(auxh,{'Spiciness','Curvature (m^3/kg)'},'FontWeight','Bold','FontSize',40)
        set(ax1, 'Position', pos +[-0.065,-0.02,0,0]);
    end
    if i==1
       xlim([0,105])
    end
    if i==1
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-450, 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+10, yLimits(2)-450, 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-450, 'c)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i ==3
        t = title('Cyclone-Anticyclone Edge','FontWeight','Bold','FontSize',36);
        t.Position(2) = t.Position(2)+10; % Move the title down (adjust the value as needed)
    end
    if i ==2
        t = title('Seamount City','FontWeight','Bold','FontSize',36);
        t.Position(2) = t.Position(2)+10; % Move the title down (adjust the value as needed)
    end    
    if i ==1
        t = title('Anticyclone','FontWeight','Bold','FontSize',36);
        t.Position(2) = t.Position(2)+10; % Move the title down (adjust the value as needed)
    end


    % Along track vertical shear
    ax2 = subplot(3,3,i+6);
    if i==2
        ttt = pcolor(max(max(data.depth.ADCP.distance./1000))-data.depth.ADCP.distance./1000,data.depth.ADCP.depth,data.depth.ADCP.along_track_shear);
    else
        ttt = pcolor(data.depth.ADCP.distance./1000,data.depth.ADCP.depth,data.depth.ADCP.along_track_shear);
    end
    % ttt = pcolor(data.depth.ADCP.distance./1000,data.depth.ADCP.depth,data.depth.ADCP.along_track_shear);
    set(ttt,'EdgeColor','none')
    set(ax2,'YDir','reverse');
    ylim(ax2,[0,400])
    ax2.XAxis.FontSize = 26;
    ax2.YAxis.FontSize = 26;
    colormap(ax2,brewermap([],'PuOr')); %flipud(slanCM('RdYlBu')))
    clim(ax2,[-4e-3,4e-3])
    xlabel('Distance (km)','FontWeight','Bold','FontSize',40)
    % xticklabels([]);
    pos = get(ax2, 'Position');
    pos(3) = 0.28;  % Set the width of all subplots to be equal
    % pos(4) = 0.21;  % Set the height of all subplots to be equal
    pos(4) = 0.28;  % Set the height of all subplots to be equal
    if i==1
        ylabel(ax2, 'Depth (m)','FontWeight','Bold','FontSize',40)
        set(ax2, 'Position', pos+[-0.072,-.05,0,0]);
    end
    if i==2
        yticklabels([]);
        set(ax2, 'Position', pos+[-0.0685,-.05,0,0]);
    end
    if i==3
        yticklabels([]);
        auxh = colorbar(ax2);
        auxh.FontSize = 26;
        ylabel(auxh,{'Along Track', 'Shear (s^{-1})'},'FontWeight','Bold','FontSize',40)
        set(ax2, 'Position', pos+[-0.065,-.05,0,0]);
    end
    if i==1
        xlim([0,105])
    end
    if i==1
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-430, 'g)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==2
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+10, yLimits(2)-430, 'h)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end
    if i==3
        xLimits = ax1.XLim;
        yLimits = ax1.YLim;
        text(xLimits(1)+8, yLimits(2)-430, 'i)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',44)
    end

    % Across track vertical shear
    % ax3 = subplot(4,3,i+9);
    % ttt = pcolor(data.depth.ADCP.distance./1000,data.depth.ADCP.depth,data.depth.ADCP.across_track_shear);
    % set(ttt,'EdgeColor','none')
    % set(ax3,'YDir','reverse');
    % ylim(ax3,[0,400])
    % ax3.XAxis.FontSize = 26;
    % ax3.YAxis.FontSize = 26;
    % colormap(ax3,brewermap([],'PuOr')); 
    % clim(ax3,[-4e-3,4e-3])
    % xlabel('Distance (km)','FontWeight','Bold','FontSize',40)
    % if i==1
    %     ylabel(ax3, 'Depth (m)','FontWeight','Bold','FontSize',40)
    % end
    % if i==2
    %     yticklabels([]);
    % end
    % if i==3
    %     yticklabels([]);
    %     auxh = colorbar(ax3);
    %     auxh.FontSize = 26;
    %     ylabel(auxh,{'Across Track', 'Shear (s^{-1})'},'FontWeight','Bold','FontSize',40)
    % end
    % pos = get(ax3, 'Position');
    % pos(3) = 0.25;  % Set the width of all subplots to be equal
    % pos(4) = 0.21;  % Set the height of all subplots to be equal
    % set(ax3, 'Position', pos+[-0.08+(i/80),-0.03,0,0]);
    % if i==3
    %     xlim([0,105])
    % end


    cd(cruise_path)
end

saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/dive_panel_depth_shear.svg');
