clear 
close all
%% Calculate spice for Argo
addpath(genpath("/Users/lindsay.grose/Documents/URI/git.nosync/"))
addpath(genpath("/Users/lindsay.grose/Documents/MATLAB/"))

load('pot_density_grid.mat')
load('all_depth_grid.mat')
load('all_density_grid.mat')


ssh_30yr = load('CapeBasin_ssh_30yr.mat');
minlon = 0;
maxlon = 35;
minlat = -45;
maxlat = -20;

latind = find(ssh_30yr.Agulhas.lat >=minlat & ssh_30yr.Agulhas.lat <=maxlat);
lonind = find(ssh_30yr.Agulhas.lon >=minlon & ssh_30yr.Agulhas.lon <=maxlon);

ssh_data.latitude = ssh_30yr.Agulhas.lat(latind);
ssh_data.longitude = ssh_30yr.Agulhas.lon(lonind);
ssh_data.adt = ssh_30yr.Agulhas.adt(lonind,latind,:);

meanssh = mean(ssh_data.adt,3)';

%% smooth
windowSize = 7; % Window size (odd number)
blackmanWindow = blackman(windowSize);
for i=1:size(argo_density_grid.CT,2)
    argo_density_grid.CT_filt(:,i) = conv(argo_density_grid.CT(:,i), blackmanWindow, 'same') / sum(blackmanWindow);
end
for i=1:size(argo_density_grid.SA,2)
    argo_density_grid.SA_filt(:,i) = conv(argo_density_grid.SA(:,i), blackmanWindow, 'same') / sum(blackmanWindow);
end


%% calculate DSC

tau = gsw_spiciness0(argo_density_grid.SA_filt,argo_density_grid.CT_filt);
dtau = diff(tau,1,1);
dpot_rho = diff(argo_density_grid.pot_rho,1);
dtau_drho = dtau./dpot_rho';
dtau2 = diff(dtau_drho,1,1);
spice_curve = dtau2./dpot_rho(:,1:end-1)';

spice_curve = [NaN(1, size(spice_curve, 2)); spice_curve(1:end, :)];
spice_curve(size(spice_curve,1)+1,:) = NaN;

%instead of depth, do density
minrho = 1026;
maxrho = 1027.5;
indmin = findClosest(minrho,argo_density_grid.pot_rho);
indmax = findClosest(maxrho,argo_density_grid.pot_rho);
DSC_mag_avg = mean(abs(spice_curve(indmin:indmax,:)),1,'omitnan');


%%
minrhos = [1026.0 1026.5 1027.0];
maxrhos = [1026.5 1027.0 1027.5];

for t =1:length(minrhos)
    minrho = minrhos(t);
    maxrho = maxrhos(t);
    indmin = findClosest(minrho,argo_density_grid.pot_rho);
    indmax = findClosest(maxrho,argo_density_grid.pot_rho);
    DSC_mag_avg_int = mean(abs(spice_curve(indmin:indmax,:)),1,'omitnan');
    
    
    
    lats = -45:.5:-20;
    lons = 0:.5:35;
    
    [lats, lons] = meshgrid(lats,lons);
    lats = lats';
    lons = lons';
    
    store_profiles = {};
    
    ctr = 1;
    for i=1:size(lats,1)
        for j=1:size(lons,2)
            store_profiles.lat(ctr) = lats(i,1);
            store_profiles.lon(ctr) = lons(1,j);
            ind = find((argo_density_grid.lat >= (lats(i,1)-.25)) & (argo_density_grid.lat <= (lats(i,1)+.25)) & (argo_density_grid.lon >= (lons(1,j)-.25)) & (argo_density_grid.lon <= (lons(1,j)+.25)));
            temp_matix = DSC_mag_avg_int(ind);
            store_profiles.DSC_mag_Avg(ctr) = mean(temp_matix,'omitnan');
            ctr = ctr +1;
            clear temp_matix
        end
    end
    
    
    stor_mean = reshape(store_profiles.DSC_mag_Avg,[71,51])';
    
    figure
    set(gcf,'Position',[1 63 2392 1274],'color','w')
    ax1 = axes;
    ttt = pcolor(ssh_data.longitude,ssh_data.latitude,meanssh);
    set(ttt,'EdgeColor','none')
    colormap(ax1,cmocean('rain'))
    clim([-.25,1.5])
    axis(ax1,[minlon,maxlon,minlat,maxlat])
    ax2 = axes;
    hold on
    tttt = pcolor(lons,lats,stor_mean);
    set(tttt,'EdgeColor','none')
    colormap(ax2,brewermap([],"BuPu"))
    clim(ax2,[0,35])
    axis(ax2,[minlon,maxlon,minlat,maxlat])
    hold on
    contour(ssh_data.longitude,ssh_data.latitude,meanssh,[0.1 0.3 0.4 0.5 0.6 2],'k',"ShowText",true,'LineWidth',4);
    latind = find(ssh_data.latitude<-39.2);
    contour(ssh_data.longitude,ssh_data.latitude(latind),meanssh(latind,:),[0.3 10],'Color','Red','LineWidth',6);
    latind = find(ssh_data.latitude>-37.5);
    lonind = find(ssh_data.longitude<14.8);
    contour(ssh_data.longitude(lonind),ssh_data.latitude(latind),meanssh(latind,lonind),[0.5 10],'Color','Green','LineWidth',6);
    contour(ssh_data.longitude,ssh_data.latitude,meanssh,[.75 10],'Color','White','LineWidth',6);

    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    set([ax1,ax2],'Position',[.08 .08 .8 .8]);
    ax1.FontWeight = 'Bold'; 
    title(ax1,['\sigma = ',num2str(minrho-1000),' - ',num2str(maxrho-1000)],'FontWeight','Bold','FontSize',80)
    set([ax1,ax2],'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
    % set(ax1,'XTick',[0 5 10.0 15 20.0 25.0 30.0],'XTickLabel',{['0',char(176),'E'],['5',char(176),'E'],['10',char(176),'E'],['15',char(176),'E'],['20',char(176),'E'],['25',char(176),'E'],['30',char(176),'E']})
    % set(ax1,'YTick',[-45 -40.0 -35.0  -30.0 -25 -20],'YTickLabel',{['45',char(176),'S'],['40',char(176),'S'],['35',char(176),'S'],['30',char(176),'S'],['25',char(176),'S'],['20',char(176),'S']})
    % set(ax1,'XTick',[0 10.0 20.0 30.0],'XTickLabel',{['0',char(176),'E'],['10',char(176),'E'],['20',char(176),'E'],['30',char(176),'E']})
    set(ax1,'YTick',[-46.0 -35.0 -25],'YTickLabel',{['46',char(176),'S'],['35',char(176),'S'],['25',char(176),'S']})
    ax1.XAxis.FontSize = 100;
    ax1.YAxis.FontSize = 100;
    if t == 1
        % text(32.5,-21,'a)','FontWeight','Bold','FontSize',60)
        title(ax1,'Mode Waters','FontWeight','Bold','FontSize',100)
        xLimits = ax2.XLim;
        yLimits = ax2.YLim;
        text(xLimits(1)-3.5, yLimits(2)+1.2, 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',100)
    end
    if t == 2
        % text(32.5,-21,'b)','FontWeight','Bold','FontSize',60)
        title(ax1,'Thermocline Waters','FontWeight','Bold','FontSize',100)
        yticklabels(ax2,[]);
        yticklabels(ax1,[]);
    end
    if t ==3 
        yticklabels(ax2,[]);
        yticklabels(ax1,[]);
        title(ax1,'Intermediate Waters','FontWeight','Bold','FontSize',100)
        % text(32.5,-21,'c)','FontWeight','Bold','FontSize',60)
        % cb1 = colorbar(ax1,'northoutside','Position',[0.081369524720487,0.874375050709371,0.760846464059317,0.03]);
        cb2 = colorbar(ax2,'Position',[0.739417798029584,0.078756975360443,0.03,0.801582487710821]);
        % cb1.FontSize = 30;
        cb2.FontSize = 80;
        cb2.FontWeight = 'Bold';
        label = {sprintf('|DSC| (m^3 kg^{-1})')};
        ylabel(cb2,label,'FontWeight','Bold','FontSize',100); % Use the formatted string as the label
        % ylabel(cb2,['Average |DSC| ',string(minrho),' - ', string(maxrho),' (m^3 kg^{-1})'],'FontWeight','Bold','FontSize',44)
        % ylabel(cb1,'SSH (m)','FontWeight','Bold','FontSize',44)
    end
    xticklabels([ax1,ax2],[])
    label = sprintf('/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_DSC_%d.svg', t);
    % label = sprintf('/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/EGU_poster/Argo_DSC_%d.svg', t);
    % saveas(gcf, label);

end

%% calculate N2
pres_grid = meshgrid(argo_depth_grid.P',argo_depth_grid.CT(1,:))';
[N2,pmid] = gsw_Nsquared(argo_depth_grid.SA,argo_depth_grid.CT,pres_grid);


%% plot fN^2/g and its gradient

lats = -45:.5:-20;
lons = 0:.5:35;

[lats, lons] = meshgrid(lats,lons);
lats = lats';
lons = lons';

minrhos = [1026.0 1026.5 1027.0];
maxrhos = [1026.5 1027.0 1027.5];

% Convert degrees to radians
Lat_rad = deg2rad(lats);
Lon_rad = deg2rad(lons);

% Earth radius (in meters)
R = 6371000;

% Compute distance between adjacent latitudes (north-south spacing)
dlat = diff(Lat_rad, 1, 1);  % latitude spacing (vertical)
dlon = diff(Lon_rad, 1, 2);  % longitude spacing (horizontal)

% Compute meridional (north-south) distance between grid points
dy = R * dlat;  % meters

% Compute zonal (east-west) distance using mean latitude between grid rows
lat_center = 0.5 * (lats(1:end,1:end-1) + lats(1:end,2:end));
dx = R * cosd(lat_center) .* deg2rad(1);  % meters

% Preallocate cumulative grids
[nrows, ncols] = size(dx);
x = zeros(nrows, ncols);  % cumulative east distance
[nrows, ncols] = size(dy);
y = zeros(nrows, ncols);  % cumulative north distance

% Build cumulative x distance from left to right (along rows)
for i = 1:nrows
    x(i,2:end) = cumsum(dx(i,1:end-1));
end

% Build cumulative y distance from bottom to top (along columns)
for j = 1:ncols
    y(2:end,j) = cumsum(dy(1:end-1,j));
end

PVx = nan(size(lats,1),size(lats,2),length(minrhos));
PVy = nan(size(lats,1),size(lats,2),length(minrhos));

for t =1:length(minrhos)
    minrho = minrhos(t);
    maxrho = maxrhos(t);
    N2_avg_int = nan(1,size(argo_depth_grid.pot_rho,2));
    for i = 1:length(N2_avg_int)
        if sum(isfinite(argo_depth_grid.pot_rho(:,i))) > 2
            indmin = findClosest(minrho,argo_depth_grid.pot_rho(:,i));
            indmax = findClosest(maxrho,argo_depth_grid.pot_rho(:,i));
            N2_avg_int(i) = mean(N2(indmin:indmax,i),'omitnan');
        end
    end
  
    PV_avg = N2_avg_int .* (gsw_f(argo_depth_grid.lat)/9.81);

    lats = -45:.5:-20;
    lons = 0:.5:35;
    
    [lats, lons] = meshgrid(lats,lons);
    lats = lats';
    lons = lons';
    
    store_profiles = {};
    
    ctr = 1;
    for i=1:size(lats,1)
        for j=1:size(lons,2)
            store_profiles.lat(ctr) = lats(i,1);
            store_profiles.lon(ctr) = lons(1,j);
            ind = find((argo_density_grid.lat >= (lats(i,1)-.25)) & (argo_density_grid.lat <= (lats(i,1)+.25)) & (argo_density_grid.lon >= (lons(1,j)-.25)) & (argo_density_grid.lon <= (lons(1,j)+.25)));
            temp_matix = PV_avg(ind);
            store_profiles.PV_Avg(ctr) = mean(temp_matix,'omitnan');
            ctr = ctr +1;
            clear temp_matix
        end
    end
    
    
    stor_mean = reshape(store_profiles.PV_Avg,[71,51])';
    
    figure
    set(gcf,'Position',[1 63 2392 1274],'color','w')
    ax1 = axes;
    ttt = pcolor(ssh_data.longitude,ssh_data.latitude,meanssh);
    set(ttt,'EdgeColor','none')
    colormap(ax1,cmocean('rain'))
    clim([-.25,1.5])
    axis(ax1,[minlon,maxlon,minlat,maxlat])
    ax2 = axes;
    hold on
    tttt = pcolor(lons,lats,(stor_mean.*-1));
    set(tttt,'EdgeColor','none')
    colormap(ax2,brewermap([],"YlOrRd"))
    if t == 1
        % clim(ax2,[2.3e-10,5e-10])
        clim(ax2,[2.5e-10,4.5e-10])
        % clim(ax2,[1e-10,4e-10])
        % clim(ax2,[1e-10,2e-10])
        % title(ax2,'Mode Waters','FontWeight','Bold','FontSize',80)
        xLimits = ax2.XLim;
        yLimits = ax2.YLim;
        text(xLimits(1)-3.5, yLimits(2), 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',100)

    end
    if t == 2
        clim(ax2,[1e-10,3e-10])
        % clim(ax2,[1e-10,2e-10])
        % title(ax2,'Thermocline Waters','FontWeight','Bold','FontSize',80)
    end
    if t == 3
        clim(ax2,[.5e-10,2.5e-10])
        % title(ax2,'Intermediate Waters','FontWeight','Bold','FontSize',80)
    end
    axis(ax2,[minlon,maxlon,minlat,maxlat])
    hold on
    contour(ssh_data.longitude,ssh_data.latitude,meanssh,[0.1 0.3 0.4 0.5 0.6 2],'k',"ShowText",true,'LineWidth',4);
    latind = find(ssh_data.latitude<-39.2);
    contour(ssh_data.longitude,ssh_data.latitude(latind),meanssh(latind,:),[0.3 10],'Color','Red','LineWidth',6);
    latind = find(ssh_data.latitude>-37.5);
    lonind = find(ssh_data.longitude<14.8);
    contour(ssh_data.longitude(lonind),ssh_data.latitude(latind),meanssh(latind,lonind),[0.5 10],'Color','Green','LineWidth',6);
    contour(ssh_data.longitude,ssh_data.latitude,meanssh,[.75 10],'Color','White','LineWidth',6);
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    set([ax1,ax2],'Position',[.08 .08 .8 .8]);
    ax1.FontWeight = 'Bold'; 
    % title(ax1,['\sigma = ',num2str(minrho-1000),' - ',num2str(maxrho-1000)],'FontWeight','Bold','FontSize',80)
    set([ax1,ax2],'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
    % set(ax1,'XTick',[0 5 10.0 15 20.0 25.0 30.0],'XTickLabel',{['0',char(176),'E'],['5',char(176),'E'],['10',char(176),'E'],['15',char(176),'E'],['20',char(176),'E'],['25',char(176),'E'],['30',char(176),'E']})
    % set(ax1,'YTick',[-45 -40.0 -35.0  -30.0 -25 -20],'YTickLabel',{['45',char(176),'S'],['40',char(176),'S'],['35',char(176),'S'],['30',char(176),'S'],['25',char(176),'S'],['20',char(176),'S']})
    set(ax1,'XTick',[0 10.0 20.0 30.0],'XTickLabel',{['0',char(176),'E'],['10',char(176),'E'],['20',char(176),'E'],['30',char(176),'E']})
    set(ax1,'YTick',[-46.0 -35.0 -25],'YTickLabel',{['46',char(176),'S'],['35',char(176),'S'],['25',char(176),'S']})

    ax1.XAxis.FontSize = 100;
    ax1.YAxis.FontSize = 100;
    if t == 1
        % text(32.5,-21,'a)','FontWeight','Bold','FontSize',60)
    end
    if t == 2
        % text(32.5,-21,'b)','FontWeight','Bold','FontSize',60)
        yticklabels(ax2,[]);
        yticklabels(ax1,[]);
    end
    if t ==3 
        yticklabels(ax2,[]);
        yticklabels(ax1,[]);
        % text(32.5,-21,'c)','FontWeight','Bold','FontSize',60)
        % cb1 = colorbar(ax1,'northoutside','Position',[0.081369524720487,0.874375050709371,0.760846464059317,0.03]);
        cb2 = colorbar(ax2,'Position',[0.739417798029584,0.078756975360443,0.03,0.801582487710821]);
        % cb1.FontSize = 30;
        cb2.FontSize = 80;
        cb2.FontWeight = 'Bold';
        label = {sprintf('Potential Vorticity (s^{-1})')};
        ylabel(cb2,label,'FontWeight','Bold','FontSize',100); % Use the formatted string as the label
        % ylabel(cb2,['Average |DSC| ',string(minrho),' - ', string(maxrho),' (m^3 kg^{-1})'],'FontWeight','Bold','FontSize',44)
        % ylabel(cb1,'SSH (m)','FontWeight','Bold','FontSize',44)
    end
    xticklabels([ax1,ax2],[])
    label = sprintf('/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_PV_%d.svg', t);
    label = sprintf('/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/EGU_poster/Argo_PV_%d.svg', t);
    % saveas(gcf, label);

end
    


%% now plot salinity


for t =1:length(minrhos)
    minrho = minrhos(t);
    maxrho = maxrhos(t);
    indmin = findClosest(minrho,argo_density_grid.pot_rho);
    indmax = findClosest(maxrho,argo_density_grid.pot_rho);
    sal_avg_int = mean(abs(argo_density_grid.SA(indmin:indmax,:)),1,'omitnan');
 

    lats = -45:.5:-20;
    lons = 0:.5:35;
    
    [lats, lons] = meshgrid(lats,lons);
    lats = lats';
    lons = lons';
    
    store_profiles = {};
    
    ctr = 1;
    for i=1:size(lats,1)
        for j=1:size(lons,2)
            store_profiles.lat(ctr) = lats(i,1);
            store_profiles.lon(ctr) = lons(1,j);
            ind = find((argo_density_grid.lat >= (lats(i,1)-.25)) & (argo_density_grid.lat <= (lats(i,1)+.25)) & (argo_density_grid.lon >= (lons(1,j)-.25)) & (argo_density_grid.lon <= (lons(1,j)+.25)));
            temp_matix = sal_avg_int(ind);
            store_profiles.sal_avg_int(ctr) = mean(temp_matix,'omitnan');
            ctr = ctr +1;
            clear temp_matix
        end
    end
    
    
    stor_mean = reshape(store_profiles.sal_avg_int,[71,51])';
    
    figure
    set(gcf,'Position',[1 63 2392 1274],'color','w')
    ax1 = axes;
    ttt = pcolor(ssh_data.longitude,ssh_data.latitude,meanssh);
    set(ttt,'EdgeColor','none')
    colormap(ax1,cmocean('rain'))
    clim([-.25,1.5])
    axis(ax1,[minlon,maxlon,minlat,maxlat])
    ax2 = axes;
    hold on
    tttt = pcolor(lons,lats,(stor_mean));
    set(tttt,'EdgeColor','none')
    colormap(ax2,cmocean('haline'))
    if t == 1
        clim(ax2,[35.1,35.6])
        xLimits = ax2.XLim;
        yLimits = ax2.YLim;
        text(xLimits(1)-3.5, yLimits(2)-.5, 'c)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',100)

    end
    if t == 2
        clim(ax2,[34.6,35.1])
    end
    if t == 3
        clim(ax2,[34.3,34.8])
    end
    axis(ax2,[minlon,maxlon,minlat,maxlat])
    hold on
    contour(ssh_data.longitude,ssh_data.latitude,meanssh,[0.1 0.3 0.4 0.5 0.6 2],'k',"ShowText",true,'LineWidth',4);
    latind = find(ssh_data.latitude<-39.2);
    contour(ssh_data.longitude,ssh_data.latitude(latind),meanssh(latind,:),[0.3 10],'Color','Red','LineWidth',6);
    latind = find(ssh_data.latitude>-37.5);
    lonind = find(ssh_data.longitude<14.8);
    contour(ssh_data.longitude(lonind),ssh_data.latitude(latind),meanssh(latind,lonind),[0.5 10],'Color','Green','LineWidth',6);
    contour(ssh_data.longitude,ssh_data.latitude,meanssh,[.75 10],'Color','White','LineWidth',6);
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    set([ax1,ax2],'Position',[.08 .11 .8 .83]);
    ax1.FontWeight = 'Bold'; 
    set([ax1,ax2],'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
    % set(ax1,'XTick',[0 5 10.0 15 20.0 25.0 30.0],'XTickLabel',{['0',char(176),'E'],['5',char(176),'E'],['10',char(176),'E'],['15',char(176),'E'],['20',char(176),'E'],['25',char(176),'E'],['30',char(176),'E']})
    % set(ax1,'YTick',[-45 -40.0 -35.0  -30.0 -25 -20],'YTickLabel',{['45',char(176),'S'],['40',char(176),'S'],['35',char(176),'S'],['30',char(176),'S'],['25',char(176),'S'],['20',char(176),'S']})
    set(ax1,'XTick',[0 10.0 20.0 30.0],'XTickLabel',{['0',char(176),'E'],['10',char(176),'E'],['20',char(176),'E'],['30',char(176),'E']})
    set(ax1,'YTick',[-46.0 -35.0 -25],'YTickLabel',{['46',char(176),'S'],['35',char(176),'S'],['25',char(176),'S']})
    ax1.XAxis.FontSize = 100;
    ax1.YAxis.FontSize = 100;
    if t == 1
        % text(32.5,-21,'a)','FontWeight','Bold','FontSize',60)
    end
    if t == 2
        % text(32.5,-21,'b)','FontWeight','Bold','FontSize',60)
        yticklabels(ax2,[]);
        yticklabels(ax1,[]);
    end
    if t ==3 
        yticklabels(ax2,[]);
        yticklabels(ax1,[]);
        % text(32.5,-21,'c)','FontWeight','Bold','FontSize',60)
        % cb1 = colorbar(ax1,'northoutside','Position',[0.081369524720487,0.874375050709371,0.760846464059317,0.03]);
        cb2 = colorbar(ax2,'Position',[0.747779002042961,0.11015414961476,0.03,0.829406289945679]);
        % cb1.FontSize = 30;
        cb2.FontSize = 80;
        cb2.FontWeight = 'Bold';
        label = {sprintf('Salinity (g/kg)')};
        ylabel(cb2,label,'FontWeight','Bold','FontSize',100); % Use the formatted string as the label
        % ylabel(cb2,['Average |DSC| ',string(minrho),' - ', string(maxrho),' (m^3 kg^{-1})'],'FontWeight','Bold','FontSize',44)
        % ylabel(cb1,'SSH (m)','FontWeight','Bold','FontSize',44)
    end
    xticklabels([])
    label = sprintf('/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_sal_%d.svg', t);
    % saveas(gcf, label);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% T&S plots, geopotential height and DSC

store_max = nan(size(argo_depth_grid.CT,2),1);
store_len = nan(size(argo_depth_grid.CT,2),1);
for i = 1:size(argo_depth_grid.CT,2)
    ind = find(isfinite(argo_depth_grid.CT(:,i)));
    if ~isempty(ind)
        store_max(i) = max(argo_depth_grid.P(ind));
        store_len(i) = length(ind);
    end
end
shallowind = find(store_max<1000);
shortind = find(store_len<10);

AgulCBlat = [-24.7123, -30.7694, -36.5926, -39.6094,-36.4991,-31.4242];  % Example latitudes of polygon vertices
AgulCBlon = [34.9138 , 34.9415 , 26.2346 , 17.9714,6.9907,15.9472];  % Example longitudes of polygon vertices
AgulCBind = find(inpolygon(argo_density_grid.lon, argo_density_grid.lat, AgulCBlon, AgulCBlat));


argo_depth_grid.gpn_1000 = nan(size(argo_depth_grid.CT));
p_ref = 1000;
for i =1:size(argo_depth_grid.CT,2)
    if ~ismember(i,shallowind) && ~ismember(i,shortind) %&& ismember(i,AgulCBind)
        argo_depth_grid.gpn_1000(:,i) = gsw_geo_strf_dyn_height(argo_depth_grid.SA(:,i),argo_depth_grid.CT(:,i),argo_depth_grid.P',p_ref);
    end
    fprintf(1,'%d \n',i);
end

%% geopotential height heat map

store_gpn10 = argo_depth_grid.gpn_1000(11,:);
mygpn_all = store_gpn10./9.81;

minlon = 0;
maxlon = 35;
minlat = -45;
maxlat = -20;

lats = -45:.5:-20;
lons = 0:.5:35;

[lats, lons] = meshgrid(lats,lons);
lats = lats';
lons = lons';

store_profiles = {};

ctr = 1;
for i=1:size(lats,1)
    for j=1:size(lons,2)
        store_profiles.lat(ctr) = lats(i,1);
        store_profiles.lon(ctr) = lons(1,j);
        ind = find((argo_density_grid.lat >= (lats(i,1)-.25)) & (argo_density_grid.lat <= (lats(i,1)+.25)) & (argo_density_grid.lon >= (lons(1,j)-.25)) & (argo_density_grid.lon <= (lons(1,j)+.25)));
        temp_matix = mygpn_all(ind);
        store_profiles.mygpn_all(ctr) = mean(temp_matix,'omitnan');
        ctr = ctr +1;
        clear temp_matix
    end
end


stor_mean = reshape(store_profiles.mygpn_all,[71,51])';

figure
set(gcf,'Position',[1 63 2392 1274],'color','w')
ax1 = axes;
ttt = pcolor(ssh_data.longitude,ssh_data.latitude,meanssh);
set(ttt,'EdgeColor','none')
colormap(ax1,cmocean('rain'))
clim([-.25,1.5])
axis(ax1,[minlon,maxlon,minlat,maxlat])
ax2 = axes;
hold on
tttt = pcolor(lons,lats,(stor_mean));
set(tttt,'EdgeColor','none')
colormap(ax2,flipud(cmocean('matter')))
axis(ax2,[minlon,maxlon,minlat,maxlat])
clim(ax2,[.8,1.6])
hold on
contour(ssh_data.longitude,ssh_data.latitude,meanssh,[0.1 0.3 0.4 0.5 0.6 2],'k',"ShowText",true,'LineWidth',4);
latind = find(ssh_data.latitude<-39.2);
contour(ssh_data.longitude,ssh_data.latitude(latind),meanssh(latind,:),[0.3 10],'Color','Red','LineWidth',6);
latind = find(ssh_data.latitude>-37.5);
lonind = find(ssh_data.longitude<14.8);
contour(ssh_data.longitude(lonind),ssh_data.latitude(latind),meanssh(latind,lonind),[0.5 10],'Color','Green','LineWidth',6);
contour(ssh_data.longitude,ssh_data.latitude,meanssh,[.75 10],'Color','White','LineWidth',6);
linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
set([ax1,ax2],'Position',[.08 .11 .8 .83]);
ax1.FontWeight = 'Bold'; 
set([ax1,ax2],'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
set(ax1,'XTick',[0 10.0 20.0 30.0],'XTickLabel',{['0',char(176),'E'],['10',char(176),'E'],['20',char(176),'E'],['30',char(176),'E']})
set(ax1,'YTick',[-46.0 -35.0 -25],'YTickLabel',{['46',char(176),'S'],['35',char(176),'S'],['25',char(176),'S']})
ax1.XAxis.FontSize = 100;
ax1.YAxis.FontSize = 100;
yticklabels(ax2,[]);
yticklabels(ax1,[]);
cb2 = colorbar(ax2,'Position',[0.747779002042961,0.11015414961476,0.03,0.829406289945679]);
cb2.FontSize = 80;
cb2.FontWeight = 'Bold';
label = {sprintf('Geopotential Height (m)')};
ylabel(cb2,label,'FontWeight','Bold','FontSize',100); % Use the formatted string as the label
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-1, yLimits(2)-.5, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',100)


label = sprintf('/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_dyn_height.svg');
% saveas(gcf, label);



%% EKE
load('Agulhas_CapeBasin_ssh_30yr.mat')

meanu = mean(SSH_data.ugos,3);
meanv = mean(SSH_data.vgos,3);


uprime = SSH_data.ugos - meanu;
vprime = SSH_data.vgos - meanv;

uprime2 = uprime.^2;
vprime2 = vprime.^2;

uprime_m =  mean(uprime2,3);
vprime_m =  mean(vprime2,3);

eke = .5*(uprime_m+vprime_m);

figure
set(gcf,'Position',[1 63 2392 1274],'color','w')
ax1 = axes;
ttt = pcolor(ssh_data.longitude,ssh_data.latitude,meanssh);
set(ttt,'EdgeColor','none')
colormap(ax1,cmocean('rain'))
clim([-.25,1.5])
axis(ax1,[minlon,maxlon,minlat,maxlat])
ax2 = axes;
hold on
tttt = pcolor(SSH_data.lon,SSH_data.lat,eke');
set(tttt,'EdgeColor','none')
colormap(ax2,cmocean('speed'))
axis(ax2,[minlon,maxlon,minlat,maxlat])
clim(ax2,[0,.3])
hold on
contour(ssh_data.longitude,ssh_data.latitude,meanssh,[0.1 0.3 0.4 0.5 0.6 2],'k',"ShowText",true,'LineWidth',4);
latind = find(ssh_data.latitude<-39.2);
contour(ssh_data.longitude,ssh_data.latitude(latind),meanssh(latind,:),[0.3 10],'Color','Red','LineWidth',6);
latind = find(ssh_data.latitude>-37.5);
lonind = find(ssh_data.longitude<14.8);
contour(ssh_data.longitude(lonind),ssh_data.latitude(latind),meanssh(latind,lonind),[0.5 10],'Color','Green','LineWidth',6);
contour(ssh_data.longitude,ssh_data.latitude,meanssh,[.75 10],'Color','White','LineWidth',6);
linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
set([ax1,ax2],'Position',[.08 .11 .8 .83]);
ax1.FontWeight = 'Bold'; 
set([ax1,ax2],'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
set(ax1,'XTick',[0 10.0 20.0 30.0],'XTickLabel',{['0',char(176),'E'],['10',char(176),'E'],['20',char(176),'E'],['30',char(176),'E']})
set(ax1,'YTick',[-46.0 -35.0 -25],'YTickLabel',{['46',char(176),'S'],['35',char(176),'S'],['25',char(176),'S']})
ax1.XAxis.FontSize = 100;
ax1.YAxis.FontSize = 100;
yticklabels(ax2,[]);
cb2 = colorbar(ax2,'Position',[0.747779002042961,0.11015414961476,0.03,0.829406289945679]);
cb2.FontSize = 80;
cb2.FontWeight = 'Bold';
label = {sprintf('EKE (m^2/s^2)')};
ylabel(cb2,label,'FontWeight','Bold','FontSize',100); % Use the formatted string as the label
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-3.5, yLimits(2)-.5, 'd)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',100)

label = sprintf('/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_EKE.svg');
% saveas(gcf, label);

%% just label lines


figure
set(gcf,'Position',[1 63 2392 1274],'color','w')
ax1 = axes;
ttt = pcolor(ssh_data.longitude,ssh_data.latitude,meanssh);
set(ttt,'EdgeColor','none')
colormap(ax1,cmocean('rain'))
clim([-.25,1.5])
axis(ax1,[minlon,maxlon,minlat,maxlat])
hold on
contour(ssh_data.longitude,ssh_data.latitude,meanssh,[0.1 0.3 0.4 0.5 0.6 2],'k',"ShowText",true,'LineWidth',6);
latind = find(ssh_data.latitude<-39.2);
contour(ssh_data.longitude,ssh_data.latitude(latind),meanssh(latind,:),[0.3 10],'Color','Red','LineWidth',14);
latind = find(ssh_data.latitude>-37.5);
lonind = find(ssh_data.longitude<14.8);
contour(ssh_data.longitude(lonind),ssh_data.latitude(latind),meanssh(latind,lonind),[0.5 10],'Color','Green','LineWidth',14);
contour(ssh_data.longitude,ssh_data.latitude,meanssh,[.75 10],'Color','White','LineWidth',14);
linkaxes(ax1)
set(ax1,'Position',[.08 .11 .8 .83]);
ax1.FontWeight = 'Bold'; 
set([ax1,ax2],'DataAspectRatio',[1/cosd(mean([minlat maxlat])) 1 1]  );
set(ax1,'XTick',[0 10.0 20.0 30.0],'XTickLabel',{['0',char(176),'E'],['10',char(176),'E'],['20',char(176),'E'],['30',char(176),'E']})
set(ax1,'YTick',[-46.0 -35.0 -25],'YTickLabel',{['46',char(176),'S'],['35',char(176),'S'],['25',char(176),'S']})
ax1.XAxis.FontSize = 100;
ax1.YAxis.FontSize = 100;
cb1 = colorbar(ax1,'Position',[0.747779002042961,0.11015414961476,0.03,0.829406289945679]);
cb1.FontSize = 80;
cb1.FontWeight = 'Bold';
label = {sprintf('SSH (m)')};
ylabel(cb1,label,'FontWeight','Bold','FontSize',100); % Use the formatted string as the label
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-3.5, yLimits(2)-.5, 'e)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',100)

label = sprintf('/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/Argo_fronts_SSH.svg');
% saveas(gcf, label);


