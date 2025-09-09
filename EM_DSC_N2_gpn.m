clear
close all

addpath(genpath("/Users/lindsay.grose/Documents/URI/git.nosync/"))
addpath(genpath("/Users/lindsay.grose/Documents/MATLAB/"))
%% load float data
load('/Users/lindsay.grose/Documents/MATLAB/Datasets/EM_Apex_Floats/Processes_data/Density/float_data_all.mat')

load('/Users/lindsay.grose/Documents/MATLAB/Datasets/EM_Apex_Floats/Processes_data/Depth/float_data_all.mat')

% minrho = 1026;
% maxrho = 1026.5;

colorem = rgb('DarkSlateGray');

%% geopotential height

p_ref = 975;
depthss = meshgrid(float_depth_data.depth,float_depth_data.absolute_salinity(1,:))';
dyn_height = gsw_geo_strf_dyn_height(float_depth_data.absolute_salinity, float_depth_data.conservative_temperature, depthss, p_ref);


%integrated value
int_dyn_height_100_float = dyn_height(51,:)./9.81;

%% buoyancy frequency
depthsss = meshgrid(float_depth_data.depth,float_depth_data.absolute_salinity(1,:))';
[N2,pmid] = gsw_Nsquared(float_depth_data.absolute_salinity(51:end,:), float_depth_data.conservative_temperature(51:end,:), depthsss(51:end,:));

%% DSC

windowSize = 41; % Window size (odd number)
blackmanWindow = blackman(windowSize);
for i=1:size(float_density_data.conservative_temperature,2)
    float_density_data.conservative_temperature_filt(:,i) = conv(float_density_data.conservative_temperature(:,i), blackmanWindow, 'same') / sum(blackmanWindow);
end
for i=1:size(float_density_data.absolute_salinity,2)
    float_density_data.absolute_salinity_filt(:,i) = conv(float_density_data.absolute_salinity(:,i), blackmanWindow, 'same') / sum(blackmanWindow);
end



%calculate DSC
tau = gsw_spiciness0(float_density_data.absolute_salinity_filt,float_density_data.conservative_temperature_filt);
dtau = diff(tau,1,1);
dpot_rho = diff(float_density_data.potential_density,1,1);
dtau_drho = dtau./dpot_rho;
dtau2 = diff(dtau_drho,1,1);
spice_curve = dtau2./dpot_rho(1:end-1,:);


blahtemp = nan(size(spice_curve,1)+2,size(spice_curve,2));    
blahtemp(2:end-1,:) = spice_curve;
float_density_data.DSC = blahtemp;


%% alternative figure with just N2

minrho = 1026.0;
maxrho = 1027.0;

ind1026 = findClosest(minrho,float_density_data.potential_density);
ind10275 = findClosest(maxrho,float_density_data.potential_density);

meanDSC_float = nan(size(float_density_data.pressure,2),1);
meantemp = nan(size(float_density_data.pressure,2),1);

for i=1:size(float_density_data.pressure,2)
    ind100 = findClosest(100,float_density_data.pressure(:,i));

    if ind100 > ind1026 && ind100 <= ind10275
        meanDSC_float(i) = mean(abs(float_density_data.DSC(ind100:ind10275,i)),'omitmissing');
        meantemp(i) = mean(abs(float_density_data.conservative_temperature(ind100:ind10275,i)),'omitmissing');
    end
    if ind100 < ind1026
        meanDSC_float(i) = mean(abs(float_density_data.DSC(ind1026:ind10275,i)),'omitmissing');
        meantemp(i) = mean(abs(float_density_data.conservative_temperature(ind1026:ind10275,i)),'omitmissing');
    end
end


meanN2_float = nan(size(meanDSC_float));
gpn_float = nan(size(meanDSC_float));

for i=1:length(meanDSC_float)
    ind1026 = findClosest(minrho,float_depth_data.potential_density(:,i));
    ind10275 = findClosest(maxrho,float_depth_data.potential_density(:,i));
    % meanN2_float(i) = mean(real(sqrt(N2(ind1026:ind10275,i))),'omitmissing');
    meanN2_float(i) = mean(N2(ind1026:ind10275,i),'omitmissing');
    if ind1026 >=51
        gpn_float(i) = dyn_height(ind1026,i)./9.81 - dyn_height(ind10275,i)./9.81;
    else
        gpn_float(i) = dyn_height(51,i)./9.81 - dyn_height(ind10275,i)./9.81;
    end
end

CBlat = [-35.9265, -38.7206, -39.6029, -41.8088,-40.2794,-28.5441,-31.1324];  % Example latitudes of polygon vertices
CBlon = [ 24.2305 , 20.2550 , 18.4067  , 16.6630,5.8175,9.9674,17.4651];  % Example longitudes of polygon vertices
CBind = find(inpolygon(float_depth_data.lon, float_depth_data.lat, CBlon, CBlat));


figure
set(gcf,'Position',[38 550 1888 787],'color','w')
% sgtitle([num2str(minrho),' - ',num2str(maxrho),' kg/m^3'],'FontWeight','Bold','FontSize',80)
ax=subplot(1,2,1);
% plot(meanN2_float(CBind),meanDSC_float(CBind),'.','MarkerSize',50,'Color',rgb('DodgerBlue'),'HandleVisibility','off')
scatter(meanN2_float(CBind),meanDSC_float(CBind),150,int_dyn_height_100_float(CBind),'filled','HandleVisibility','off')
colormap(ax,flipud(cmocean('matter')))
clim(ax,[.7,1.4])
hold on
ax.XAxis.FontSize = 28;
ax.YAxis.FontSize = 28;
xlabel('N^2 (1/s^2)','FontWeight','Bold','FontSize',50)
ylabel('|DSC| (m^3/kg)','FontWeight','Bold','FontSize',50)
xLimits = ax.XLim;
yLimits = ax.YLim;
text(xLimits(1)-.7e-5, yLimits(2)+10, 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)
ylim([0,150])
xlim([.6e-5,4e-5])
title('Mode and Thermocline','FontWeight','Bold','FontSize',50)

% intermediate 
minrho = 1027.0;
maxrho = 1027.5;

ind1026 = findClosest(minrho,float_density_data.potential_density);
ind10275 = findClosest(maxrho,float_density_data.potential_density);

meanDSC_float = nan(size(float_density_data.pressure,2),1);
meantemp = nan(size(float_density_data.pressure,2),1);

for i=1:size(float_density_data.pressure,2)
    ind100 = findClosest(100,float_density_data.pressure(:,i));

    if ind100 > ind1026 && ind100 <= ind10275
        meanDSC_float(i) = mean(abs(float_density_data.DSC(ind100:ind10275,i)),'omitmissing');
        meantemp(i) = mean(abs(float_density_data.conservative_temperature(ind100:ind10275,i)),'omitmissing');
    end
    if ind100 < ind1026
        meanDSC_float(i) = mean(abs(float_density_data.DSC(ind1026:ind10275,i)),'omitmissing');
        meantemp(i) = mean(abs(float_density_data.conservative_temperature(ind1026:ind10275,i)),'omitmissing');
    end
end


meanN2_float = nan(size(meanDSC_float));
gpn_float = nan(size(meanDSC_float));

for i=1:length(meanDSC_float)
    ind1026 = findClosest(minrho,float_depth_data.potential_density(:,i));
    ind10275 = findClosest(maxrho,float_depth_data.potential_density(:,i));
    % meanN2_float(i) = mean(real(sqrt(N2(ind1026:ind10275,i))),'omitmissing');
    meanN2_float(i) = mean(N2(ind1026:ind10275,i),'omitmissing');
    if ind1026 >=51
        gpn_float(i) = dyn_height(ind1026,i)./9.81 - dyn_height(ind10275,i)./9.81;
    else
        gpn_float(i) = dyn_height(51,i)./9.81 - dyn_height(ind10275,i)./9.81;
    end
end

CBlat = [-35.9265, -38.7206, -39.6029, -41.8088,-40.2794,-28.5441,-31.1324];  % Example latitudes of polygon vertices
CBlon = [ 24.2305 , 20.2550 , 18.4067  , 16.6630,5.8175,9.9674,17.4651];  % Example longitudes of polygon vertices
CBind = find(inpolygon(float_depth_data.lon, float_depth_data.lat, CBlon, CBlat));

ax=subplot(1,2,2);
% plot(meanN2_float(CBind),meanDSC_float(CBind),'.','MarkerSize',50,'Color',rgb('DodgerBlue'),'HandleVisibility','off')
scatter(meanN2_float(CBind),meanDSC_float(CBind),150,int_dyn_height_100_float(CBind),'filled','HandleVisibility','off')
colormap(ax,flipud(cmocean('matter')))
clim(ax,[.7,1.4])
cb1 = colorbar(ax);
hold on
ax.XAxis.FontSize = 28;
ax.YAxis.FontSize = 28;
cb1.FontSize = 28;
xlabel('N^2 (1/s^2)','FontWeight','Bold','FontSize',50)
ylabel('|DSC| (m^3/kg)','FontWeight','Bold','FontSize',50)
ylabel(cb1,'Geopotential height (m)','FontWeight','Bold','FontSize',50)
xLimits = ax.XLim;
yLimits = ax.YLim;
text(xLimits(1)-.9e-6, yLimits(2), 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)
ylim([0,150])
xlim([5.8e-6,11.2e-6])
title('Intermediate','FontWeight','Bold','FontSize',50)

% saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/EM_DSC_N2.svg');
