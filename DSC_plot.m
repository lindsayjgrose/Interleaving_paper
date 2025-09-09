clear
close all

%%
load('pot_density_grid.mat')
load('/Users/lindsay.grose/Documents/MATLAB/Datasets/WireFlyer/QUICCHE/Dives/Processed/20230318_070616/20230318_070616_derived_quantities.mat');


%% plot
minrho = 26.35;
maxrho = 27.1;

mean_DSC = mean(abs(data.density.DSC),1,'omitmissing');

minind = find(min(mean_DSC)==mean_DSC);
% maxind = find(max(mean_DSC)==mean_DSC);

maxind = findClosest(80,data.density.dist(1,:));

C = contour(data.density.dist,data.density.potential_density,data.density.DSC,[-70 70]);

x = data.density.dist(:,310:end);                  % size [1 x N]
y = data.density.potential_density(:,310:end);     % size [1 x M]
z = data.density.DSC(:,310:end);                        % size [M x N], matches meshgrid(y,x)

figure
set(gcf,'Position',[1 159 1726 1177],'color','w')

% ax2 = subplot(2,4,2:3);
ax2 = subplot(2,4,2:3);
ttt = pcolor(data.density.dist,data.density.potential_density,data.density.DSC);
set(ttt,'EdgeColor','none')
set(ax2,'YDir','reverse');
colormap(ax2,(cmocean('balance')))
cb = colorbar(ax2,'Position',[0.806363479473117,0.549702633814784,0.01808611496489,0.399807756893538]);
ax2.FontSize = 24;
hold on
xline(data.density.dist(1,minind),'LineWidth',5)
xline(data.density.dist(1,maxind),'LineWidth',5)
minPoints = 500;
% Parse the contour matrix
i = 1;
while i < size(C, 2)
    level = C(1, i);
    numPoints = C(2, i);
    contourSegment = C(:, i+1:i+numPoints);
    
    % Create a polygon from the contour and count points inside
    in = inpolygon(x(:), y(:), contourSegment(1,:), contourSegment(2,:));
    numInside = sum(in);
    
    % Plot only if it passes the threshold
    if numInside >= minPoints
        plot(contourSegment(1,:), contourSegment(2,:), 'Color', rgb('Cyan'), 'LineWidth', 3);
    end
    
    i = i + numPoints + 1;
end
xlabel('Distance (km)','FontWeight','Bold','FontSize',36)
ylabel('\sigma (kg/m^3)','FontWeight','Bold','FontSize',44)
ylabel(cb,'m^3/kg','FontWeight','Bold','FontSize',36)
title('Spiciness Curvature','FontWeight','Bold','FontSize',44)
clim([-150,150])
ylim(ax2,[minrho,maxrho])
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-10, yLimits(2)-.7, 'a)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)
pos = get(ax2, 'Position');
pos(3) = 0.6;  % Set the width of all subplots to be equal
pos(4) = 0.4;  % Set the height of all subplots to be equal
set(ax2, 'Position', pos+[-.14,-.035,0,0]);
xlim(ax2,[5,106])


% ax1 = subplot(2,4,5:6);
ax1 = subplot(2,4,[5:6]);
plot(data.density.SA(:,minind),data.density.CT(:,minind),'k');
hold on;
scatter(data.density.SA(:,minind),data.density.CT(:,minind),140,data.density.DSC(:,minind), 'filled');
colormap(ax1,(cmocean('balance'))); %flipud(slanCM('RdYlBu')))
clim(ax1,[-150,150])
ax1.XAxis.FontSize = 34;
ax1.YAxis.FontSize = 34;
hold on
[contour_handle,~] = contour(pot_density_grid.sal,pot_density_grid.temp,pot_density_grid.pot_rho-1000,[24,24.5,25,25.5,26,26.25,26.5,26.75,27,27.25,27.5,27.75,28,28.5],'k','HandleVisibility','off','color',rgb('DarkGray'));
labelhand = clabel(contour_handle,'FontSize',30);
axis(ax1,[34.55,35.1,6.8,11.6])
% axis(ax1,[34.5,35.4,6,14])
% for k = 1:length(labelhand)
%     if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(27))% Check if it's a text object
%         labelhand(k).Position(1) = labelhand(k).Position(1) + 0.5; % Shift left
%         labelhand(k).Position(2) = labelhand(k).Position(2) + 1.5; % Shift left
%     end
% end
ylabel(['Temperature (',char(176),'C)'],'FontWeight','Bold','FontSize',44);
xlabel("Salinity (g/kg)",'FontWeight','Bold','FontSize',44);
% plot(meanprof_raw.SA_smooth,meanprof_raw.CT_smooth,'color',rgb('Cyan'),'linewidth',4)
xLimits = ax1.XLim;
yLimits = ax1.YLim;
text(xLimits(1)-.1, yLimits(2)-.5, 'b)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)
pos = get(ax1, 'Position');
pos(3) = 0.36;  % Set the width of all subplots to be equal
pos(4) = 0.38;  % Set the height of all subplots to be equal
set(ax1, 'Position', pos+[-0.02,-.02,0,0]);

% ax2 = subplot(2,4,7:8);
ax2 = subplot(2,4,[7:8]);
plot(data.density.SA(:,maxind),data.density.CT(:,maxind),'k');
hold on;
scatter(data.density.SA(:,maxind),data.density.CT(:,maxind),140,data.density.DSC(:,maxind), 'filled');
colormap(ax2,(cmocean('balance'))); %flipud(slanCM('RdYlBu')))
clim(ax2,[-150,150])
ax2.XAxis.FontSize = 34;
ax2.YAxis.FontSize = 34;
% cb1 = colorbar;
% cb1.FontSize = 30;
% cb1.FontWeight = 'Bold';
% cb1 = colorbar(ax2,'Position',[0.918762088974855,0.129371851712021,0.019342359767892,0.819822695995538]);
% ylabel(cb1,'Diapycnal Spiciness Curvature (m^3/kg)','FontWeight','Bold','FontSize',40); % Use the formatted string as the label
hold on
[contour_handle,~] = contour(pot_density_grid.sal,pot_density_grid.temp,pot_density_grid.pot_rho-1000,[24,24.5,25,25.5,26,26.25,26.5,26.75,27,27.25,27.5,27.75,28,28.5],'k','HandleVisibility','off','color',rgb('DarkGray'));
labelhand = clabel(contour_handle,'FontSize',30);
% labelhand = clabel(contour_handle,h,'FontSize',30,'LabelSpacing', 900);
axis(ax2,[34.75,35.35,8.5,13.5])
% axis(ax2,[34.5,35.4,6,14])
for k = 1:length(labelhand)
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(27))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) + 0.6; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) + 2.9; % Shift left
    end
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(26.5))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) + 0.02; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) - 0; % Shift left
    end    
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(26))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) + 0.01; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) - .2; % Shift left
    end        
end
% ylabel(['Temperature (',char(176),'C)'],'FontWeight','Bold','FontSize',44);
xlabel("Salinity (g/kg)",'FontWeight','Bold','FontSize',44);
% plot(meanprof_raw.SA_smooth,meanprof_raw.CT_smooth,'color',rgb('Cyan'),'linewidth',4)
xLimits = ax2.XLim;
yLimits = ax2.YLim;
text(xLimits(1)-.06, yLimits(2)-.5, 'c)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontWeight','Bold','FontSize',50)
pos = get(ax2, 'Position');
pos(3) = 0.36;  % Set the width of all subplots to be equal
pos(4) = 0.38;  % Set the height of all subplots to be equal
set(ax2, 'Position', pos+[0.01,-.02,0,0]);



saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/explain_DSC_plot.svg');


%% insets for main figure

minrho = 26.35;
maxrho = 27.1;

mean_DSC = mean(abs(data.density.DSC),1,'omitmissing');

minind = find(min(mean_DSC)==mean_DSC);
% maxind = find(max(mean_DSC)==mean_DSC);

maxind = findClosest(80,data.density.dist(1,:));


figure
set(gcf,'Position',[1 159 1726 1177],'color','w')

% ax2 = subplot(2,4,7:8);
ax2 = subplot(2,4,[7:8]);
plot(data.density.SA(:,maxind),data.density.CT(:,maxind),'k');
hold on;
scatter(data.density.SA(:,maxind),data.density.CT(:,maxind),400,data.density.DSC(:,maxind), 'filled');
colormap(ax2,(cmocean('balance'))); %flipud(slanCM('RdYlBu')))
clim(ax2,[-150,150])
ax2.XAxis.FontSize = 34;
ax2.YAxis.FontSize = 34;
hold on
[contour_handle,~] = contour(pot_density_grid.sal,pot_density_grid.temp,pot_density_grid.pot_rho-1000,[24,24.5,25,25.5,26,26.25,26.5,26.75,27,27.25,27.5,27.75,28,28.5],'k','HandleVisibility','off','color',rgb('DarkGray'));
labelhand = clabel(contour_handle,'FontSize',30);
axis(ax2,[35.07,35.12,11.3,11.6])
yticklabels([])
xticklabels([])
for k = 1:length(labelhand)
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(27))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) + 0.5; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) + 2.5; % Shift left
    end
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(26.5))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) - 0.2; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) - .9; % Shift left
    end    
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(26))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) + 0.01; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) - .2; % Shift left
    end        
end
xLimits = ax2.XLim;
yLimits = ax2.YLim;
pos = get(ax2, 'Position');
pos(3) = 0.36;  % Set the width of all subplots to be equal
pos(4) = 0.38;  % Set the height of all subplots to be equal
set(ax2, 'Position', pos+[0.01,-.02,0,0]);

saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/explain_DSC_window1.svg');


figure
set(gcf,'Position',[1 159 1726 1177],'color','w')

ax2 = subplot(2,4,[7:8]);
plot(data.density.SA(:,maxind),data.density.CT(:,maxind),'k');
hold on;
scatter(data.density.SA(:,maxind),data.density.CT(:,maxind),400,data.density.DSC(:,maxind), 'filled');
colormap(ax2,(cmocean('balance'))); %flipud(slanCM('RdYlBu')))
clim(ax2,[-150,150])
ax2.XAxis.FontSize = 34;
ax2.YAxis.FontSize = 34;
hold on
[contour_handle,~] = contour(pot_density_grid.sal,pot_density_grid.temp,pot_density_grid.pot_rho-1000,[24,24.5,25,25.5,26,26.25,26.5,26.75,27,27.25,27.5,27.75,28,28.5],'k','HandleVisibility','off','color',rgb('DarkGray'));
labelhand = clabel(contour_handle,'FontSize',30);
axis(ax2,[34.88,34.98,9.9,10.5])
yticklabels([])
xticklabels([])
for k = 1:length(labelhand)
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(27))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) + 0.5; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) + 2.5; % Shift left
    end
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(26.5))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) - 0.2; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) - .9; % Shift left
    end    
    if isa(labelhand(k), 'matlab.graphics.primitive.Text') && (labelhand(k).String == string(26))% Check if it's a text object
        labelhand(k).Position(1) = labelhand(k).Position(1) + 0.01; % Shift left
        labelhand(k).Position(2) = labelhand(k).Position(2) - .2; % Shift left
    end        
end
xLimits = ax2.XLim;
yLimits = ax2.YLim;
pos = get(ax2, 'Position');
pos(3) = 0.36;  % Set the width of all subplots to be equal
pos(4) = 0.38;  % Set the height of all subplots to be equal
set(ax2, 'Position', pos+[0.01,-.02,0,0]);


saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/explain_DSC_window2.svg');



C = contour(data.density.dist,data.density.potential_density,data.density.DSC,[-70 70]);

x = data.density.dist(:,310:end);                  % size [1 x N]
y = data.density.potential_density(:,310:end);     % size [1 x M]
z = data.density.DSC(:,310:end);                        % size [M x N], matches meshgrid(y,x)

figure
set(gcf,'Position',[1 159 1726 1177],'color','w')

% ax2 = subplot(2,4,2:3);
ax2 = subplot(2,4,2:3);
ttt = pcolor(data.density.dist,data.density.potential_density,data.density.DSC);
set(ttt,'EdgeColor','none')
set(ax2,'YDir','reverse');
colormap(ax2,(cmocean('balance')))
hold on
xline(data.density.dist(1,minind),'LineWidth',5)
xline(data.density.dist(1,maxind),'LineWidth',5)
minPoints = 500;
% Parse the contour matrix
i = 1;
while i < size(C, 2)
    level = C(1, i);
    numPoints = C(2, i);
    contourSegment = C(:, i+1:i+numPoints);
    
    % Create a polygon from the contour and count points inside
    in = inpolygon(x(:), y(:), contourSegment(1,:), contourSegment(2,:));
    numInside = sum(in);
    
    % Plot only if it passes the threshold
    if numInside >= minPoints
        scatter(contourSegment(1,:), contourSegment(2,:),80, rgb('Cyan'), 'filled');
    end
    
    i = i + numPoints + 1;
end
clim([-150,150])
ylim(ax2,[minrho+.01,26.6])
xlim(ax2,[94.5,104])
daspect([45 1 1])
yticklabels([])
xticklabels([])

saveas(gcf, '/Users/lindsay.grose/Documents/MATLAB/Saved_plots/QUICCHE/Interleaving_paper/explain_DSC_window3.svg');



