function plotNodeHeatmap(FN, Ephys, channels, maxVal, Params, coords, ...
    metricVarName, metricLabel, cmap, figFolder, figName, oneFigureHandle, subsetChannelName)
%PLOTNODEHEATMAP Summary of this function goes here
%   Detailed explanation goes here

% subsetChannelName : not impelemented yet, but the idea is to allow
% plotting a subset of available channels
%% plot
p = [50 100 1150 570];

if Params.showOneFig
    if isgraphics(oneFigureHandle)
        set(oneFigureHandle, 'OuterPosition', p);
    else 
        oneFigureHandle = figure;
        set(oneFigureHandle, 'OuterPosition', p);
    end 
else 
    F1 = figure;
    F1.OuterPosition = p;
end 

tiledlayout(1,2)
aesthetics; axis off; hold on

%% coordinates

% Perform transpose if not column vector 
if size(channels, 1) == 1
    channels = channels'; 
end 

xc = coords(:, 1);
yc = coords(:, 2); 

%% Assign metric variable 

metricVals = Ephys.(metricVarName);

%% plot electrodes
% TODO: I think a lot of rectangle coloring can be simplified
numCbarTicks = 5;


%% Left electrode plot (scaled to individual recording)
nexttile
uniqueXc = sort(unique(xc));
nodeScaleF = 2/3; 

makeMinSpikeCountZero = 1;

if makeMinSpikeCountZero == 1
    minSpikeCountToPlot = 0;
else 
    minSpikeCountToPlot = min(spikeCount);
end 

numChannels = length(channels);
for i = 1:numChannels

    pos = [xc(i)-(0.5*nodeScaleF) yc(i)-(0.5*nodeScaleF) nodeScaleF nodeScaleF];
        try
            colorToUse = cmap(ceil(length(cmap) * ((metricVals(i) - minSpikeCountToPlot)/(prctile(metricVals,99,'all')-minSpikeCountToPlot))),1:3);
            rectangle('Position',pos,'Curvature',[1 1],'FaceColor',colorToUse,'EdgeColor','w','LineWidth',0.1) 
        catch
            if (metricVals(i) - minSpikeCountToPlot) / (prctile(metricVals,95,'all') - minSpikeCountToPlot) == 0
                rectangle('Position',pos,'Curvature',[1 1],'FaceColor', ...
                    cmap(ceil(length(cmap)*((metricVals(i)- minSpikeCountToPlot)/(prctile(metricVals,99,'all')-minSpikeCountToPlot))+0.00001),1:3),'EdgeColor','w','LineWidth',0.1)
            else
                rectangle('Position',pos,'Curvature',[1 1],'FaceColor',cmap(length(cmap),1:3),'EdgeColor','w','LineWidth',0.1) 
            end
        end
    if Params.includeChannelNumberInPlots 
        text(pos(1) + 0.5 * nodeScaleF, pos(2) + 0.5 * nodeScaleF, ...
            sprintf('%.f', channels(i)), 'HorizontalAlignment','center')
    end 
end
ylim([min(yc) - 1, max(yc) + 1])
xlim([min(xc) - 1, max(xc) + 1])
axis off

cb = colorbar;
cb.Box = 'off';
cb.Ticks = linspace(0, 1, numCbarTicks);

tickLabels = cell(numCbarTicks, 1);
for nTick = 1:numCbarTicks
    if nTick == 1
        tickLabels{nTick} = num2str(minSpikeCountToPlot);
    else 
        tickLabels{nTick} = num2str(round((nTick-1) / numCbarTicks * prctile(metricVals,99,'all'),2));
    end 
end 

cb.TickLabels = tickLabels;


cb.TickDirection = 'out';
cb.Label.String = metricLabel;
title({strcat(regexprep(FN,'_','','emptymatch'),' Electrode heatmap scaled to recording'),' '});

%% Right electrode plot (scaled to all recordings)

nexttile
for i = 1:numChannels
    pos = [xc(i)-(0.5*nodeScaleF) yc(i)-(0.5*nodeScaleF) nodeScaleF nodeScaleF];
        try
            rectangle('Position', pos, 'Curvature', [1 1], 'FaceColor', ...
                cmap(ceil(length(cmap)*((metricVals(i) - minSpikeCountToPlot) / (maxVal-minSpikeCountToPlot))),1:3),'EdgeColor','w','LineWidth',0.1)
        catch
            if (metricVals(i)-minSpikeCountToPlot)/(maxVal - minSpikeCountToPlot) == 0
                rectangle('Position',pos,'Curvature',[1 1],'FaceColor', ...
                    cmap(ceil(length(cmap)*((metricVals(i) - minSpikeCountToPlot) / (maxVal-minSpikeCountToPlot))+0.00001),1:3),'EdgeColor','w','LineWidth',0.1)
            else
                 rectangle('Position',pos,'Curvature',[1 1],'FaceColor',cmap(length(cmap),1:3),'EdgeColor','w','LineWidth',0.1) 
            end
        end
end
ylim([min(yc)-1 max(yc)+1])
xlim([min(xc)-1 max(xc)+1])
axis off

cb = colorbar;
cb.Box = 'off';
cb.Ticks = linspace(0, 1, numCbarTicks);

tickLabels = cell(numCbarTicks, 1);
for nTick = 1:numCbarTicks
    if nTick == 1
        tickLabels{nTick} = num2str(minSpikeCountToPlot);
    else 
        tickLabels{nTick} = num2str(round((nTick-1) / numCbarTicks * maxVal, 2));
    end 
end 

cb.TickLabels = tickLabels;
cb.TickDirection = 'out';
cb.Label.String = metricLabel;
title({strcat(regexprep(FN,'_','','emptymatch'),' Electrode heatmap scaled to entire data batch'),' '});

% save figure
figPath = fullfile(figFolder, figName);

if Params.showOneFig
    pipelineSaveFig(figPath, Params.figExt, Params.fullSVG, oneFigureHandle);
else
    pipelineSaveFig(figPath, Params.figExt, Params.fullSVG)
end 

if Params.showOneFig
    clf(oneFigureHandle, 'reset')
else 
    close(F1);
end 


end

