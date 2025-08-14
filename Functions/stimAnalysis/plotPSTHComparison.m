function plotPSTHComparison(allNetworkResponses, outputDir)
% PLOTPOSTSTIMCOMPARISON - Generates comparison plots for PSTH metrics across conditions.
%
% Syntax:
%   plotPSTHComparison(allNetworkResponses, outputDir)
%
% Description:
%   This function uses HalfViolinPlot to maintain visual consistency with
%   the rest of the MEANAP pipeline.
%
% Inputs:
%   allNetworkResponses - A cell array where each cell contains the
%                         'networkResponse' structure from a single recording.
%
%   outputDir           - The path to the top-level output directory. Figures
%                         will be saved in a '2A_IndividualNeuronalAnalysis' subfolder.

    if isempty(allNetworkResponses)
        disp('No network response data provided for comparison plots.');
        return;
    end
    
    figFolder = fullfile(outputDir, '2A_IndividualNeuronalAnalysis');
    if ~exist(figFolder, 'dir')
        mkdir(figFolder);
    end

    % Pool data from all recordings
    pooledData = vertcat(allNetworkResponses{:});
    
    if ~isfield(pooledData, 'condition')
        disp('No "condition" field found in networkResponse data. Skipping comparison plots.');
        return;
    end

    % Extract unique conditions
    conditions = unique({pooledData.condition});
    numConditions = length(conditions);

    % --- Figure 1: Response Magnitude (AUC) ---
    fig1 = figure('Position', [100, 100, 800, 600], 'Color', 'w');
    sgtitle('Response Magnitude Comparison by Condition');
    
    % Subplot 1: Response AUC
    ax1 = subplot(1, 2, 1); hold on;
    title('Response AUC');
    ylabel('AUC');
    
    % Subplot 2: Baseline-Corrected AUC
    ax2 = subplot(1, 2, 2); hold on;
    title('Baseline-Corrected AUC');
    ylabel('AUC');

    for i = 1:numConditions
        conditionData = pooledData(strcmp({pooledData.condition}, conditions{i}));
        color = lines(numConditions);

        HalfViolinPlot([conditionData.auc_response], i, color(i,:), 0.4, 'auto');
        HalfViolinPlot([conditionData.auc_corrected], i, color(i,:), 0.4, 'auto');
    end
    set([ax1, ax2], 'xtick', 1:numConditions, 'xticklabel', conditions);
    
    % --- Figure 2: Response Duration (Half-Max Time) ---
    fig2 = figure('Position', [150, 150, 600, 500], 'Color', 'w');
    ax3 = gca; hold on;
    title('Response Duration Comparison by Condition');
    ylabel('Half-Max Response Time (ms)');

    for i = 1:numConditions
        conditionData = pooledData(strcmp({pooledData.condition}, conditions{i}));
        % Filter for positive, finite values only for time metric
        time_data = [conditionData.halfRmax_time_ms];
        time_data = time_data(time_data > 0 & isfinite(time_data));
        color = lines(numConditions);
        
        HalfViolinPlot(time_data, i, color(i,:), 0.4, 'auto');
    end
    set(ax3, 'xtick', 1:numConditions, 'xticklabel', conditions);

    % --- Save Figures ---
    pipelineSaveFig(fullfile(figFolder, '12_response_magnitude_comparison_by_condition'), 'png', false, fig1);
    pipelineSaveFig(fullfile(figFolder, '13_response_duration_comparison_by_condition'), 'png', false, fig2);
    
    close(fig1);
    close(fig2);

end
