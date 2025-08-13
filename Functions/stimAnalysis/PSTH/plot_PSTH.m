%% SCRIPT TO COMPARE GROUPS OF NETWORK RESPONSES AND DECAY 

clear; clc; close all;

%% --- USER PARAMETERS ---

% --- Define the parent folders for the two groups ---
% This script provisionally requires you to place all output folders into a
% parent folder called HIGH AC and LOW AC
folder_group1 = 'HIGH AC';
folder_group2 = 'LOW AC';

% --- Define short labels for the x-axis of the plots ---
label1 = 'High AC';
label2 = 'Low AC';

% --- Analysis Option: Filter out non-positive values ---
% Set to 'true' to only include values > 0 in calculations and plots.
% Set to 'false' to include all finite values (including negative ones).
filter_positive_only = false;

%% --- HELPER FUNCTIONS ---

% Function to load and pool all electrode data from a group folder
function electrode_data = load_all_electrode_data(parent_folder)
    all_files = dir(parent_folder);
    sub_folders = all_files([all_files.isdir]);
    sub_folders = sub_folders(~ismember({sub_folders.name},{'.','..'}));

    electrode_data = struct('auc_response', [], 'auc_corrected', [], 'halfRmax_time_ms', []);
    for i = 1:length(sub_folders)
        folder_name = fullfile(parent_folder, sub_folders(i).name);
        matFile_info = dir(fullfile(folder_name, 'networkResponse*.mat'));
        if isempty(matFile_info), continue; end

        file_path = fullfile(folder_name, matFile_info(1).name);
        loaded_data = load(file_path, 'networkResponse');
        d = loaded_data.networkResponse;

        electrode_data.auc_response = [electrode_data.auc_response, [d.auc_response]];
        electrode_data.auc_corrected = [electrode_data.auc_corrected, [d.auc_corrected]];
        electrode_data.halfRmax_time_ms = [electrode_data.halfRmax_time_ms, [d.halfRmax_time_ms]];
    end
end

% Function to load and pair file-level data based on culture ID
function paired_file_data = load_paired_data(folder1, folder2, filter_positive)

    function id_map = create_culture_id_map(parent_folder)
        id_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
        dir_info = dir(parent_folder);
        sub_folders = {dir_info([dir_info.isdir]).name};
        sub_folders = sub_folders(~ismember(sub_folders,{'.','..'}));

        for i = 1:length(sub_folders)
            folder_name = sub_folders{i};
            token = regexp(folder_name, 'OWT220207_(\d+[A-Z])', 'tokens');
            if ~isempty(token) && ~isempty(token{1})
                id_map(token{1}{1}) = fullfile(parent_folder, folder_name);
            end
        end
    end

    map1 = create_culture_id_map(folder1);
    map2 = create_culture_id_map(folder2);

    common_ids = intersect(keys(map1), keys(map2));
    num_pairs = length(common_ids);

    if num_pairs == 0, error('No common culture IDs found to perform paired analysis.'); end
    fprintf('Found %d paired recordings based on culture ID.\n', num_pairs);

    paired_file_data = struct('auc_response', nan(num_pairs, 2), 'auc_corrected', nan(num_pairs, 2), 'halfRmax_time_ms', nan(num_pairs, 2));

    for i = 1:num_pairs
        folder_paths = {map1(common_ids{i}), map2(common_ids{i})};
        for j = 1:2
            matFile_info = dir(fullfile(folder_paths{j}, 'networkResponse*.mat'));
            if isempty(matFile_info), continue; end

            d = load(fullfile(folder_paths{j}, matFile_info(1).name), 'networkResponse').networkResponse;

            % Apply filter conditionally before taking the mean
            filter_vec = @(vec) vec( (filter_positive && vec > 0) | (~filter_positive & isfinite(vec)) );

            paired_file_data.auc_response(i, j) = mean(filter_vec([d.auc_response]));
            paired_file_data.auc_corrected(i, j) = mean(filter_vec([d.auc_corrected]));
            paired_file_data.halfRmax_time_ms(i, j) = mean(filter_vec([d.halfRmax_time_ms]));
        end
    end
end

%% --- LOAD AND PROCESS DATA FOR BOTH LEVELS ---

% Electrode-level (unpaired)
electrode_data1 = load_all_electrode_data(folder_group1);
electrode_data2 = load_all_electrode_data(folder_group2);

% File-level (paired)
paired_file_data = load_paired_data(folder_group1, folder_group2, filter_positive_only);

%% --- PLOTTING SETUP ---

color1 = [0.8500, 0.3250, 0.0980]; % Reddish-orange
color2 = [0.0, 0.4470, 0.7410];   % Blue

% Helper function for violin plots with conditional filtering
function plot_violin_scatter(ax, data, x_pos, group_color, filter_positive)
    if isempty(data), return; end

    % Apply filter based on the master flag
    if filter_positive
        data = data(data > 0 & isfinite(data));
    else
        data = data(isfinite(data));
    end
    if isempty(data), return; end

    plot_width = 0.4;
    [f, yi] = ksdensity(data);
    f = f / max(f) * plot_width;
    fill(ax, [f, zeros(1, numel(f))] + x_pos, [yi, fliplr(yi)], group_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    jitter = (rand(size(data)) - 0.5) * plot_width;
    scatter(ax, x_pos + jitter - (plot_width/2), data, 20, group_color, 'filled', 'MarkerFaceAlpha', 0.6);

    mean_val = mean(data);
    std_val = std(data); % Calculate Standard Deviation

    % Plot the Standard Deviation as a vertical bar
    plot(ax, [x_pos, x_pos], [mean_val - std_val, mean_val + std_val], 'k', 'LineWidth', 2.5);
    scatter(ax, x_pos, mean_val, 80, 'k', 'filled');

    n_count = numel(data);
    text(ax, x_pos, max(ylim(ax)), sprintf('n = %d', n_count), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
end

%% --- FIGURES 1 & 2: ELECTRODE-LEVEL COMPARISON ---

fprintf('\nGenerating electrode-level comparison figures...\n');
fig1 = figure('Position', [100, 100, 1000, 500], 'Color', 'w', 'Name', 'Electrode-Level AUC Comparison');
ax1 = subplot(1, 2, 1); hold on;
plot_violin_scatter(ax1, electrode_data1.auc_response, 1, color1, filter_positive_only);
plot_violin_scatter(ax1, electrode_data2.auc_response, 2, color2, filter_positive_only);
title(ax1, 'Response AUC'); ylabel(ax1, 'AUC');
ax2 = subplot(1, 2, 2); hold on;
plot_violin_scatter(ax2, electrode_data1.auc_corrected, 1, color1, filter_positive_only);
plot_violin_scatter(ax2, electrode_data2.auc_corrected, 2, color2, filter_positive_only);
title(ax2, 'Baseline-Corrected AUC'); ylabel(ax2, 'AUC');
sgtitle(fig1, 'Electrode-Level: Magnitude of Response Comparison, DIV 63, 6uA Stimulation', 'FontSize', 16, 'FontWeight', 'bold');

fig2 = figure('Position', [150, 150, 600, 500], 'Color', 'w', 'Name', 'Electrode-Level Decay Comparison');
ax3 = gca; hold on;
plot_violin_scatter(ax3, electrode_data1.halfRmax_time_ms, 1, color1, filter_positive_only);
plot_violin_scatter(ax3, electrode_data2.halfRmax_time_ms, 2, color2, filter_positive_only);
ylabel(ax3, 'Half-Max Response Time (ms)');
sgtitle(fig2, 'Electrode-Level: Response Decay Comparison, DIV 63, 6uA Stimulation', 'FontSize', 16, 'FontWeight', 'bold');

for ax = [ax1, ax2, ax3]
    set(ax, 'xtick', [1, 2], 'xticklabel', {label1, label2}, 'FontSize', 12, 'box', 'off', 'TickDir', 'out');
    xlim(ax, [0.5, 2.5]); grid(ax, 'on');
    if filter_positive_only
        ylim(ax, [0, max(ylim(ax)) * 1.1]);
    else
        current_ylim = ylim(ax);
        ylim(ax, [current_ylim(1), current_ylim(2) * 1.1]);
    end
end

%% --- FIGURES 3 & 4: FILE-LEVEL COMPARISON ---

fprintf('Generating file-level comparison figures...\n');
fig3 = figure('Position', [200, 200, 1000, 500], 'Color', 'w', 'Name', 'File-Level AUC Comparison');
ax4 = subplot(1, 2, 1); hold on;
plot(ax4, [1, 2], paired_file_data.auc_response', '-', 'Color', [0.7 0.7 0.7]);
plot_violin_scatter(ax4, paired_file_data.auc_response(:,1), 1, color1, filter_positive_only);
plot_violin_scatter(ax4, paired_file_data.auc_response(:,2), 2, color2, filter_positive_only);
title(ax4, 'Mean Response AUC'); ylabel(ax4, 'AUC');

ax5 = subplot(1, 2, 2); hold on;
plot(ax5, [1, 2], paired_file_data.auc_corrected', '-', 'Color', [0.7 0.7 0.7]);
plot_violin_scatter(ax5, paired_file_data.auc_corrected(:,1), 1, color1, filter_positive_only);
plot_violin_scatter(ax5, paired_file_data.auc_corrected(:,2), 2, color2, filter_positive_only);
title(ax5, 'Mean Baseline-Corrected AUC'); ylabel(ax5, 'AUC');
sgtitle(fig3, 'File-Level: Magnitude of Response Comparison, DIV 63, 6uA Stimulation', 'FontSize', 16, 'FontWeight', 'bold');

fig4 = figure('Position', [250, 250, 600, 500], 'Color', 'w', 'Name', 'File-Level Decay Comparison');
ax6 = gca; hold on;
plot(ax6, [1, 2], paired_file_data.halfRmax_time_ms', '-', 'Color', [0.7 0.7 0.7]);
plot_violin_scatter(ax6, paired_file_data.halfRmax_time_ms(:,1), 1, color1, filter_positive_only);
plot_violin_scatter(ax6, paired_file_data.halfRmax_time_ms(:,2), 2, color2, filter_positive_only);
ylabel(ax6, 'Mean Half-Max Response Time (ms) per File');
sgtitle(fig4, 'File-Level: Response Decay Comparison, DIV 63, 6uA Stimulation', 'FontSize', 16, 'FontWeight', 'bold');

for ax = [ax4, ax5, ax6]
    set(ax, 'xtick', [1, 2], 'xticklabel', {label1, label2}, 'FontSize', 12, 'box', 'off', 'TickDir', 'out');
    xlim(ax, [0.5, 2.5]); grid(ax, 'on');
    if filter_positive_only
        ylim(ax, [0, max(ylim(ax)) * 1.1]);
    else
        current_ylim = ylim(ax);
        ylim(ax, [current_ylim(1), current_ylim(2) * 1.1]);
    end
end
linkaxes([ax1, ax2, ax4, ax5], 'y');

fprintf('\nAnalysis complete. Four figures generated.\n');