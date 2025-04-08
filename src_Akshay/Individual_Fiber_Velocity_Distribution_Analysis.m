%% Individual Fiber Velocity Distribution Analysis
% This script analyzes the distribution of individual fiber velocities (Vzc)
% from the previously saved data files.

clear; close all;

% Create output directory for plots if it doesn't exist
if ~exist('./plots/fiber_distributions', 'dir')
    mkdir('./plots/fiber_distributions');
end

% Define aspect ratios to analyze
aspect_ratios = [2, 10, 15, 20];

% Load all AR data (change the filename to match your saved file)
load('mean_velocity_all_AR_range_id_2.mat');

%% 1. Statistical Analysis of Fiber Velocities
% Initialize structure to store statistics
fiber_stats = struct();

% Process each aspect ratio
for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    ar_field = sprintf('AR%d', ar);
    
    % Check if data exists for this AR
    if ~isfield(results, ar_field)
        fprintf('No data found for AR = %d. Skipping...\n', ar);
        continue;
    end
    
    % Get data
    ar_data = results.(ar_field);
    all_vzc = ar_data.all_Vzc;
    
    % Calculate time-averaged statistics
    n_timesteps = length(all_vzc);
    
    % Initialize arrays to store statistics for each timestep
    mean_values = zeros(n_timesteps, 1);
    median_values = zeros(n_timesteps, 1);
    std_values = zeros(n_timesteps, 1);
    min_values = zeros(n_timesteps, 1);
    max_values = zeros(n_timesteps, 1);
    iqr_values = zeros(n_timesteps, 1);  % Interquartile range
    
    % Analyze each timestep
    for t = 1:n_timesteps
        vzc_t = all_vzc{t};
        
        % Skip if empty
        if isempty(vzc_t)
            continue;
        end
        
        % Calculate statistics
        mean_values(t) = mean(vzc_t, 'omitnan');
        median_values(t) = median(vzc_t, 'omitnan');
        std_values(t) = std(vzc_t, 'omitnan');
        min_values(t) = min(vzc_t);
        max_values(t) = max(vzc_t);
        q = quantile(vzc_t, [0.25, 0.75]);
        iqr_values(t) = q(2) - q(1);
    end
    
    % Store statistics
    fiber_stats.(ar_field) = struct(...
        'mean', mean_values, ...
        'median', median_values, ...
        'std', std_values, ...
        'min', min_values, ...
        'max', max_values, ...
        'iqr', iqr_values, ...
        'time_steps', ar_data.time_steps);
    
    % Print average statistics
    fprintf('AR = %d Statistics:\n', ar);
    fprintf('  Average Mean Velocity: %.4f\n', mean(mean_values, 'omitnan'));
    fprintf('  Average Median Velocity: %.4f\n', mean(median_values, 'omitnan'));
    fprintf('  Average Standard Deviation: %.4f\n', mean(std_values, 'omitnan'));
    fprintf('  Average IQR: %.4f\n', mean(iqr_values, 'omitnan'));
    fprintf('  Min Velocity: %.4f\n', min(min_values));
    fprintf('  Max Velocity: %.4f\n', max(max_values));
    fprintf('\n');
end

% Save statistics
save('./plots/fiber_distributions/fiber_velocity_statistics.mat', 'fiber_stats');

%% 2. Visualization: Box Plots of Velocity Distribution Over Time
for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    ar_field = sprintf('AR%d', ar);
    
    % Check if data exists for this AR
    if ~isfield(results, ar_field)
        continue;
    end
    
    % Get data
    ar_data = results.(ar_field);
    all_vzc = ar_data.all_Vzc;
    time_steps = ar_data.time_steps;
    
    % Select sample points (e.g., every 10th timestep)
    % Adjust this based on your data size
    sample_interval = max(1, floor(length(all_vzc)/10));
    sample_indices = 1:sample_interval:length(all_vzc);
    
    % Create figure for box plots
    figure('Position', [100, 100, 1200, 600]);
    
    % Prepare data for boxplot
    boxplot_data = [];
    group_labels = [];
    
    for i = 1:length(sample_indices)
        t_idx = sample_indices(i);
        
        % Skip if empty
        if t_idx > length(all_vzc) || isempty(all_vzc{t_idx})
            continue;
        end
        
        vzc_t = all_vzc{t_idx};
        boxplot_data = [boxplot_data; vzc_t];
        group_labels = [group_labels; repmat(t_idx, length(vzc_t), 1)];
    end
    
    % Create boxplot
    boxplot(boxplot_data, group_labels, 'Labels', time_steps(sample_indices));
    title(sprintf('Distribution of Fiber Velocities Over Time (AR = %d)', ar), 'FontSize', 14);
    xlabel('Time Step', 'FontSize', 12);
    ylabel('Vertical Velocity (V_z)', 'FontSize', 12);
    grid on;
    
    % Save figure
    saveas(gcf, sprintf('./plots/fiber_distributions/boxplot_AR%d.png', ar));
    saveas(gcf, sprintf('./plots/fiber_distributions/boxplot_AR%d.fig', ar));
end

%% 3. Visualization: Alternative Distribution Plots
% Create separate histogram plots for each timestep instead of violin plots
for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    ar_field = sprintf('AR%d', ar);
    
    % Check if data exists for this AR
    if ~isfield(results, ar_field)
        continue;
    end
    
    % Get data
    ar_data = results.(ar_field);
    all_vzc = ar_data.all_Vzc;
    time_steps = ar_data.time_steps;
    
    % Select 4-5 representative timesteps
    % Adjust this based on your data
    n_samples = min(5, length(all_vzc));
    sample_indices = round(linspace(1, length(all_vzc), n_samples));
    
    % Create figure with subplots for each timestep
    figure('Position', [100, 100, 1200, 800]);
    
    for i = 1:n_samples
        t_idx = sample_indices(i);
        
        % Skip if index out of range or data empty
        if t_idx > length(all_vzc) || isempty(all_vzc{t_idx})
            continue;
        end
        
        % Get data for this timestep
        vzc_t = all_vzc{t_idx};
        
        % Create subplot
        subplot(2, ceil(n_samples/2), i);
        
        % Create histogram
        histogram(vzc_t, 20, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.7);
        
        % Add kernel density estimate if available
        if exist('ksdensity', 'file') == 2
            hold on;
            [f, xi] = ksdensity(vzc_t);
            plot(xi, f, 'r-', 'LineWidth', 2);
        end
        
        % Add label
        if t_idx <= length(time_steps)
            title(sprintf('t = %d', time_steps(t_idx)), 'FontSize', 12);
        end
        
        xlabel('Vertical Velocity (V_z)', 'FontSize', 10);
        ylabel('Frequency', 'FontSize', 10);
        grid on;
    end
    
    sgtitle(sprintf('Distribution of Fiber Velocities Over Time (AR = %d)', ar), 'FontSize', 14);
    
    % Save figure
    saveas(gcf, sprintf('./plots/fiber_distributions/distribution_over_time_AR%d.png', ar));
    saveas(gcf, sprintf('./plots/fiber_distributions/distribution_over_time_AR%d.fig', ar));
end

%% 4. Visualization: Histograms of Fiber Velocities
% Create a figure with histograms for each AR at a specific timestep
% Choose a representative timestep (e.g., halfway through the simulation)
timestep_idx = cell(length(aspect_ratios), 1);

for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    ar_field = sprintf('AR%d', ar);
    
    % Check if data exists for this AR
    if ~isfield(results, ar_field)
        continue;
    end
    
    % Get data
    ar_data = results.(ar_field);
    all_vzc = ar_data.all_Vzc;
    
    % Select middle timestep
    timestep_idx{ar_idx} = floor(length(all_vzc)/2);
end

% Create a multi-panel histogram figure
figure('Position', [100, 100, 1200, 800]);

for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    ar_field = sprintf('AR%d', ar);
    
    % Skip if data is missing
    if ~isfield(results, ar_field) || isempty(timestep_idx{ar_idx})
        continue;
    end
    
    % Get data for the selected timestep
    ar_data = results.(ar_field);
    vzc = ar_data.all_Vzc{timestep_idx{ar_idx}};
    
    % Skip if data is empty
    if isempty(vzc)
        continue;
    end
    
    % Create subplot
    subplot(2, 2, ar_idx);
    
    % Create histogram with kernel density estimate
    histogram(vzc, 20, 'Normalization', 'probability', 'FaceColor', 'auto', 'FaceAlpha', 0.7);
    hold on;
    
    % Add kernel density estimate if Statistics Toolbox is available
    if exist('ksdensity', 'file') == 2
        [f, xi] = ksdensity(vzc);
        plot(xi, f, 'r-', 'LineWidth', 2);
    end
    
    title(sprintf('AR = %d', ar), 'FontSize', 12);
    xlabel('Vertical Velocity (V_z)', 'FontSize', 10);
    ylabel('Frequency', 'FontSize', 10);
    grid on;
end

sgtitle('Distribution of Individual Fiber Velocities', 'FontSize', 14);

% Save figure
saveas(gcf, './plots/fiber_distributions/histograms_all_AR.png');
saveas(gcf, './plots/fiber_distributions/histograms_all_AR.fig');

%% 5. Visualization: Velocity Standard Deviation Over Time
% This shows how the spread of fiber velocities changes over time
figure('Position', [100, 100, 1000, 600]);

% Colors for different ARs
colors = {'b', 'k', 'r', 'm'};

hold on;
for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    ar_field = sprintf('AR%d', ar);
    
    % Skip if data is missing
    if ~isfield(fiber_stats, ar_field)
        continue;
    end
    
    % Get statistics
    stats = fiber_stats.(ar_field);
    
    % Plot standard deviation over time
    time_steps = stats.time_steps - stats.time_steps(1); % Normalize to start at 0
    plot(time_steps, stats.std, colors{ar_idx}, 'LineWidth', 1.5);
end

title('Standard Deviation of Fiber Velocities Over Time', 'FontSize', 14);
xlabel('Non-dimensional Time', 'FontSize', 12);
ylabel('Standard Deviation of V_z', 'FontSize', 12);
legend(arrayfun(@(x) sprintf('AR = %d', x), aspect_ratios, 'UniformOutput', false), 'Location', 'best');
grid on;

% Save figure
saveas(gcf, './plots/fiber_distributions/velocity_std_over_time.png');
saveas(gcf, './plots/fiber_distributions/velocity_std_over_time.fig');

%% 6. Visualization: Min-Max Range Plot
% This shows the range of fiber velocities over time
figure('Position', [100, 100, 1000, 600]);

for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    ar_field = sprintf('AR%d', ar);
    
    % Skip if data is missing
    if ~isfield(fiber_stats, ar_field)
        continue;
    end
    
    % Get statistics
    stats = fiber_stats.(ar_field);
    
    % Create subplot
    subplot(2, 2, ar_idx);
    
    % Normalize time to start at 0
    time_steps = stats.time_steps - stats.time_steps(1);
    
    % Plot mean as a line
    plot(time_steps, stats.mean, 'k-', 'LineWidth', 1.5);
    hold on;
    
    % Plot min and max as separate lines instead of using fill
    plot(time_steps, stats.min, '--', 'Color', colors{ar_idx}, 'LineWidth', 1);
    plot(time_steps, stats.max, '--', 'Color', colors{ar_idx}, 'LineWidth', 1);
    
    % Add shading between min and max manually if vectors are same length
    if all(~isnan(stats.min)) && all(~isnan(stats.max)) && length(stats.min) == length(stats.max)
        x_fill = [time_steps; flip(time_steps)];
        y_fill = [stats.min; flip(stats.max)];
        
        % Check if dimensions are compatible
        if size(x_fill, 1) == size(y_fill, 1) && size(x_fill, 2) == size(y_fill, 2)
            fill(x_fill, y_fill, colors{ar_idx}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
    
    title(sprintf('Velocity Range (AR = %d)', ar), 'FontSize', 12);
    xlabel('Non-dimensional Time', 'FontSize', 10);
    ylabel('Vertical Velocity (V_z)', 'FontSize', 10);
    grid on;
end

sgtitle('Range of Fiber Velocities Over Time', 'FontSize', 14);

% Save figure
saveas(gcf, './plots/fiber_distributions/velocity_range_over_time.png');
saveas(gcf, './plots/fiber_distributions/velocity_range_over_time.fig');

%% 7. Clustering Analysis: Find Velocity-Based Clusters
% This section attempts to detect clusters of fibers with similar velocities
% using k-means clustering

% Try this analysis for one timestep for each AR
for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    ar_field = sprintf('AR%d', ar);
    
    % Skip if data is missing
    if ~isfield(results, ar_field)
        continue;
    end
    
    % Get data
    ar_data = results.(ar_field);
    all_vzc = ar_data.all_Vzc;
    
    % Select a timestep with good data (middle of simulation)
    t_idx = floor(length(all_vzc)/2);
    if t_idx < 1 || t_idx > length(all_vzc) || isempty(all_vzc{t_idx})
        continue;
    end
    
    vzc = all_vzc{t_idx};
    
    % Skip if not enough data points
    if length(vzc) < 10
        continue;
    end
    
    % Determine optimal number of clusters using silhouette method
    max_clusters = min(5, floor(length(vzc)/10)); % Don't try too many clusters
    silhouette_vals = zeros(max_clusters-1, 1);
    
    for k = 2:max_clusters
        % Perform k-means clustering
        [idx, ~] = kmeans(vzc, k);
        
        % Calculate silhouette value
        if exist('silhouette', 'file') == 2
            s = silhouette(vzc, idx);
            silhouette_vals(k-1) = mean(s);
        else
            % If Statistics Toolbox is not available, just use k=3
            silhouette_vals = [0 1 0];
            break;
        end
    end
    
    % Find optimal number of clusters
    [~, best_k_idx] = max(silhouette_vals);
    best_k = best_k_idx + 1; % Add 1 because we started at k=2
    
    % Perform final clustering with optimal k
    [idx, centroids] = kmeans(vzc, best_k);
    
    % Create figure to visualize clusters
    figure('Position', [100, 100, 800, 600]);
    
    % Create histogram with different colors for each cluster
    hold on;
    cluster_colors = hsv(best_k); % Generate distinct colors
    
    for k = 1:best_k
        % Get velocities for this cluster
        cluster_vzc = vzc(idx == k);
        
        % Create histogram for this cluster
        histogram(cluster_vzc, 15, 'FaceColor', cluster_colors(k,:), 'FaceAlpha', 0.6, ...
            'EdgeColor', 'none', 'DisplayName', sprintf('Cluster %d (n=%d)', k, sum(idx == k)));
        
        % Plot centroid
        line([centroids(k) centroids(k)], ylim, 'Color', cluster_colors(k,:), 'LineWidth', 2, 'LineStyle', '--');
    end
    
    title(sprintf('Velocity Clusters for AR = %d (k = %d)', ar, best_k), 'FontSize', 14);
    xlabel('Vertical Velocity (V_z)', 'FontSize', 12);
    ylabel('Count', 'FontSize', 12);
    legend('Location', 'best');
    grid on;
    
    % Save figure
    saveas(gcf, sprintf('./plots/fiber_distributions/velocity_clusters_AR%d.png', ar));
    saveas(gcf, sprintf('./plots/fiber_distributions/velocity_clusters_AR%d.fig', ar));
    
    % Report cluster statistics
    fprintf('AR = %d Velocity Clustering (k = %d):\n', ar, best_k);
    for k = 1:best_k
        cluster_vzc = vzc(idx == k);
        fprintf('  Cluster %d: %d fibers (%.1f%%), Mean Vz = %.4f\n', ...
            k, sum(idx == k), 100*sum(idx == k)/length(vzc), mean(cluster_vzc));
    end
    fprintf('\n');
end

fprintf('Analysis complete. All plots saved in ./plots/fiber_distributions/\n');