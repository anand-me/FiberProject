%% Fiber Streaming Analysis for Sedimentation Data
% This script analyzes the streaming behavior of flexible fibers during sedimentation
% by examining velocity patterns among individual fibers and clusters.
%
% The script produces:
%   1. Mean velocity vs time plots for different aspect ratios and flexibilities
%   2. Cluster velocity analysis showing how clusters stream differently
%   3. Comparison of isolated fibers vs clustered fibers streaming behavior
%
% Uses output from the cluster analysis code to examine how flexibility and
% aspect ratio affect streaming dynamics in sedimentation.

%% Configuration Parameters
clear; close all;

% Store original directory to return to later
original_dir = pwd;

% Base data directory path - update this with your actual base path
base_data_path = 'E:/FibrePaper';

% Domain size
Lx = 4*pi;
Lz = 8*pi;

% Simulation parameters
N = 359+11; % # of nodes per element depending on aspect ratio
N_proc = 32; % Number of processors used in simulation
timestep = 0.00001; % Simulation timestep
N_save = 5000; % Steps between data dumps

% Time steps to analyze
N_start = 2081;
N_end = 2095;

% Aspect ratios to study (can be expanded)
% aspect_ratios = [2, 10, 15];
% cauchy_numbers = [0.10, 0.42, 0.84, 1.40, 2.10] * 1e-4; % Flexibility values

aspect_ratios = 15;
cauchy_numbers = 0.10 * 1e-4;

% Output settings
outputFolder = 'Streaming_Analysis_Results';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Initialize Arrays for Velocity Data
time_steps = N_start:N_end;
n_steps = length(time_steps);

% Arrays to store results across all AR and Ca combinations
all_mean_velocities = cell(length(aspect_ratios), length(cauchy_numbers));
all_cluster_velocities = cell(length(aspect_ratios), length(cauchy_numbers));
all_isolated_velocities = cell(length(aspect_ratios), length(cauchy_numbers));

%% Loop Through Aspect Ratios and Flexibility Values
for ar_idx = 1:length(aspect_ratios)
    AR = aspect_ratios(ar_idx);
    
    for ca_idx = 1:length(cauchy_numbers)
        Ca = cauchy_numbers(ca_idx);
        
        % Initialize arrays for this specific combination
        mean_Vz = zeros(n_steps, 1);
        cluster_vel_history = cell(n_steps, 1);
        isolated_Vz = cell(n_steps, 1);
        clustered_Vz = cell(n_steps, 1);
        
        % Data folder path based on aspect ratio
        dataFolder = fullfile(base_data_path, sprintf('AR%d/n3/out_v3/', AR));
        
        % Ensure the folder exists
        if ~exist(dataFolder, 'dir')
            warning('Data folder does not exist: %s. Skipping this AR-Ca combination.', dataFolder);
            continue;
        end
        
        % Change to data folder
        cd(dataFolder);
        
        % Format for data files
        fmt = 'tecp.dat%03i.%05i';
        
        % Loop through timesteps
        for t_idx = 1:n_steps
            istep = time_steps(t_idx);
            
            % Initialize arrays for this timestep
            X = []; Z = []; Vx = []; Vz = [];
            Xc = []; Zc = []; Vxc = []; Vzc = [];
            
            % Load data from all processors
            Ng_fiber = 0;
            for i_proc = 0:N_proc-1
                tmp = sprintf(fmt, i_proc, istep);
                filename = fullfile('out_v3', tmp);
                
                % Skip if file doesn't exist
                if ~exist(filename, 'file')
                    fprintf('File not found: %s\n', filename);
                    continue;
                end
                
                % Read data
                fid_r = fopen(filename, 'r');
                if fid_r < 0
                    fprintf('Cannot open file: %s\n', filename);
                    continue;
                end
                
                % Read header
                for i = 1:3
                    str = fgets(fid_r);
                end
                
                % Parse node and element count
                tmp = sscanf(str, 'ZONE N= %i, E= %i,F=FEPOINT, ET=QUADRILATERAL');
                N_tot = tmp(1);
                E_tot = tmp(2);
                
                % Read data
                data = csvread(filename, 3);
                fclose(fid_r);
                
                % Check data consistency
                if (rem(N_tot, N) ~= 0)
                    warning('ERROR: total number of nodes is not divisible by nodes per element!');
                    continue;
                end
                
                % Process fiber data
                N_fiber = N_tot/N; % fibers in local processor
                for i = 1:N_fiber
                    s1 = 1 + N*(i-1);
                    s2 = N*i;
                    i_fiber = i + Ng_fiber;
                    
                    % Store fiber coordinates and velocities
                    X(i_fiber,:) = data(s1:s2, 1);
                    Z(i_fiber,:) = data(s1:s2, 2);
                    Vx(i_fiber,:) = data(s1:s2, 3);
                    Vz(i_fiber,:) = data(s1:s2, 4);
                end
                
                Ng_fiber = Ng_fiber + N_fiber;
            end
            
            % Calculate center positions and velocities
            for i = 1:Ng_fiber
                % Calculate center position
                Xc(i) = mean(X(i,:));
                Zc(i) = mean(Z(i,:));
                
                % Adjust for periodic boundary conditions
                if (Xc(i) < 0)
                    X(i,:) = X(i,:) + Lx;
                    Xc(i) = mean(X(i,:));
                elseif (Xc(i) > Lx)
                    X(i,:) = X(i,:) - Lx;
                    Xc(i) = mean(X(i,:));
                end
                
                if (Zc(i) < 0)
                    Z(i,:) = Z(i,:) + Lz;
                    Zc(i) = mean(Z(i,:));
                elseif (Zc(i) > Lz)
                    Z(i,:) = Z(i,:) - Lz;
                    Zc(i) = mean(Z(i,:));
                end
                
                % Calculate center velocities
                Vxc(i) = mean(Vx(i,:));
                Vzc(i) = mean(Vz(i,:));
            end
            
            % Store mean vertical velocity for this timestep
            mean_Vz(t_idx) = mean(Vzc);
            
            %% Perform clustering analysis to identify fiber clusters
            % Calculate pairwise distances with periodic boundaries
            dist = zeros(Ng_fiber, Ng_fiber);
            for i = 1:Ng_fiber
                for j = 1:Ng_fiber
                    tmpx = abs(Xc(i) - Xc(j));
                    if (tmpx > Lx/2)
                        tmpx = Lx - tmpx;
                    end
                    
                    tmpz = abs(Zc(i) - Zc(j));
                    if (tmpz > Lz/2)
                        tmpz = Lz - tmpz;
                    end
                    
                    dist(i,j) = sqrt(tmpx^2 + tmpz^2);
                end
            end
            
            % Compute similarity matrix
            S = exp(-dist.^2);
            
            % Apply threshold to disconnect distant fibers
            S_eps = S;
            S_eps(S_eps < 0.7) = 0;
            
            % Create graph and perform spectral clustering
            G_eps = graph(S_eps);
            
            % Check if the graph has any edges before proceeding
            if numedges(G_eps) > 0
                k = max(unique(conncomp(G_eps)));
                % Make sure k is at least 1 to avoid empty clusters
                k = max(k, 1);
                idx3 = spectralcluster(S_eps, k, 'Distance', 'precomputed');
            else
                % If no connections, assign all fibers to the same cluster
                fprintf('Warning: No connections found at timestep %d. All fibers considered isolated.\n', istep);
                idx3 = -ones(Ng_fiber, 1);
                k = 0;
            end
            
            % Filter small clusters
            if k > 0  % Only process if clusters were found
                B = unique(idx3);
                Ncount = histc(idx3, B);
                
                for i = 1:Ng_fiber
                    if idx3(i) ~= -1  % Only check valid cluster assignments
                        tmp = idx3(i);
                        if (Ncount(tmp) < 11)  % Ignore clusters with < 11 fibers
                            idx3(i) = -1;
                        end
                    end
                end
            else
                B = -1;  % If no clusters, only have the "isolated" group
            end
            
            % Get updated cluster IDs after filtering
            B = unique(idx3);
            
            % Calculate velocities for each cluster
            cluster_vel = zeros(size(B, 1), 1);
            cluster_count = zeros(size(B, 1), 1);
            
            for i = 1:size(B, 1)
                for j = 1:Ng_fiber
                    if (idx3(j) == B(i))
                        cluster_vel(i) = cluster_vel(i) + Vzc(j);
                        cluster_count(i) = cluster_count(i) + 1;
                    end
                end
            end
            
            % Average velocity for each cluster
            cluster_vel = cluster_vel ./ cluster_count;
            
            % Store cluster velocities for this timestep
            cluster_vel_history{t_idx} = cluster_vel;
            
            % Separate velocities for isolated vs clustered fibers
            isolated = (idx3 == -1);
            isolated_Vz{t_idx} = Vzc(isolated);
            clustered_Vz{t_idx} = Vzc(~isolated);
            
            fprintf('Processed AR=%d, Ca=%g, Timestep=%d: %d fibers, %d clusters\n', ...
                AR, Ca, istep, Ng_fiber, size(B, 1)-1);
        end
        
        % Store results for this AR and Ca combination
        all_mean_velocities{ar_idx, ca_idx} = mean_Vz;
        all_cluster_velocities{ar_idx, ca_idx} = cluster_vel_history;
        all_isolated_velocities{ar_idx, ca_idx} = isolated_Vz;
        
        % Return to original directory
        cd(original_dir);
        
        % Create output directory if it doesn't exist
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end
        
        % Save data for this combination
        save(fullfile(outputFolder, sprintf('velocity_data_AR%d_Ca%g.mat', AR, Ca*1e4)), ...
            'mean_Vz', 'cluster_vel_history', 'isolated_Vz', 'clustered_Vz', 'time_steps');
    end
end

%% Generate Plots

% Plot 1: Mean Velocity vs Time for Different Aspect Ratios and Flexibilities
figure('Position', [100, 100, 1000, 600]);
colors = lines(length(aspect_ratios));
markers = {'o', 's', 'd', '^', 'v'};

for ca_idx = 1:length(cauchy_numbers)
    subplot(2, 3, ca_idx);
    hold on;
    
    for ar_idx = 1:length(aspect_ratios)
        AR = aspect_ratios(ar_idx);
        mean_Vz = all_mean_velocities{ar_idx, ca_idx};
        
        % Skip if data is empty
        if isempty(mean_Vz)
            continue;
        end
        
        plot(time_steps, mean_Vz, 'Color', colors(ar_idx,:), 'Marker', markers{ar_idx}, ...
            'LineWidth', 1.5, 'DisplayName', sprintf('AR=%d', AR));
    end
    
    title(sprintf('Ca = %g × 10^{-4}', cauchy_numbers(ca_idx)*1e4));
    xlabel('Time Step');
    ylabel('Mean Vertical Velocity');
    grid on;
    legend('Location', 'Best');
end

sgtitle('Mean Sedimentation Velocity vs Time for Different Aspect Ratios and Flexibilities');
saveas(gcf, fullfile(outputFolder, 'mean_velocity_comparison.png'));
saveas(gcf, fullfile(outputFolder, 'mean_velocity_comparison.fig'));

% Plot 2: Cluster Velocity Distribution for Different Flexibilities
figure('Position', [100, 100, 1000, 600]);

for ar_idx = 1:length(aspect_ratios)
    AR = aspect_ratios(ar_idx);
    
    subplot(1, length(aspect_ratios), ar_idx);
    hold on;
    
    all_cluster_vels = cell(length(cauchy_numbers), 1);
    
    for ca_idx = 1:length(cauchy_numbers)
        Ca = cauchy_numbers(ca_idx);
        cluster_vels = [];
        
        cluster_vel_history = all_cluster_velocities{ar_idx, ca_idx};
        if isempty(cluster_vel_history)
            continue;
        end
        
        % Combine all cluster velocities across timesteps
        for t_idx = 1:length(cluster_vel_history)
            if ~isempty(cluster_vel_history{t_idx})
                % Skip the first entry (idx = -1, isolated fibers)
                cluster_vels = [cluster_vels; cluster_vel_history{t_idx}(2:end)];
            end
        end
        
        all_cluster_vels{ca_idx} = cluster_vels;
    end
    
    % Create boxplot for this aspect ratio
    boxplot(cell2mat(all_cluster_vels), repelem(1:length(cauchy_numbers), cellfun(@length, all_cluster_vels)), ...
        'Labels', arrayfun(@(x) sprintf('%g', x*1e4), cauchy_numbers, 'UniformOutput', false), ...
        'Colors', lines(length(cauchy_numbers)));
    
    title(sprintf('AR = %d', AR));
    xlabel('Cauchy Number (× 10^{-4})');
    ylabel('Cluster Vertical Velocity');
    grid on;
end

sgtitle('Distribution of Cluster Velocities for Different Flexibilities');
saveas(gcf, fullfile(outputFolder, 'cluster_velocity_distribution.png'));
saveas(gcf, fullfile(outputFolder, 'cluster_velocity_distribution.fig'));

% Plot 3: Comparison of Isolated vs Clustered Fibers
figure('Position', [100, 100, 1000, 600]);

for ar_idx = 1:length(aspect_ratios)
    AR = aspect_ratios(ar_idx);
    
    subplot(1, length(aspect_ratios), ar_idx);
    hold on;
    
    isolated_means = zeros(length(cauchy_numbers), 1);
    clustered_means = zeros(length(cauchy_numbers), 1);
    
    for ca_idx = 1:length(cauchy_numbers)
        Ca = cauchy_numbers(ca_idx);
        
        isolated_Vz = all_isolated_velocities{ar_idx, ca_idx};
        if isempty(isolated_Vz)
            continue;
        end
        
        % Calculate mean velocities for isolated and clustered fibers
        isolated_all = [];
        clustered_all = [];
        
        for t_idx = 1:length(isolated_Vz)
            if ~isempty(isolated_Vz{t_idx})
                isolated_all = [isolated_all; isolated_Vz{t_idx}];
                clustered_all = [clustered_all; all_cluster_velocities{ar_idx, ca_idx}{t_idx}];
            end
        end
        
        isolated_means(ca_idx) = mean(isolated_all);
        clustered_means(ca_idx) = mean(clustered_all);
    end
    
    % Plot comparison
    bar([isolated_means, clustered_means]);
    
    set(gca, 'XTick', 1:length(cauchy_numbers), ...
        'XTickLabel', arrayfun(@(x) sprintf('%g', x*1e4), cauchy_numbers, 'UniformOutput', false));
    
    title(sprintf('AR = %d', AR));
    xlabel('Cauchy Number (× 10^{-4})');
    ylabel('Mean Vertical Velocity');
    legend('Isolated Fibers', 'Clustered Fibers');
    grid on;
end

sgtitle('Comparison of Isolated vs Clustered Fiber Velocities');
saveas(gcf, fullfile(outputFolder, 'isolated_vs_clustered_comparison.png'));
saveas(gcf, fullfile(outputFolder, 'isolated_vs_clustered_comparison.fig'));

%% Plot 4: Streaming Enhancement Factor vs Flexibility
figure('Position', [100, 100, 800, 600]);
hold on;

for ar_idx = 1:length(aspect_ratios)
    AR = aspect_ratios(ar_idx);
    
    enhancement_factor = zeros(length(cauchy_numbers), 1);
    
    for ca_idx = 1:length(cauchy_numbers)
        isolated_Vz = all_isolated_velocities{ar_idx, ca_idx};
        if isempty(isolated_Vz)
            continue;
        end
        
        % Calculate mean velocities for isolated and clustered fibers
        isolated_all = [];
        clustered_all = [];
        
        for t_idx = 1:length(isolated_Vz)
            if ~isempty(isolated_Vz{t_idx})
                isolated_all = [isolated_all; isolated_Vz{t_idx}];
                clustered_all = [clustered_all; all_cluster_velocities{ar_idx, ca_idx}{t_idx}];
            end
        end
        
        % Enhancement factor = ratio of clustered to isolated velocity
        if ~isempty(isolated_all) && mean(isolated_all) ~= 0
            enhancement_factor(ca_idx) = mean(clustered_all) / mean(isolated_all);
        else
            enhancement_factor(ca_idx) = NaN;
        end
    end
    
    % Plot enhancement factor vs flexibility
    plot(cauchy_numbers*1e4, enhancement_factor, 'o-', 'LineWidth', 2, ...
        'DisplayName', sprintf('AR=%d', AR), 'Color', colors(ar_idx,:));
end

title('Streaming Enhancement Factor vs Flexibility');
xlabel('Cauchy Number (× 10^{-4})');
ylabel('Streaming Enhancement Factor (V_{clustered}/V_{isolated})');
legend('Location', 'Best');
grid on;

saveas(gcf, fullfile(outputFolder, 'streaming_enhancement_factor.png'));
saveas(gcf, fullfile(outputFolder, 'streaming_enhancement_factor.fig'));

%% Save Final Data for Paper Discussion
save(fullfile(outputFolder, 'streaming_analysis_complete.mat'), ...
    'aspect_ratios', 'cauchy_numbers', 'all_mean_velocities', ...
    'all_cluster_velocities', 'all_isolated_velocities', 'time_steps');

fprintf('Analysis complete. Results saved in %s\n', outputFolder);

% Return to original directory if not already there
if ~strcmp(pwd, original_dir)
    cd(original_dir);
end