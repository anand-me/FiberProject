%% Fiber Cluster Analysis for Sedimentation Simulations
% This code analyzes fiber clustering patterns during sedimentation using spectral clustering
% with support for different aspect ratios and time filtering.

clear; close all;

%% Configuration Parameters
% Domain dimensions
Lx = 4*pi;
Lz = 8*pi;

% Simulation parameters
N_proc = 32;
fmt = 'tecp.dat%03i.%05i';

% Time filtering settings
timestep = 0.00001;  % Simulation timestep
N_save = 5000;       % Steps between data dumps
dt = N_save * timestep;

% Define filter widths (in real time units)
filter_width_0 = 0;        % Instantaneous (no filtering)
filter_width_50 = 50*dt;   % Filter width 50
filter_width_100 = 100*dt; % Filter width 100

% Select which filter width to use
selected_filter = filter_width_50; % Change to filter_width_50 or filter_width_100 as needed

% Store original directory (where the script is)
originalDir = pwd;
fprintf('Original directory: %s\n', originalDir);

% Main output directory (in the original directory)
mainFolder = fullfile(originalDir, 'Fiber_Clustering_Results');
if ~exist(mainFolder, 'dir')
    mkdir(mainFolder);
end

%% Aspect Ratio Settings
% Define aspect ratios to analyze
aspect_ratios = [2, 10, 15, 20];

% Define file ranges for each aspect ratio
file_ranges = {
    [22001, 22203],  % AR 2 (range 3)
    [25001, 25204],  % AR 10 (range 3)
    [25001, 25228],  % AR 15 (range 3)
    [501, 600]       % AR 20 (range 2)
};

% Node counts for each aspect ratio
nodes_per_AR = [366, 359, 370, 405];

% Data directories for each aspect ratio
data_folders = {
    'E:/FibrePaper/AR2/n3/out_v3/',
    'E:/FibrePaper/AR10/n3/out_v3/',
    'E:/FibrePaper/AR15/n3/out_v3/',
    'D:/OneDrive - Florida State University/Vahid_Material/FFT_IBM/suspension/medium/DEN_10_g_1000/S1_more_flexible/n3/out_v3/'
};

%% Process Each Aspect Ratio
for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    fprintf('Processing Aspect Ratio %d...\n', ar);
    
    % Get node count for this aspect ratio
    N = nodes_per_AR(ar_idx);
    
    % Get data folder
    dataFolder = data_folders{ar_idx};
    
    % Get file range and calculate mid-point
    range = file_ranges{ar_idx};
    mid_step = round(mean(range));
    
    % Create output folder for this AR
    outputFolder = fullfile(mainFolder, sprintf('AR%d_Filter%.0f', ar, selected_filter));
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    fprintf('  Using timestep %d (midpoint of range)...\n', mid_step);
    
    % Change to data directory for reading
    cd(dataFolder);
    fprintf('  Changed to data directory: %s\n', pwd);
    
    % Load fiber data for this timestep
    X = []; Z = []; Vx = []; Vz = [];
    Ng_fiber = 0;
    
    % Process each processor
    for i_proc = 0:N_proc-1
        tmp = sprintf(fmt, i_proc, mid_step);
        
        % Check if file exists
        if ~exist(tmp, 'file')
            continue;
        end
        
        % Open file
        fid_r = fopen(tmp, 'r');
        if fid_r < 0
            continue;
        end
        
        % Read header
        for i = 1:3
            str = fgets(fid_r);
        end
        
        % Parse header
        tmp_header = sscanf(str, 'ZONE N= %i, E= %i,F=FEPOINT, ET=QUADRILATERAL');
        N_tot = tmp_header(1);
        E_tot = tmp_header(2);
        
        % Read data
        tmp_data = csvread(tmp, 3);
        fclose(fid_r);
        
        % Extract node data
        data = tmp_data(1:N_tot, :);
        
        % Check if total nodes divisible by nodes per element
        if (rem(N_tot, N) ~= 0)
            fprintf('ERROR: total number of nodes is not divisible by node per element!\n');
            continue;
        end
        
        % Process each fiber
        N_fiber = N_tot/N;
        for i = 1:N_fiber
            s1 = 1 + N*(i-1);
            s2 = N*i;
            i_fiber = i + Ng_fiber;
            
            % Allocate arrays if first processor
            if i_proc == 0 && i == 1
                X = zeros(1000, N);
                Z = zeros(1000, N);
                Vx = zeros(1000, N);
                Vz = zeros(1000, N);
            end
            
            % Expand arrays if needed
            if i_fiber > size(X, 1)
                X = [X; zeros(500, N)];
                Z = [Z; zeros(500, N)];
                Vx = [Vx; zeros(500, N)];
                Vz = [Vz; zeros(500, N)];
            end
            
            % Store data
            X(i_fiber, 1:s2-s1+1) = data(s1:s2, 1);
            Z(i_fiber, 1:s2-s1+1) = data(s1:s2, 2);
            Vx(i_fiber, 1:s2-s1+1) = data(s1:s2, 3);
            Vz(i_fiber, 1:s2-s1+1) = data(s1:s2, 4);
        end
        is there osmthi g else that I need to chnage 
        Ng_fiber = i_fiber;
    end
    
    % Return to original directory for saving results
    cd(originalDir);
    fprintf('  Returned to original directory: %s\n', pwd);
    
    % Trim arrays to actual number of fibers
    X = X(1:Ng_fiber, :);
    Z = Z(1:Ng_fiber, :);
    Vx = Vx(1:Ng_fiber, :);
    Vz = Vz(1:Ng_fiber, :);
    
    fprintf('  Loaded data for %d fibers\n', Ng_fiber);
    
    % Calculate fiber centers
    Xc = zeros(Ng_fiber, 1);
    Zc = zeros(Ng_fiber, 1);
    Vxc = zeros(Ng_fiber, 1);
    Vzc = zeros(Ng_fiber, 1);
    
    for i = 1:Ng_fiber
        Xc(i) = mean(X(i, :));
        Zc(i) = mean(Z(i, :));
        
        % Adjust for periodic boundaries
        if (Xc(i) < 0)
            X(i, :) = X(i, :) + Lx;
            Xc(i) = mean(X(i, :));
        elseif (Xc(i) > Lx)
            X(i, :) = X(i, :) - Lx;
            Xc(i) = mean(X(i, :));
        end
        
        if (Zc(i) < 0)
            Z(i, :) = Z(i, :) + Lz;
            Zc(i) = mean(Z(i, :));
        elseif (Zc(i) > Lz)
            Z(i, :) = Z(i, :) - Lz;
            Zc(i) = mean(Z(i, :));
        end
        
        % Calculate center velocities
        Vxc(i) = mean(Vx(i, :));
        Vzc(i) = mean(Vz(i, :));
    end
    
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
            
            dist(i, j) = sqrt(tmpx^2 + tmpz^2);
        end
    end
    
    % Create similarity matrix (Eq. 12)
    S = exp(-dist.^2);
    
    % Apply threshold (Eq. 13)
    S_eps = S;
    S_eps(S_eps < 0.7) = 0;
    
    % Perform spectral clustering
    G_eps = graph(S_eps);
    k = max(unique(conncomp(G_eps)));  % Number of clusters
    idx3 = spectralcluster(S_eps, k, 'Distance', 'precomputed');
    
    % Filter small clusters
    B = unique(idx3);
    Ncount = histc(idx3, B);
    
    % Print Ncount details
    fprintf('  Number of clusters before filtering: %d\n', length(Ncount));
    fprintf('  Cluster sizes before filtering:\n');
    for i = 1:length(Ncount)
        fprintf('    Cluster %d: %d fibers\n', B(i), Ncount(i));
    end
    
    for i = 1:Ng_fiber
        tmp = idx3(i);
        if (Ncount(tmp) < 11)  % If a cluster has less than 11 particles, ignore it
            idx3(i) = -1;
        end
    end
    
    % Get updated cluster IDs
    B = unique(idx3);  % Number of clusters after ignoring small ones
    Ncount_filtered = histc(idx3, B);
    
    % Print filtered Ncount details
    fprintf('  Number of clusters after filtering: %d\n', length(Ncount_filtered));
    fprintf('  Cluster sizes after filtering:\n');
    for i = 1:length(Ncount_filtered)
        fprintf('    Cluster %d: %d fibers\n', B(i), Ncount_filtered(i));
    end
    
    % Calculate cluster velocities
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
    
    % Calculate cluster centroids
    cluster_centroids = zeros(size(B, 1), 2);
    
    for i = 1:size(B, 1)
        cluster_fibers = find(idx3 == B(i));
        if ~isempty(cluster_fibers)
            cluster_centroids(i, 1) = mean(Xc(cluster_fibers));
            cluster_centroids(i, 2) = mean(Zc(cluster_fibers));
        end
    end
    
    % Plot clusters
    figure('Position', [400, 558, 260, 420]);
    
    % Generate colors for clusters
    try
        % Try to use linspecer if available
        color = linspecer(size(B, 1));
    catch
        % Fall back to default colormap if linspecer not available
        cmap = colormap(jet(size(B, 1)));
        color = cmap;
    end
    
    color(1, :) = 0;  % Ignored clusters will be shown in black
   

    %%
    % Plot fibers colored by cluster
for i = 1:size(B, 1)
    for j = 1:Ng_fiber
        if (idx3(j) == B(i))
            plot(X(j, :), Z(j, :), '.', 'MarkerSize', 1, 'Color', color(i, :));
            hold on;
        end
    end
end

% Set plot properties
axis equal;
axis([0, Lx, 0, Lz]);
title(sprintf('Fiber Clusters (AR = %d)', ar));

% Create custom tick marks and labels for x-axis (0, 2π, 4π)
xticks([0, 2*pi, 4*pi]);
xticklabels({'0', '2\pi', '4\pi'});

% Create custom tick marks and labels for y-axis (0, 2π, 4π, 6π, 8π)
yticks([0, 2*pi, 4*pi, 6*pi, 8*pi]);
yticklabels({'0', '2\pi', '4\pi', '6\pi', '8\pi'});

xlabel('x', 'Interpreter','latex');
ylabel('z', 'Interpreter','latex');

% Optional: Adjust font size if needed
set(gca, 'FontSize', 12);
    %%
    
    % Save figure
    saveas(gcf, fullfile(outputFolder, sprintf('clusters_AR%d_%d.png', ar, mid_step)));
    saveas(gcf, fullfile(outputFolder, sprintf('clusters_AR%d_%d.fig', ar, mid_step)));
    
    % Save cluster data
    clusters = struct();
    clusters.centroids = cluster_centroids;
    clusters.velocities = cluster_vel;
    clusters.counts = cluster_count;
    clusters.fiber_clusters = idx3;
    clusters.fiber_positions = [Xc, Zc];
    clusters.fiber_velocities = [Vxc, Vzc];
    clusters.Ncount = Ncount;          % Original cluster counts
    clusters.Ncount_filtered = Ncount_filtered; % Filtered cluster counts
    
    save(fullfile(outputFolder, sprintf('clusters_AR%d_%d.mat', ar, mid_step)), 'clusters');
    
    fprintf('  Analysis complete for AR = %d\n', ar);
    fprintf('  Found %d clusters after filtering\n', size(B, 1) - (B(1) == -1));  % Subtract 1 if isolated fibers exist
    fprintf('  Results saved to %s\n', outputFolder);
end

fprintf('All clustering analysis complete.\n');