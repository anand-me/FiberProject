%% Mean Vertical Velocity Analysis for Fiber Sedimentation
% This script extracts mean fiber velocity as a function of time from sedimentation data
% for multiple aspect ratios and creates comparative visualizations

%{
find /mnt/e/FibrePaper/AR10/n3/out_v3/ -name "tecp.dat031.*" | sort -V | head -n 1000
/mnt/d/OneDrive - Florida State University/Vahid_Material/FFT_IBM/suspension/medium/DEN_10_g_1000/S1_more_flexible/n3/out
%% AR 2 (done by Akshay)
%_____________________________________
tecp.dat031.00101 tecp.dat031.00200
tecp.dat031.02001 tecp.dat031.02200
tecp.dat031.22001 tecp.dat031.22203
%_____________________________________

%% AR 10 (done by Akshay)
%_____________________________________
tecp.dat031.00101 tecp.dat031.00200
tecp.dat031.02001 tecp.dat031.02500
tecp.dat031.25001 tecp.dat031.25204
%_____________________________________

%% AR 15 (done by Akshay)
%_____________________________________
tecp.dat031.00101 tecp.dat031.00200
tecp.dat031.02001 tecp.dat031.02500
tecp.dat031.25001 tecp.dat031.25228
%_____________________________________

%% AR 20 (done by Vahid)
%_____________________________________
tecp.dat031.00001 tecp.dat031.00050
tecp.dat031.00501 tecp.dat031.00600
%_____________________________________
%}




clear; close all;

% Domain dimensions
Lx = 4*pi;
Lz = 8*pi;

% Define aspect ratios to analyze
aspect_ratios = [2, 10, 15, 20];

% Node counts for each aspect ratio (as a cell array)
ar_nodes = [2, 10, 15, 20];
nodes_counts = [366, 359, 370, 405];

% Data directories for each aspect ratio (as a cell array)
data_folders = cell(length(aspect_ratios), 1);
data_folders{1} = 'E:/FibrePaper/AR2/n3/out_v3/';
data_folders{2} = 'E:/FibrePaper/AR10/n3/out_v3/';
data_folders{3} = 'E:/FibrePaper/AR15/n3/out_v3/';
data_folders{4} = 'D:\OneDrive - Florida State University\Vahid_Material\FFT_IBM\suspension\medium\DEN_10_g_1000\S1_more_flexible\n3\out_v3\';

% Define file ranges for each aspect ratio
% Format: [start1, end1; start2, end2; start3, end3]
file_ranges = cell(length(aspect_ratios), 1);
file_ranges{1} = [101, 200; 2001, 2200; 22001, 22203];  % AR 2
file_ranges{2} = [101, 200; 2001, 2500; 25001, 25204];  % AR 10
file_ranges{3} = [101, 200; 2001, 2500; 25001, 25228];  % AR 15
file_ranges{4} = [1, 50; 501, 600; NaN, NaN];          % AR 20

% Number of processors
N_proc = 32;

% File format
fmt = 'tecp.dat%03i.%05i';

% Store results
results = struct();
currentDir = pwd;

%% Process each aspect ratio
for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    fprintf('Processing Aspect Ratio %d...\n', ar);
    
    % Select the desired file range (using the middle range for consistency)
    range_idx = 3; % Using the middle range
    ar_idx = find(aspect_ratios == ar);
    range = file_ranges{ar_idx};
    N_start = range(range_idx, 1);
    N_end = range(range_idx, 2);
    
    % Create arrays to store velocities
    time_steps = N_start:N_end;
    mean_Vz = zeros(size(time_steps));
    
    % Cell array to store all individual fiber velocities at each timestep
    all_Vzc = cell(size(time_steps));
    
    % Data path
    ar_idx = find(aspect_ratios == ar);
    dataFolder = data_folders{ar_idx};
    cd(dataFolder);
    
    % Process each timestep
    for i = 1:length(time_steps)
        istep = time_steps(i);
        fprintf('AR=%d, Processing timestep %d (%d/%d)...\n', ar, istep, i, length(time_steps));
        
        Ng_fiber = 0;
        for i_proc = 0:N_proc-1
            % Create filename
            tmp = sprintf(fmt, i_proc, istep);
            
            % Check if file exists
            if ~exist(tmp, 'file')
                %fprintf('File not found: %s\n', tmp);
                continue;
            end
            
            % Open file
            fid_r = fopen(tmp, 'r');
            if fid_r < 0
                %warning('Cannot open file: %s', tmp);
                continue;
            end
            
            % Read header
            for j = 1:3
                str = fgets(fid_r);
            end
            
            % Parse header
            tmp_header = sscanf(str, 'ZONE N= %i, E= %i,F=FEPOINT, ET=QUADRILATERAL');
            if isempty(tmp_header) || length(tmp_header) < 2
                fclose(fid_r);
                continue;
            end
            
            N_tot = tmp_header(1);
            E_tot = tmp_header(2);
            
            % Read data
            try
                tmp_data = csvread(tmp, 3);
                fclose(fid_r);
            catch
                warning('Error reading data from file: %s', tmp);
                fclose(fid_r);
                continue;
            end
            
            % Check data length
            if size(tmp_data, 1) < N_tot
                warning('Incomplete data in file: %s', tmp);
                continue;
            end
            
            % Extract data for nodes
            data = tmp_data(1:N_tot, :);
            
            % Get nodes per element for this aspect ratio
            ar_idx = find(ar_nodes == ar);
            N = nodes_counts(ar_idx);
            
            % Check if total number of nodes is divisible by nodes per element
            if (rem(N_tot, N) ~= 0)
                warning('ERROR: total number of nodes is not divisible by node per element!');
                continue;
            end
            
            % Process each fiber
            N_fiber = N_tot/N; % Number of fibers in local processor
            for j = 1:N_fiber
                s1 = 1 + N*(j-1);
                s2 = N*j;
                i_fiber = j + Ng_fiber;
                
                % Allocate arrays if first processor
                if i_proc == 0 && j == 1
                    X = zeros(1000, N); % Pre-allocate with assumption of max 1000 fibers
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
                if s2 <= size(data, 1) && size(data, 2) >= 4
                    X(i_fiber, 1:s2-s1+1) = data(s1:s2, 1);
                    Z(i_fiber, 1:s2-s1+1) = data(s1:s2, 2);
                    Vx(i_fiber, 1:s2-s1+1) = data(s1:s2, 3);
                    Vz(i_fiber, 1:s2-s1+1) = data(s1:s2, 4);
                end
            end
            
            Ng_fiber = i_fiber;
        end
        
        % Calculate center values
        Vxc = zeros(Ng_fiber, 1);
        Vzc = zeros(Ng_fiber, 1);
        
        for j = 1:Ng_fiber
            Vxc(j) = mean(Vx(j, :), 'omitnan');
            Vzc(j) = mean(Vz(j, :), 'omitnan');
        end
        
        % Store mean velocities for this timestep
        mean_Vz(i) = mean(Vzc, 'omitnan');
        
        % Store individual fiber velocities
        all_Vzc{i} = Vzc;
        
        fprintf('AR=%d, Timestep %d: %d fibers, Mean Vz = %.4f\n', ar, istep, Ng_fiber, mean_Vz(i));
        
        % Clear large arrays for next iteration
        clear X Z Vx Vz Vxc Vzc;
    end
    
    % Store results for this aspect ratio
    results.(sprintf('AR%d', ar)) = struct('time_steps', time_steps, 'mean_Vz', mean_Vz, 'all_Vzc', {all_Vzc});
    
    % Plot and save for this aspect ratio
    figure('Position', [100, 100, 800, 600]);
    plot(time_steps, mean_Vz, 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
    title(sprintf('Mean Vertical Velocity of Fibers vs Time (AR=%d)', ar));
    xlabel('Time Step');
    ylabel('Mean Vertical Velocity');
    grid on;
    
    % Save figure
    saveas(gcf, sprintf('mean_vz_AR%d.png', ar));
    saveas(gcf, sprintf('mean_vz_AR%d.fig', ar));
    
    % Save data
    save(sprintf('mean_velocity_AR%d.mat', ar), 'time_steps', 'mean_Vz', 'all_Vzc');
end

% Return to original directory
cd(currentDir);

%% Create combined plot
figure('Position', [100, 100, 1000, 800]);

% Line style and color for each aspect ratio
line_styles = {'-r*', '-g*', '-b*', '-m*'};

% Plot vertical velocity for all aspect ratios
hold on;
legend_str = {};

for ar_idx = 1:length(aspect_ratios)
    ar = aspect_ratios(ar_idx);
    
    % Get data
    data = results.(sprintf('AR%d', ar));
    time_steps = data.time_steps;
    mean_Vz = data.mean_Vz;
    
    % Normalize time steps to start from 0
    norm_time_steps = time_steps - time_steps(1);
    
    % Plot
    plot(norm_time_steps, mean_Vz, line_styles{ar_idx}, 'LineWidth', 2, 'MarkerSize', 8);
    legend_str{end+1} = sprintf('AR %d', ar);
end

title('Mean Vertical Velocity of Fibers vs Time');
xlabel('Time Step');
ylabel('Mean Vertical Velocity');
legend(legend_str, 'Location', 'best');
grid on;

% Save combined figure
saveas(gcf, 'mean_vz_all_AR.png');
saveas(gcf, 'mean_vz_all_AR.fig');

% Save all results
save('mean_velocity_all_AR.mat', 'results');

fprintf('Analysis complete. Plots and data saved.\n');

%% plots

%% range_id = 3 for AR 2, 10, 15
figure(6)
load mean_velocity_all_AR_range_id_3.mat
plot(results.AR2.mean_Vz(1,end-200:end), '-b', 'LineWidth', 1.5); % 203 points
hold on;
plot(results.AR10.mean_Vz(1,end-200:end), '-k', 'LineWidth', 1.5); % 224 points
hold on;
plot(results.AR15.mean_Vz(1,end-200:end), '-r', 'LineWidth', 1.5); % 228 points
xlabel('Non-dimensional time', 'FontSize', 12);
ylabel('$\overline{V}_z$', 'Interpreter', 'latex', 'FontSize', 24);
legend('AR = 2', 'AR = 10', 'AR = 15', 'Location', 'best', 'FontSize', 12);
grid on;
xticks(0:20:200);
xlim([0 200]);
% Save outputs
saveas(gcf, './plots/mean_velocity_AR_2_10_15_range_id_3.jpg');
save('./plots/mean_velocity_AR_2_10_15_range_id_3.mat');
savefig('./plots/mean_velocity_AR_2_10_15_range_id_3.fig');

%% range_id = 2 for AR 2, 10, 15
figure(7)
load mean_velocity_all_AR_range_id_2.mat
plot(results.AR2.mean_Vz, '-b', 'LineWidth', 1.5); % 200 points
hold on;
plot(results.AR10.mean_Vz(1,301:end), '-k', 'LineWidth', 1.5); % 500 points
hold on;
plot(results.AR15.mean_Vz(1,301:end), '-r', 'LineWidth', 1.5); % 500 points
xlabel('Non-dimensional time', 'FontSize', 12);
ylabel('$\overline{V}_z$', 'Interpreter', 'latex', 'FontSize', 24);
legend('AR = 2', 'AR = 10', 'AR = 15', 'Location', 'best', 'FontSize', 12);
grid on;
% Save outputs
saveas(gcf, './plots/mean_velocity_AR_2_10_15_range_id_2.jpg');
save('./plots/mean_velocity_AR_2_10_15_range_id_2.mat');
savefig('./plots/mean_velocity_AR_2_10_15_range_id_2.fig');

%% range_id = 2, AR = 20
figure(8)
load mean_velocity_all_AR_range_id_2.mat
plot(results.AR20.mean_Vz, '-m', 'LineWidth', 1.5); % 100 points
xlabel('Non-dimensional time', 'FontSize', 12);
ylabel('$\overline{V}_z$', 'Interpreter', 'latex', 'FontSize', 24);
grid on;
% Save outputs
saveas(gcf, './plots/mean_velocity_AR_20_range_id_2.jpg');
save('./plots/mean_velocity_AR_20_range_id_2.mat');
savefig('./plots/mean_velocity_AR_20_range_id_2.fig');
