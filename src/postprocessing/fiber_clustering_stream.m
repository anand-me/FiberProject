%% Mean Velocity Analysis for Fiber Sedimentation


%% Configuration
clear; close all;
Lx = 4*pi;
Lz = 8*pi;
N = 359+11; % # of nodes per element AR=20-> 405  | AR=15-> 370 | AR=10-> 359 | AR=2-> 366
N_proc = 32;
N_start = 2081;
N_end = 2100;

% Create array to store mean velocities
time_steps = N_start:N_end;
mean_Vz = zeros(size(time_steps));
mean_Vx = zeros(size(time_steps));

%% Data path
dataFolder = 'E:/FibrePaper/AR15/n3/out_v3/';
currentDir = pwd;
cd(dataFolder);

%% DATA FORMAT  
fmt = 'tecp.dat%03i.%05i';

%% Process each timestep
for i = 1:length(time_steps)
    istep = time_steps(i);
    fprintf('Processing timestep %d...\n', istep);
    
    Ng_fiber = 0;
    for i_proc = 0:N_proc-1
        tmp = sprintf(fmt, i_proc, istep);
        filename1 = fullfile('out_v3', tmp);
        
        % Skip if file doesn't exist
        if ~exist(filename1, 'file')
            fprintf('File not found: %s\n', filename1);
            continue;
        end
        
        fid_r = fopen(filename1, 'r');
        if fid_r < 0
            warning('Cannot open file: %s', filename1);
            continue;
        end
        
        % Read header
        for j = 1:3
            str = fgets(fid_r);
        end
        
        tmp = sscanf(str, 'ZONE N= %i, E= %i,F=FEPOINT, ET=QUADRILATERAL');
        N_tot = tmp(1);
        E_tot = tmp(2);
        
        try
            tmp = csvread(filename1, 3);
            fclose(fid_r);
        catch
            warning('Error reading data from file: %s', filename1);
            fclose(fid_r);
            continue;
        end
        
        data = tmp(1:N_tot, :);
        if (rem(N_tot, N) ~= 0)
            warning('ERROR: total number of nodes is not divisible by node per element!');
            continue;
        end
        
        N_fiber = N_tot/N; % in local proc
        for j = 1:N_fiber
            s1 = 1 + N*(j-1);
            s2 = N*j;
            i_fiber = j + Ng_fiber;
            
            X(i_fiber, :) = data(s1:s2, 1);
            Z(i_fiber, :) = data(s1:s2, 2);
            Vx(i_fiber, :) = data(s1:s2, 3);
            Vz(i_fiber, :) = data(s1:s2, 4);    
        end
        
        Ng_fiber = i_fiber;
    end
    
    % Calculate center values
    for j = 1:Ng_fiber
        Vxc(j) = mean(Vx(j, :));
        Vzc(j) = mean(Vz(j, :));
    end
    
    % Store mean velocities for this timestep
    mean_Vx(i) = mean(Vxc);
    mean_Vz(i) = mean(Vzc);
    
    fprintf('Timestep %d: %d fibers, Mean Vz = %.4f\n', istep, Ng_fiber, mean_Vz(i));
    
    % Clear arrays for next iteration
    clear X Z Vx Vz Vxc Vzc;
end

% Return to original directory
cd(currentDir);

%% Plot mean fiber velocity vs time
figure('Position', [100, 100, 800, 600]);

% Plot vertical velocity
subplot(2, 1, 1);
plot(time_steps, mean_Vz, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
title('Mean Vertical Velocity of Fibers vs Time');
xlabel('Time Step');
ylabel('Mean Vertical Velocity');
grid on;

% Plot horizontal velocity
subplot(2, 1, 2);
plot(time_steps, mean_Vx, 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
title('Mean Horizontal Velocity of Fibers vs Time');
xlabel('Time Step');
ylabel('Mean Horizontal Velocity');
grid on;

% Save figure
saveas(gcf, 'mean_velocity_vs_time.png');
saveas(gcf, 'mean_velocity_vs_time.fig');

% Save data
save('mean_velocity_data.mat', 'time_steps', 'mean_Vx', 'mean_Vz');

fprintf('Analysis complete. Plots and data saved.\n');