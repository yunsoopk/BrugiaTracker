% Clear all variables before starting the script
clear batch_analyze_data; close all;

%% Get list of files from a folders and subfolders
DefaultPath = 'Z:\Brugia Videos (TRS) (Adult Male)\';
root_folder_path = uigetdir(DefaultPath);
root_folder_path = [root_folder_path filesep];
files = dir([root_folder_path '**' filesep '*.xls']);
disp(['Total number of files : ' num2str(size(files,1))]);

%% Global params
% save_plots = 1;
fig_resolution = 600;
time_difference = 1;
paths = split(root_folder_path, filesep);
path = strjoin(paths(1:3), filesep);
xlsfile = strjoin([path filesep paths(3) '_all' '.xlsx'], '');

% Create templates for root concentrations
cen_vel_all = zeros(1020, size(files,1)*3);
cen_vel_header = {};
ang_vel_all = zeros(1020, size(files,1)*2);
ang_vel_header = {};
eul_rate_all = zeros(1020, size(files,1)*2);
eul_rate_header = {};
ext_rate_all = zeros(1020, size(files,1)*2);
ext_rate_header = {};
ecc_rate_all = zeros(1020, size(files,1)*2);
ecc_rate_header = {};
curv_all = zeros(1020, size(files,1)*3);
curv_header = {};

f = figure;
%% Loop through each file
for loop = 1: size(files, 1)
    %% Read video file
    file_name = files(loop).name; % Get file name
    folder_path = [files(loop).folder filesep];
    disp(['File ' num2str(loop) '/' num2str(size(files,1)) ': ' file_name]);
    
     %% Add headers to the excel sheets
    cen_vel_header = [cen_vel_header {'X' ,'Y', file_name}];
    ang_vel_header = [ang_vel_header {'Angle', file_name}]; %#ok<*AGROW>
    eul_rate_header = [eul_rate_header {'Euler No', file_name}];
    ext_rate_header = [ext_rate_header {'Extent', file_name}];
    ecc_rate_header = [ecc_rate_header {'Eccentricity', file_name}];
    curv_header = [curv_header {'X' ,'Y', file_name}];
    
    %% Save centroid velocity plot
    data = xlsread([folder_path file_name], 'Centroids');
    end_index_xls = size(data, 1);
    x = transpose(1:end_index_xls);
    derivatives = calc_derivative(data, time_difference);
    % Calculate euclid distance of the derivatives
    centroid_velocities = sqrt(derivatives(:,1).^2+derivatives(:,2).^2);
    
    % copy to main data
    cen_vel_all(1:end_index_xls, 3*loop-2:3*loop-1) = data;
    cen_vel_all(1:1:size(centroid_velocities,1), 3*loop) = centroid_velocities;
    
    save_plot([folder_path file_name(1:end-4) '_cen_vel.png'], ...
               x(1:end-time_difference), centroid_velocities, ...
               'Centroid Velocity', ...
               fig_resolution, 0);
           
   % Calculate centroid curvature
   curv_data = transpose(menger_curve(data));
   curv_data = curv_data(2:end-1);
   
   % copy to main data
    curv_all(1:end_index_xls, 3*loop-2:3*loop-1) = data;
    curv_all(1:1:size(curv_data,1), 3*loop) = curv_data;
           
    %% Save Ellipse info plots
    data = xlsread([folder_path file_name], 'Axis Info');
    
    % Eccentricity
    ecc_data = data(:,4);
    save_plot([folder_path file_name(1:end-4) '_ecc.png'], ...
               x, ecc_data, 'Eccentricity', ...
               fig_resolution, 0);
           
    % Eccentricity derivative
    ecc_derv = calc_derivative(ecc_data, time_difference);
    save_plot([folder_path file_name(1:end-4) '_ecc_derv.png'], ...
               x(1:end-time_difference), ecc_derv, ...
               'Change in eccentricity', ...
               fig_resolution, 0);
           
    % Copy to main excel sheet
    ecc_rate_all(1:end_index_xls, 2*loop-1) = ecc_data;
    ecc_rate_all(1:size(ecc_derv,1), 2*loop) = ecc_derv;
    
    % Worm angle plot (Scatter)
    angle_data = data(:,3);
    save_plot([folder_path file_name(1:end-4) '_angles.png'], ...
               x, angle_data, 'Angle', ...
               fig_resolution, 1);
    
    % Worm angle derivative
    angle_derv = calc_derivative(angle_data, time_difference);
    save_plot([folder_path file_name(1:end-4) '_angles_derv.png'], ...
               x(1:end-time_difference), angle_derv, ...
               'Change in angles', ...
               fig_resolution, 0);
           
    % Copy to main excel sheet
    ang_vel_all(1:end_index_xls, 2*loop-1) = angle_data;
    ang_vel_all(1:size(angle_derv,1), 2*loop) = angle_derv;
           
    %% Save Misc info plots
    data = xlsread([folder_path file_name], 'Misc Info');
    
    % Euler number
    eul_data = data(:, 1);
    save_plot([folder_path file_name(1:end-4) '_euler.png'], ...
               x, eul_data, 'Euler Number', ...
               fig_resolution, 0);
           
    % Euler number derivative
    eul_derv = calc_derivative(eul_data, time_difference);
    save_plot([folder_path file_name(1:end-4) '_euler_derv.png'], ...
               x(1:end-time_difference), eul_derv, ...
               'Change in Euler Number', ...
               fig_resolution, 0);
           
    % Copy to main excel sheet
    eul_rate_all(1:end_index_xls, 2*loop-1) = eul_data;
    eul_rate_all(1:size(eul_derv,1), 2*loop) = eul_derv;
           
    % Extent
    ext_data = data(:, 2);
    save_plot([folder_path file_name(1:end-4) '_extent.png'], ...
               x, ext_data, 'Extent', ...
               fig_resolution, 0);
           
    % Euler number derivative
    ext_derv = calc_derivative(ext_data, time_difference);
    save_plot([folder_path file_name(1:end-4) '_extent_derv.png'], ...
               x(1:end-time_difference), ext_derv, 'Change in Extent', ...
               fig_resolution, 0);
           
    % Copy to main excel sheet
    ext_rate_all(1:end_index_xls, 2*loop-1) = ext_data;
    ext_rate_all(1:size(ext_derv,1), 2*loop) = ext_derv;
        
end

%% Save excel sheets
% Add headers to the sheets
cen_vel_all = [cen_vel_header; num2cell(cen_vel_all)];
ang_vel_all = [ang_vel_header; num2cell(ang_vel_all)];
eul_rate_all = [eul_rate_header; num2cell(eul_rate_all)];
ext_rate_all = [ext_rate_header; num2cell(ext_rate_all)];
ecc_rate_all = [ecc_rate_header; num2cell(ecc_rate_all)];
curv_all = [curv_header; num2cell(curv_all)];

% Write to excel sheets
xlswrite(xlsfile,cen_vel_all,'Centroid Velocity');
xlswrite(xlsfile,ang_vel_all,'Angular Velocity');
xlswrite(xlsfile,eul_rate_all,'Rate of Euler Number');
xlswrite(xlsfile,ext_rate_all,'Rate of Extent');
xlswrite(xlsfile,ecc_rate_all,'Rate of Eccentricity');
xlswrite(xlsfile,curv_all,'Curvature');

disp([num2str(size(files,1)) ' file(s) processed']);