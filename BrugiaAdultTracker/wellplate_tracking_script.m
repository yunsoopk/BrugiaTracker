%% This script is used to run on Brugia Worm videos to extract data
% Magnifications used are 7.8x and 7.1x.
% The videos were recorded at 696 x 520 resolution for 4 minutes (~5fps).
% This corresponds to 53 pixels/mm for 7.8x and 47 pixels/mm for 7.1x.

% Clear all variables before starting the script
clear wellplate_tracking_script; close all;

%% Read video file
DefaultPath = 'C:\';
[file_name,folder_path] = uigetfile([DefaultPath filesep ...
                                    '*.avi;*.mp4']);
vid = VideoReader([folder_path,file_name]);

%% Save results to a folder
results_path = [folder_path 'Results'];

% Create results folder if it does not exist
if exist(results_path,'dir') ~= 7
    mkdir(results_path);
end

%% Save output data
% Save tracks data as a video file
out_vid_file = VideoWriter([results_path filesep ...
                           file_name(1:end-4) '_out']);

% Save parameters data into an excel file
xlsfile = [results_path filesep file_name(1:end-4) '_out' '.xls'];
worm_pos_file = [results_path filesep file_name(1:end-4) '_pos.txt']
well_pos_file = [results_path filesep file_name(1:end-4) '_well.txt'];

%% Global params
batch = 0;
percent_thresh = 0.9;
smoothing_coefficient = 8;
start_index = 1;
end_index = vid.NumberOfFrames; %#ok<*VIDREAD>
end_index_xls = end_index + 1; % To include header on the top
manual = 0;
show_video = 1;
worm_pos = [];
ellipse_points = linspace(0,2*pi, 50);

% Get worm locations from text file
if isfile(worm_pos_file)
    params = dlmread(worm_pos_file);
    worm_pos = params(1:2);
    percent_thresh = params(3);
end

% Create arrays to save data for each parameter
centroids = cell(end_index_xls,2); % Save centroid information
% Save ellipse info (major, minor, orientation, eccentricity)
axis_info = cell(end_index_xls,4);
misc_info = cell(end_index_xls,2); % Save Euler no and extent
areas = cell(end_index_xls,1); % Save area information

% Add header information to the data
centroids(1,:) = {'X', 'Y'};
misc_info(1,:) = {'Euler Number', 'Extent'};
areas(1) = {'Areas'};
axis_info(1,:) = {'Major', 'Minor', 'Angle', 'Eccentricity'};

included_frames = [];
% f = figure;
% Open the video file for writing
open(out_vid_file);
%% Loop through each frame and extract the parameters
for i = start_index:end_index
    img = read(vid,i);
    disp(['Frame ' num2str(i)]);
    if size(img,3) == 3
        img = rgb2gray(img);
    end
    if(i>start_index)
        prev = bin; % Store previous binary image
    end
    threshold = mean2(img)*percent_thresh;
%     % Arguments: image, radius_range, tolerance, manual
%     if i == start_index
%         % Read the well plate position from a file
%         if isfile(well_pos_file)
%             dims = dlmread(well_pos_file);
%             center = dims(1:2);
%             radius = dims(3);
%         else
%             [center, radius] = get_circular_well(img, ...
%                                             [180 250], 10, manual, batch);
%             % Write the well plate position to a file
%             dlmwrite(well_pos_file, [center, radius]);
%         end
%         [rows,cols] = size(img);
%         [xx, yy] = ndgrid((1:rows)-center(2), (1:cols)-center(1));
%         mask = (xx.^2 + yy.^2) > radius^2;
%     end
%     % Set background to zero
%     img(mask) = uint8(0);
% %     if i == start_index
% %         threshold = (sumabs(img)/nnz(img))*percent_thresh;
% %     end
    bin = img < threshold; % Apply percent threshold
    bin = bwareaopen(bin,700);
    if i == start_index
        % store borders to handle worm touching edge cases
        borders = bwboundaries(bin);
        % Get sizes of each boundary and sort them
        sizes_all = cellfun(@(x) size(x,1), borders); 
        [~, indices] = sort(sizes_all(:), 'descend');
        % Get outer image and well boundaries to handle worm touching
        % the well cases
        inner_boundary = [borders{indices(1)}; borders{indices(2)}];
    end
    % Set boundaries pixels to zero to separate well from worm blobs
    for j = 1 : size(inner_boundary,1)
        bin(inner_boundary(j,1), inner_boundary(j,2)) = 0;
    end
    
    lab = bwlabel(bin);
    if(i==start_index)
        if isempty(worm_pos)
            RGB = insertText(lab,[10 10],'Select the worm','FontSize',...
                            18,'BoxOpacity',0.4,'TextColor','black');
            figure; imshow(img);
            fig2 = figure;
            imshow(RGB);
            worm_pos = ginput(1);
            close all;
            % Save for future use
            dlmwrite(worm_pos_file, [worm_pos, percent_thresh]);
        end
        bin = lab==lab(round(worm_pos(2)),round(worm_pos(1)));
    else
        U = lab.*double(prev);
        U = unique(U(:));
        U(U==0) = [];
        if size(U,1) > 1
            vals = [];
            for k = 1: size(U,1)
                vals = [vals nnz(lab == U(k))];
            end
            [~, idx] = max(vals);
            U(U~=U(idx)) = [];
        end
        if(isempty(U))
            disp(['NoWormFoundException: Skipped frame' num2str(i)]);
            continue;
        end
        bin = lab==U;
    end
    results = regionprops(bin, uint8(bin), 'Centroid', 'Area',...
                          'Extent', 'EulerNumber',...
                          'MajorAxisLength', 'MinorAxisLength',...
                          'Orientation', 'Eccentricity');
    if i == start_index
        worm_area = results.Area;
    end
    if results.Area > 4*worm_area
        disp(['WormAreaException: Skipped frame' num2str(i)]);
        continue
    end
    
    % Store information
    included_frames = [included_frames i+1]; %#ok<*AGROW>
    centroids(i+1, :) = {results.Centroid(1), results.Centroid(2)};
    misc_info(i+1, :) = {results.EulerNumber, results.Extent};
    areas(i+1) = {results.Area};
    axis_info(i+1, :) = {results.MajorAxisLength,...
                        results.MinorAxisLength,...
                        results.Orientation, results.Eccentricity};
    
    % Plot centroid tracks
    if ~show_video
        f.Visible = 'off';
    end
    
    % Get minEllipse coordinates to draw on the image
    a = results.MajorAxisLength/2;
    b = results.MinorAxisLength/2;
    x_c = results.Centroid(1);
    y_c = results.Centroid(2);
    theta = deg2rad(-results.Orientation);
    x = x_c + a*cos(ellipse_points)*cos(theta) ...
            - b*sin(ellipse_points)*sin(theta);
    y = y_c + a*cos(ellipse_points)*sin(theta) ...
            + b*sin(ellipse_points)*cos(theta);
    
    imshow(bin);
    hold on;
    % Plot centroid position
    plot(cell2mat(centroids(2:i+1,1)),cell2mat(centroids(2:i+1,2)), ...
         'LineWidth', 1, 'Color', 'b');
     % Draw minEllipse of the worm
    plot(x,y,'r','Linewidth',1)
    hold off;
    pause(0.01);
    frame = getframe(gcf);
    % Write the out frame
    writeVideo(out_vid_file,frame);
end

%% Close video file
close(out_vid_file);

%% Store extracted information to files
% Remove skipped frames from the data
if size(included_frames,2) ~= end_index
    centroids = centroids(included_frames,:);
    misc_info = misc_info(included_frames,:);
    areas = areas(included_frames);
    axis_info = axis_info(included_frames,:);
end

% Write the lists to excel sheets
xlswrite(xlsfile,centroids,'Centroids');
xlswrite(xlsfile,misc_info,'Misc Info');
xlswrite(xlsfile,areas,'Areas');
xlswrite(xlsfile,axis_info,'Axis Info');
disp(['Processed file: ' file_name]);
