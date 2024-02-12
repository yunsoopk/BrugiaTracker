% Clear all variables before starting the script
clear batch_processing_script; close all;

%% Get list of files from a folder
DefaultPath = 'C:\work\BrugiaAdultTracker\WormAssay_Video\';
root_folder_path = uigetdir(DefaultPath);
root_folder_path = [root_folder_path filesep];
files = dir([root_folder_path '**' filesep '*.avi']); % Ouputs struct with file details
files2 = dir([root_folder_path '**' filesep '*.mp4']); % Ouputs struct with file details
files = [files; files2];
disp(['Total number of files : ' num2str(size(files,1))]);

%% Global params
batch = 1;
percent_thresh = 0.9;
smoothing_coefficient = 8;
start_index = 1;
manual = 0;
show_video = 0;

not_processed = {};
counter = 1;
ellipse_points = linspace(0,2*pi, 50);

%% Loop through each file
for l = 1: size(files, 1)
    %% Read video file
    file_name = files(l).name; % Get file name
    if contains(file_name, '_out')
        continue
    end
    folder_path = [files(l).folder filesep];
    %% Save results to a folder
    results_path = [folder_path 'Results'];

    % Create results folder if it does not exist
    if exist(results_path,'dir') ~= 7
        mkdir(results_path);
    end
    disp(['File ' num2str(l) '/' num2str(size(files,1)) ': ' file_name]);
    vid = VideoReader([folder_path,file_name]);  %#ok<TNMLP>
    
    %% Save output data
    % Save tracks data as a video file
    out_vid_file = VideoWriter([results_path filesep ...
                           file_name(1:end-4) '_out']); %#ok<TNMLP>

    % Save parameters data into an excel file
    xlsfile = [results_path filesep file_name(1:end-4) '_out' '.xls'];
    worm_pos_file = [results_path filesep file_name(1:end-4) '_pos.txt'];
    well_pos_file = [results_path filesep file_name(1:end-4) '_well.txt'];
    worm_pos = [];
    
    % Get worm locations from text file
    if isfile(worm_pos_file)
        params = dlmread(worm_pos_file);
        worm_pos = params(1:2);
        percent_thresh = params(3);
    end
    
    % Worm location not stored for the file
    % Not possible to get input from the user in batch mode
    if isempty(worm_pos)
        not_processed(counter, :) = {[folder_path file_name]}; %#ok<*SAGROW>
        counter = counter + 1;
        continue;
    end
    
    % Set end index based on no of frames
    end_index = vid.NumberOfFrames; %#ok<*VIDREAD>
    end_index_xls = end_index + 1; % To include header on the top
    
    
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
    f = figure;
    % Open the video file for writing
    open(out_vid_file);
    
    %% Loop through each frame and extract the parameters
    for i = start_index:end_index
        img = read(vid,i);
%         disp(['Frame ' num2str(i)]);
        if size(img,3) == 3
            img = rgb2gray(img);
        end
        if(i>start_index)
            prev = bin; % Store previous binary image
        end
        threshold = mean2(img)*percent_thresh;
        % Arguments: image, radius_range, tolerance, manual
        if i == start_index
            % Read the well plate position from a file
            if isfile(well_pos_file)
                dims = dlmread(well_pos_file);
                center = dims(1:2);
                radius = dims(3);
            else
                [center, radius] = get_circular_well(img, ...
                                            [180 250], 10, manual, batch);
                % Write the well plate position to a file
                dlmwrite(well_pos_file, [center, radius]);
            end
            [rows,cols] = size(img);
            [xx, yy] = ndgrid((1:rows)-center(2), (1:cols)-center(1));
            mask = (xx.^2 + yy.^2) > radius^2;
        end
        if nnz(mask) == 0 % Handle circle not found cases
            not_processed(counter, :) = {[folder_path file_name]}; %#ok<*SAGROW>
            counter = counter + 1;
            continue;
        end
        % Set background to zero
        img(mask) = uint8(0);

        bin = img < threshold; % Apply percent threshold

        bin = bwareaopen(bin,700);
        if i == start_index
            % store borders to handle worm touching edge cases
            borders = bwboundaries(bin);
            % Get sizes of each boundary and sort them
            sizes_all = cellfun(@(x) size(x,1), borders); 
            [~, indices] = sort(sizes_all(:), 'descend');
            inner_boundary = [borders{indices(1)}; borders{indices(2)}];
        end
        % Set inner boundary pixels to zero
        for j = 1 : size(inner_boundary,1)
            bin(inner_boundary(j,1), inner_boundary(j,2)) = 0;
        end

        lab = bwlabel(bin);
        if(i==start_index) 
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
                break;
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
            disp(['Skipped frame' num2str(i)]);
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

        imshow(img);
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
end
% Save not processed files to text file
if ~isempty(not_processed)
    fid = fopen([folder_path 'not_processed.txt'], 'w');
    fprintf(fid, '%s\n', not_processed{:});
    fclose(fid);
end
disp([num2str(size(files,1)-size(not_processed,1)) ' file(s) processed']);
