clear all;
close all;
%Find files to analyze
folder = uigetdir;
files = dir( fullfile(folder,'*_out.xls') );
filenames = {};
for i = 1:length(files)
   filenames{i} = fullfile(folder,files(i).name); 
end
% filenames = find_files_recursive(folder,'\w*success+\w*.xls');
data_points = 2040;

% %Read in xy data
% %Transpose xy data
% curvature = menger_curve_array(data);
% 
% bends = zeros(1,size(curvature,2));
% for i = 1:size(curvature,2)
%     bends(i) = number_of_bends(curvature,hystresis_window_size);
% end
h = waitbar(0,'Processing files...');
emp = cell(10000,100);
count = zeros(data_points/2-1,1);
for i = 1:length(filenames)
    disp(filenames{i});
    data = xlsread(filenames{i},'skel');
    if(size(data,1)>data_points)
        data = data(1:data_points,:);
    end
    data_length = size(data,1)/2;
    waitbar(i/length(filenames),h);
    %calc velocity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = data(2:2:end,:);
    y = data(1:2:end-1,:);
    xd = diff(x,1,1);
    yd = diff(y,1,1);
    v = sqrt(xd.^2+yd.^2)*(400/150); %400um/pixel
    if(i==1)
        vt = zeros(data_length-1,size(v,2));
    end
    if(size(vt,1)~=size(v,1))
        disp([filenames{i},' does not have correct number of data points']);
        if(size(vt,1) > size(v,1))
            vt(1:size(v,1),:) = vt(1:size(v,1),:)+v;
        else
            vt = v(1:size(vt,1),:);
        end
    else
        vt = vt+v;
    end
    
    vf = [transpose(1:size(v,1)),v];
    xlswrite(filenames{i},emp,'Velocity');
    xlswrite(filenames{i},vf,'Velocity');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %calc curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = menger_curve_XY(transpose(x),transpose(y)).*(150/400);
    xlswrite(filenames{i},emp,'Curvature');
    k = transpose(k);
    xlswrite(filenames{i},k,'Curvature');
    [path,name,~] = fileparts(filenames{i});
%     data = transpose(k);
%     size(data)
%     saveas(plot_curvature(k(1:10:end,:),0,0,0), fullfile(path,[name,'_plot.png']));
    saveas(plot_curvature(k(1:10:end,:),0,0,0), fullfile(path,[name,'_plot_2.png']));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calc fourier transform for Curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ser = transpose(1:size(k,1));
    L = length(ser);
%     data_point = 2;
%     Fs = 10;            % Sampling frequency (frames per second)
%     T = 1/Fs;           % Sampling period
%     L = size(k,2);      % Length of signal
%     t = (0:L-1)*T;      % Time vector

    fft_data = fft(k,[],1);
    P2 = abs(fft_data/L);
%     Y = fft(k(data_point,1:L));
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = Fs*(0:(L/2))/L;
    xlswrite(filenames{i},emp,'FT Curvature');
    xlswrite(filenames{i},[ser,P2],'FT Curvature');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calc fourier transform for velocity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ser = transpose(1:size(v,1));
    L = length(ser);
%     data_point = 2;
%     Fs = 10;            % Sampling frequency (frames per second)
%     T = 1/Fs;           % Sampling period
%     L = size(k,2);      % Length of signal
%     t = (0:L-1)*T;      % Time vector

    fft_data = fft(v,[],1);
    P2 = abs(fft_data/L);
%     Y = fft(k(data_point,1:L));
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = Fs*(0:(L/2))/L;
    xlswrite(filenames{i},emp,'FT Velocity');
    xlswrite(filenames{i},[ser,P2],'FT Velocity');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %calc number of bends
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = transpose(k);
    bends = 1;
    sign_change = 0;
    hystresis_window_size = 2;
    bend_total = zeros(size(k,2),1);
    %initialize state values
    for z = 1:size(k,2)
        curr_state = k(1,z)/abs(k(1,z));
        prev_state = curr_state;
        for j = 2:size(k,1)
            curr_state = k(j,z)/abs(k(j,z));
            if(prev_state ~= curr_state)
                sign_change = sign_change + 1;
            end
            if(sign_change >= hystresis_window_size)
                bends = bends + 1;
                prev_state = curr_state;
                sign_change = 0;
            end
        end
        bend_total(z) = bends;
        bends = 1;
    end
    k = transpose(k);
    xlswrite(filenames{i},emp,'Number of bends');
    xlswrite(filenames{i},bend_total,'Number of bends');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %calc change in curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kdiff = abs(diff(k,1,1));
    xlswrite(filenames{i},emp,'Change Curvature');
    xlswrite(filenames{i},kdiff,'Change Curvature');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %calc total change in curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kdiff = sum(kdiff,2);
    if(i==1)
        kt = zeros(data_length-1,1);
    end
    if(size(kt,1)~=size(kdiff,1))
        kt(1:size(kdiff,1),:) = kt(1:size(kdiff,1),:)+kdiff;
        count(1:size(kdiff,1),:) = count(1:size(kdiff,1),:) + ones(size(kdiff,1),1);
    else
        kt = kt+kdiff;
        count = count + 1;
    end
    xlswrite(filenames{i},emp,'Change Curvature Total');
    xlswrite(filenames{i},kdiff,'Change Curvature Total');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%calc average velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(vt,2)
    vt(:,i) = vt(:,i)./count;
end
vt = [transpose(1:size(vt,1)),vt];
[path,name,~] = fileparts(folder);
xlswrite(fullfile(folder,[name,'.xls']),vt,'Average velocity');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calc average curvature change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kt = kt./count;
kt = [transpose(1:size(kt,1)),kt];
xlswrite(fullfile(folder,[name,'.xls']),kt,'Average curvature change');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(h)