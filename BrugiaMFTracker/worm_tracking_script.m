% Clear variables and close all windows
clear all; close all;

% Read video file
[file,path] = uigetfile('*');
vid = VideoReader([path,file]);

%visualize evolution
visual_img = 1;
visual_cost = 0;

%number of points
num_points = 75;
smoothing_coefficient = 8;

% initialize params
cost = zeros(vid.NumberOfFrames,1);
thresh = 100;
start_index = 1;%round(vid.NumberOfFrames/4);
end_index = vid.NumberOfFrames;

% Initialize video file for writing
out=regexp(path(1:end-1),filesep,'split');
results_path = '';
for i = 1:length(out)
    if strcmp(results_path,'') == 1
        results_path = out{i};
    else
        if strcmp(out{i},'Brugia') == 1
            results_path = [results_path filesep 'Results'];
        else
            results_path = [results_path filesep out{i}];
        end
    end
end

% Create results folder if it does not exist
if exist(results_path,'dir') ~= 7
    mkdir(results_path);
end

% Set out video and excel file names and open them for writing
out_vid_file = [results_path filesep file(1:end-4) '_out' '.mp4'];
vido=VideoWriter(out_vid_file,'MPEG-4');
vido.FrameRate = vid.FrameRate;
open(vido);

% Initialize the excel file for data logging
xlsfile = [results_path filesep file(1:end-4) '_out' '.xls'];
xlsdatacenter = [];
xlsdataleft = [];
xlsdataright = [];

for i = start_index:end_index
    
    img = read(vid,i);
    if(i>start_index)
        prev = bin;
        skel_prev = skel;
    end
    %bin = (background(:,:,1)-img(:,:,1))>thresh;
    %thresh = mean(mean(img(:,:,1)))/1.5;
    %bin = ~bradley(img(:,:,1), [45 45], 10);
    bin = img(:,:,1)<mean(mean(img(:,:,1)))*.9;
%     if(i==start_index)
%         bin = imfill(bin,'holes');
%     end
    bin = imclearborder(bin);
    bin = bwareaopen(bin,700);
    bin = ~bwareaopen(~bin,30);
    %bin = imfill(bin,'holes');
    lab = bwlabel(bin);
    if(i==start_index) 
        RGB = insertText(lab,[10 10],'Select the worm','FontSize',18,...
            'BoxOpacity',0.4,'TextColor','white');
        imshow(RGB)
        [col,row] = ginput(1);
        bin = lab==lab(round(row),round(col));
    else
        U = lab.*double(prev);
        U = unique(U(:));
        U(U==0) = [];
        if(isempty(U))
            break;
        end
        bin = lab==U;
    end
    borders = bwboundaries(bin,8);
    
    borders = smoothBorders(borders,smoothing_coefficient);
    d = borders{1};
    k = menger_curve_bound(d);
    if(i==start_index)
        RGB = insertText(img,[10 10],'Select tail and head',...
            'FontSize',18, 'BoxOpacity',0.4,'TextColor','white');
        imshow(RGB)
        [ct,rt] = ginput(1);
        [ch,rh] = ginput(1);
        tic;
    end
    %k = transpose(-k);
    %kt = abs(k-1).^-1/max(abs(k-1).^-1);
    %kh = abs(k-.5).^-1/max(abs(k-.5).^-1);
    %k = k/max(k);
    %{
    dist_tail = sqrt((d(:,2)-ct).^2+(d(:,1)-rt).^2);
    if(min(dist_tail)==0)
        dist_tail(dist_tail==0) = .000001;
    end
    mt = dist_tail<20;
    dist_tail = dist_tail.^-1;
    dist_tail = dist_tail/max(dist_tail);
    dist_head = sqrt((d(:,2)-ch).^2+(d(:,1)-rh).^2);
    if(min(dist_head)==0)
        dist_head(dist_head==0) = .000001;
    end
    mh = dist_head<20;
    dist_head = dist_head.^-1;
    dist_head = dist_head/max(dist_head);
    %}
    %cost_tail = (transpose(k)+1)*100 + sqrt((d(:,2)-ct).^2+(d(:,1)-rt).^2);
    %cost_head = (transpose(k)+.5)*100 + sqrt((d(:,2)-ch).^2+(d(:,1)-rh).^2);
    %didn't fix anything
    if(i==start_index)
        D = sqrt((d(:,2)-ct).^2+(d(:,1)-rt).^2);
        cost_tail = abs((transpose(k)+1))*100 + D;
        cost_tail = cost_tail + 100000*double(D>40);
        D = sqrt((d(:,2)-ch).^2+(d(:,1)-rh).^2);
        cost_head = abs((transpose(k)+.5))*100 + D;
        cost_head = cost_head + 100000*double(D>40);
    else 
        D = sqrt((d(:,2)-ct).^2+(d(:,1)-rt).^2);
        cost_tail = abs((transpose(k)+1))*100 + D;
        cost_tail = cost_tail + 100000*double(D>60);
        %D = sqrt((d(:,2)-skel(2,2)).^2+(d(:,1)-skel(2,1)).^2)*2;
        %cost_tail = cost_tail + D;
        D = sqrt((d(:,2)-ch).^2+(d(:,1)-rh).^2);
        cost_head = abs((transpose(k)+.5))*100 + D;
        cost_head = cost_head + 100000*double(D>60);
        %D = sqrt((d(:,2)-skel(end-1,2)).^2+(d(:,1)-skel(end-1,1)).^2)*2;
        %cost_head = cost_head + D;
    end
    %cost_tail = (k + .05*dist_tail).*double(mt);
    %cost_head = (k + .05*dist_head).*double(mh);
    ct = d(cost_tail==min(cost_tail),2);
    rt = d(cost_tail==min(cost_tail),1);
    ch = d(cost_head==min(cost_head),2);
    rh = d(cost_head==min(cost_head),1);
    
    if(rh==rt&&ch==ct)
        disp('trouble')
        skel(end,:) = skel(end-1,:);
        skel(1,:) = skel(2,:);
        ct = (ct+skel(1,2)*4)/5;
        rt = (rt+skel(1,1)*4)/5;
        ch = (ch+skel(end,2)*4)/5;
        rh = (rh+skel(end,1)*4)/5;
    end
    
%     if(visual_cost)
%         clf('reset');
%         plot(cost_tail);
%         hold on;
%         plot(cost_head,'g');
%         drawnow;
%     end
    
    %if(length(borders)>1)
    if(i>start_index)
        prev_skel = skel;
        skel(1,:) = [rt,ct];
        skel(end,:) = [rh,ch];
        %{
        [side_1,side_2] = splitBound(d,[rh,ch],[rt,ct]);
       
        t1 = side_2;
        for j = 2:length(borders)
            t1 = [t1;borders{j}];
        end
        dist1 = sqrt(sum((skel(1,:)-side_1(1,:)).^2,2));
        dist2 = sqrt(sum((skel(end,:)-side_1(1,:)).^2,2));
        if(dist1<dist2 && dist1<50)
            skel(1,:) = side_1(1,:);
            skel(end,:) = side_1(end,:);
        elseif(dist2<50)
            skel(end,:) = side_1(1,:);
            skel(1,:) = side_1(end,:);
        else
            disp('trouble');
        end
        if(sum(sum(isnan(left)))>0||sum(sum(isnan(right)))>0)
            disp('error');
        end
        %}
        if(i>3060)
            [left,right] = moveSkel(skel,widths,borders,bin,1);
        else
            [left,right] = moveSkel(skel,widths,borders,bin,0);
        end
        skel = (left+right)/2;
        skel(2,:) = (skel(1,:)+skel(3,:))/2;
        skel(end-1,:) = (skel(end,:)+skel(end-2,:))/2;
        %%
        [rows,cols] = find(bin);
        X = [min(cols), max(cols)];
        Y = [min(rows), max(rows)];
        %{
        subplot(2,1,1);
        imshow(img(Y(1):Y(2),X(1):X(2),:))%,'border','tight');
        hold on;
        plot(skel(:,2)-X(1),skel(:,1)-Y(1),'y-','linewidth',2);
        %}
        %%
        [skel,interval] = distributePoints(skel);
        skel = smoothSkel(skel,.3);
        [left,right] = genEdges(skel,widths/2);
    else
        %{
        %if(i>200)
        [side_1,side_2] = splitBound(d,[rh,ch],[rt,ct]);
        dist1 = sqrt(sum((skel(1,:)-side_1(1,:)).^2,2));
        dist2 = sqrt(sum((skel(end,:)-side_1(1,:)).^2,2));
        if(dist1<dist2 && dist1<50)
            skel(1,:) = side_1(1,:);
            skel(end,:) = side_1(end,:);
        elseif(dist2<50)
            skel(end,:) = side_1(1,:);
            skel(1,:) = side_1(end,:);
        end
        [left,right] = moveSkel2(skel,widths,borders,bin,0);
        skel = (left+right)/2;
        skel = distributePoints(skel);
        %else
        %}
        if(i==235)
            disp('stopping');
        end
        [side_1,side_2] = splitBound(d,[rh,ch],[rt,ct]);
        
        if(i>start_index)
            prev_skel = skel;
            [skel,widths] = buildSkel(side_1,side_2,num_points);
            d1 = sum(sqrt(sum((skel-prev_skel).^2,2)));
            d2 = sum(sqrt(sum((flipud(skel)-prev_skel).^2,2)));
            if(d2<d1)
                skel = flipud(skel);
            end
        else
            [skel,widths] = buildSkel(side_1,side_2,num_points);
            %widths(2) = widths(2)*2;
            %widths(end-1) = widths(end-1)*2;
        end
        [skel,interval1] = distributePoints(skel);

        %end

        %skel = smoothSkel(skel,.25);
        [left,right] = genEdges(skel,widths/2);
    end
    
    %img(:,:,1) = img(:,:,1)+uint8(perim)*255;
%     if(sd>.9)
%         img(:,:,3) = img(:,:,3)+uint8(bwmorph(occ,'dilate',1))*255;
%     end
    
    %optical flow
    if(i<start_index)
        img(:,:,1) = img(:,:,1)+uint8(prev-bin)*255;
        img(:,:,3) = img(:,:,3)+uint8(bin-prev)*255;
    end
    if(i==start_index)
        set(gcf, 'Position', get(0,'Screensize'));
    end
    if(visual_img)
        clf;
        %[rows,cols] = find(bin);
        %X = [min(cols), max(cols)];
        %Y = [min(rows), max(rows)];
        %subplot(2,1,1);
        %imshow(img(Y(1):Y(2),X(1):X(2),:))%,'border','tight');
        % Display frame no on the image
        img = insertText(img,[10 10],['Frame No :' num2str(i)],...
            'FontSize',18,'BoxOpacity',0.4,'TextColor','white');
        imshow(img,'border','tight');
        hold on;
        %{
        norm = calcNorms(borders);
        for j = 1:length(borders)
            t = borders{j};
            plot(t(:,2),t(:,1),'r-');
            for z = 1:10:length(t)
                plot([t(z,2);t(z,2)+10*norm{j}(z,1)],[t(z,1);t(z,1)+10*norm{j}(z,2)],'b-');
            end
        end
        %}
        %plot(side_1(:,2),side_1(:,1),'c-','linewidth',2);
        plot(left(:,2),left(:,1),'c-','linewidth',1);
        plot(right(:,2),right(:,1),'c-','linewidth',1);
        plot(skel(:,2),skel(:,1),'y-','linewidth',2);
        plot(skel(1,2),skel(1,1),'ro','linewidth',1);
        plot(skel(end,2),skel(end,1),'bo','linewidth',1);
        %{
        clf('reset');
        plot(cost_head);
        hold on;
        plot(cost_tail,'r-');
        %}
        %plot(side_2(:,2),side_2(:,1),'g-','linewidth',2);
        %plot(skel(:,2)-X(1),skel(:,1)-Y(1),'y-','linewidth',2);
        %plot(ct-X(1),rt-Y(1),'ro');
        %plot(ch-X(1),rh-Y(1),'ko');
        %curvature = abs(menger_curve(skel));
        %[v,ind] = max(curvature(:));
        %plot(skel(ind,2),skel(ind,1),'go');
        %plot(skel(ind,2),skel(ind,1),'g*');
        title(i)
        %title([num2str(min(cost_tail)),' ',num2str(min(cost_head))])
        %set(gcf, 'Position', get(0,'Screensize'));
        %title(sum(sum(sqrt(diff(skel,1,1).^2))))
        %title(length(borders))
        frame = getframe(gcf);
        writeVideo(vido,frame);
        drawnow;
        %{
        vects = diff(skel,1,1);
        subplot(2,1,2);
        plot(vects(:,1)/interval1,'bo');
        hold on;
        plot(vects(:,2)/interval1,'ro');
        pause(.01);
        %}
        xlsdatacenter = [xlsdatacenter;[{'X'},num2cell(transpose(skel(:,2)))];[{'Y'},num2cell(transpose(skel(:,1)))]];
        xlsdataleft = [xlsdataleft;{'X'},num2cell(transpose(left(:,2)));{'Y'},num2cell(transpose(left(:,1)))];
        xlsdataright = [xlsdataright;{'X'},num2cell(transpose(right(:,2)));{'Y'},num2cell(transpose(right(:,1)))];
    end
    cost(i) = min(cost_head);

end
toc;
close(vido);
xlswrite(xlsfile,xlsdatacenter,'skel');
xlswrite(xlsfile,xlsdataleft,'left');
xlswrite(xlsfile,xlsdataright,'right');