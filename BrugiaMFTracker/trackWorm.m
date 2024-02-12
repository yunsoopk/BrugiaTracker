function [ points,left,right ] = trackWorm( points, binary )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%clear all
%close all

%visualize evolution
visual = 0;
max_distance = 3;

%number of points
seperation = 15;
max_iterations = 300;

grav = bwperim(binary,8);
[row,col] = find(grav);
pointsg = [row,col];

%weights
alpha = .2;
beta = 0.5;
gamma = 0.0001;

%constraints
distance_factor = 1;
width = 3;

%start timer
tic;
for i = 1:max_iterations
    %record previous points
    prev_points = points;
    
    %this calculates the distances
    D = points(1:end-1,:)-points(2:end,:);
    dist = sqrt(sum(D.^2,2))-seperation;
    %tring to stop oscillations
    dist = dist/sum(abs(dist)) .* double(abs(dist)>distance_factor);
    dist = [0,0;D.*[dist,dist]];
    seperation = mean(sqrt(sum(D.^2,2)));
    
    %this calculates the smoothness
    smooth = points(1:end-2,:)-2*points(2:end-1,:)+points(3:end,:);
    S = sqrt(sum(smooth.^2,2));
    S = S/sum(S);
    smooth = smooth .* [S,S];
    smooth = [0,0;smooth;0,0];
    
    %this calculates the left and right points
    [left,right] = genEdges(points,width);
    
    %this calculates the pull on each point
    pull = zeros(length(points),2);
    for j = 1:length(points)
        p1 = [pointsg(:,1)-left(j,1),pointsg(:,2)-left(j,2)];
        p2 = [pointsg(:,1)-right(j,1),pointsg(:,2)-right(j,2)];
        P1 = abs(sqrt(sum(p1.^2,2)));
        P2 = abs(sqrt(sum(p2.^2,2)));
        P1 = (1-(P1/max(P1))).*double(P1==min(P1)).*double(P1<max_distance);
        P2 = (1-(P2/max(P2))).*double(P2==min(P2)).*double(P2<max_distance);
        P = (P1+P2)/2;
        if(binary(round(points(j,1)),round(points(j,2))))
            pull(j,:) = [sum((p1(:,1)+p2(:,1)).*P),sum((p1(:,2)+p2(:,2)).*P)];
        else
            pull(j,:) = 10000*[sum((p1(:,1)+p2(:,1)).*P),sum((p1(:,2)+p2(:,2)).*P)];
        end
    end
    %check for division by zero errors
    pull(isnan(pull)) = 0;

    %pull
    %calculate new points
    points = points + alpha.*dist + beta.*smooth +gamma.*pull;
    
    %measure change from previous points
    C = sum(sqrt(sum((prev_points - points).^2,2)));
    if(C<1)
        %beep;
        %disp(i)
        %disp('slowing');
        break;
    end
    
    if(visual)
        clf('reset');
        imshow(blank,'border','tight');
        hold on;
        plot(pointsg(:,2),pointsg(:,1),'c*');
        plot(pointso(:,2),pointso(:,1),'g-o');
        plot(points(:,2),points(:,1),'r-o');
        plot(left(:,2),left(:,1),'k-o');
        plot(right(:,2),right(:,1),'b-o');
        pause(.1)
    end
end
[left,right] = genEdges(points,width);
end