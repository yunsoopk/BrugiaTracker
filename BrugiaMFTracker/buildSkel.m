function [ skel,widths ] = buildSkel( side_1, side_2, num_points )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
skel = zeros(num_points,2);
if(length(side_1)<length(side_2))
    t = side_1;
    side_1 = side_2;
    side_2 = t;
end

widths = skel(:,1);
skel(1,:) = side_1(1,:);

dist = sqrt(sum(diff(side_1,1).^2,2));
for i = 2:length(dist)
    dist(i) = dist(i-1)+dist(i);
end
interval = dist(end)/(num_points-1);

%  plot(side_1(:,2),side_1(:,1))
%  hold on;
%  plot(side_2(:,2),side_2(:,1))

for i = 2:num_points
    side1_dist = abs(dist-(i-1)*interval);
    [~,side1_idx] = min(side1_dist);
    %plot(side_1(side1_idx,2),side_1(side1_idx,1),'g*');  
    dist2 = sqrt(((side_2(:,1)-side_1(side1_idx,1)).^2 ...
        +(side_2(:,2)-side_1(side1_idx,2)).^2));
    [widths(i),side2_idx] = min(dist2);
    %plot(side_2(side2_idx,2),side_2(side2_idx,1),'b*');
    skel(i,:) = (side_1(side1_idx,:)+side_2(side2_idx,:))/2;
    %if(skel(i,1)==side_1(side1_idx,1)&&skel(i,2)==side_1(side1_idx,2))
    if(i>10&&i<num_points-10&&widths(i)<5)
        disp('error');
    end
end

X = transpose(spline(1:num_points,transpose(skel(:,2)),1:.2:num_points));
Y = transpose(spline(1:num_points,transpose(skel(:,1)),1:.2:num_points));

dist = [0;sqrt(sum(diff([Y,X]).^2,2))];
for i = 2:length(dist)
    dist(i) = dist(i-1)+dist(i);
end
interval = dist(end)/(num_points-1);

for i = 2:num_points-1
    side1_dist = abs(dist-(i-1)*interval);
    [~,side1_idx] = min(side1_dist);
    skel(i,:) = [Y(side1_idx),X(side1_idx)];
end


end