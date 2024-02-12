function [ skel, interval ] = distributePoints( skel )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
num_points = length(skel);
%{
tmp = skel(1,:);
for i = 1:num_points-2
    dist = sqrt(sum(diff(skel(i:i+2,:),1,1).^2));
    if(dist(1)<dist(2))
        tmp = [tmp;skel(i:i+2)
end
%}
%[ skel ] = smoothSkel( skel, .3 );
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

