function [ left, right ] = genEdges( points,width )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
one = diff(points,1,1);
two = points(2:end,:)-points(1:end-1,:);
t = one(1:end-1,:)+two(2:end,:);
a = sqrt(sum(t.^2,2));
t(:,1) = t(:,1)./a;
t(:,2) = t(:,2)./a;
t = [0,0;t;0,0];
left = [points(:,1)-t(:,2).*width,points(:,2)+t(:,1).*width];
right = [points(:,1)+t(:,2).*width,points(:,2)-t(:,1).*width];
end

