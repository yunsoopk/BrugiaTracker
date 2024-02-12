function [ side_1, side_2 ] = splitBound( bound, head, tail )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~,head_index] = min(sqrt(sum((bound(:,1)-head(:,1)).^2+(bound(:,2)-head(:,2)).^2,2)));
%head_index = find(bound(:,1)==head(1)&bound(:,2)==head(2));
%tail_index = find(bound(:,1)==tail(1)&bound(:,2)==tail(2));
[~,tail_index] = min(sqrt(sum((bound(:,1)-tail(:,1)).^2+(bound(:,2)-tail(:,2)).^2,2)));

if(head_index==tail_index)
    head_index = head_index+1;
    if(head_index>length(bound))
        head_index = 1;
    end
end

if(head_index<tail_index);
    side_1 = bound(head_index:tail_index,:);
    side_2 = [bound(tail_index:end,:);bound(1:head_index,:)];
else
    side_1 = bound(tail_index:head_index,:);
    side_2 = [bound(head_index:end,:);bound(1:tail_index,:)];
end

end

