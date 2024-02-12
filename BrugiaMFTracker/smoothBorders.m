function [ borders ] = smoothBorders( borders, smoothing_coefficient )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
for i = 1:length(borders)
    d = borders{i};
    for j = 1:smoothing_coefficient
        d = (d+[d(2:end,:);d(1,:)]+[d(end,:);d(1:end-1,:)])/3;
    end
    borders(i) = {d};
end

end

