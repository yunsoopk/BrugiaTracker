function [] = save_plot(file_name, x, data, y_label, resolution, scatter)
%SAVE_PLOT Saves plots
%   Detailed explanation goes here
if ~scatter
    plot(x,data);
else
    plot(x, data, 'o');
end
title(['Mean : ' num2str(mean(data)) ...
        ' std\_dev : ' num2str(std(data))]);
xlabel('Time');
ylabel(y_label);
xlim([1 size(x,1)]);
set(gca, 'box', 'off');

% Save plot
print(file_name, '-dpng', ...
    ['-r' num2str(resolution)]);
end

