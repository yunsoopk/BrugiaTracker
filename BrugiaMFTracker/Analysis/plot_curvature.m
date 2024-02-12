function [ figure1 ] = plot_curvature( plotting_data, end_effects, tmin, tmax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Handles setting up the tmin and tmax values

if(~tmin)
    tmin = 1;
end

if(~tmax)
   tmax = size(plotting_data,1);
end
figure1 = figure(1);
ax = gca;
color_array = [0 0 1;0 1 1; 1 1 0; 1 0 0];
c = new_colormap(color_array,10);
colormap(c);
if(size(plotting_data,2)>1)
    if(end_effects)
        pcolor(fix(size(plotting_data,2)*.1):size(plotting_data,2)-...
            fix(size(plotting_data,2)*.1),tmin:tmax,...
            plotting_data(tmin:tmax,fix(size(plotting_data,2)*.1):...
            size(plotting_data,2)-fix(size(plotting_data,2)*.1)));
        xlabel('Skeletal body coordinate, 10 - 90%');
    else
        hand = pcolor(tmin:tmax,1:size(plotting_data,2),transpose(plotting_data(tmin:tmax,:)));
        shading interp;
        set(hand, 'linestyle', 'none');
        set(ax,'Ydir','reverse');
        xlabel('Time (seconds)');
    end
    ylabel('Body coordinate, I/L');
    h = colorbar;
    set(h, 'ylim', [min(plotting_data(:)) max(plotting_data(:))])
else
    plot(transpose(plotting_data(tmin:tmax,:)));
    
    ylabel('Curvature value');
    
% saveas(figure1,'plot.png');
end

