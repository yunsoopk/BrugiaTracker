function [ new_colors ] = new_colormap( c_array, gradient)
% new_colormap: This function generates a new custom colormap to be used 
%   in plotting functions.
%
%   c_array: This is the array which you will store the colors you want to 
%           range from. The highest value will be assigned to the last
%           color.

%   gradient: This is the number of intermediate colors you want when going
%           from color to color.

num_colors = size( c_array, 1);

for i = 1 : num_colors - 1;
	diff( i, :) = ( c_array( i, :) - c_array( i+1, :)) / gradient;
end

j = 1;
index = 1;
new_colors( 1, :) = c_array( 1, :);

for i = 1 : ( gradient * (num_colors-1))
	new_colors( i+1, :) = new_colors( i, :) - diff( j, :);
    
    if( ~( i - gradient * index))
		j = j+1;
		index = index + 1;
    end
end

%This corrects any rounding errors that may happen in the computer
new_colors = abs(new_colors);

for i = 1:size(new_colors,1)
    for j = 1:size(new_colors,2)
        if(new_colors(i,j) > 1)
           new_colors(i,j) = 1; 
        end
        
    end
end

end


