function [center, radius] = get_circular_well(img, radius_range, tolerance, ...
                                    manual, batch)
%get_circular_mask Finds the well plate circle and outputs the mask
%   Uses Hough Circle algorithm to find circles and creates a mask for the
%   inside of the well plate.

    % Image needs to be grayscale
    radii = [];
    if size(img,3) == 3
        img = rgb2gray(img);
    end
    if nargin ~= 5
        throw(MException('MATLAB:notEnoughInputs',...
                        'Not enough input arguments.'));
    end

    if manual == 0
        [centers, radii, ~] = imfindcircles(img,radius_range, ...
                                'ObjectPolarity','dark');
    end
    if (isempty(radii) && batch)
        mask = zeros(size(img));
        return
    end
    if (manual || isempty(radii))
        RGB = insertText(img,[10 10],...
            'Select the well and DBL-CLICK to exit','FontSize',18,...
            'BoxOpacity',0.4,'TextColor','white');
        fig_h = figure;
        imshow(RGB);
        h = imellipse(gca,[size(img,2)/4 size(img,1)/4 ...
            sum(radius_range) sum(radius_range)]);
        wait(h);
        pos = h.getPosition;
        close(fig_h);
        radii = (max(pos(3), pos(4))/2);
        centers = [pos(1) + radii, pos(2) + radii];
    end
    % Find maximum circle (Assuming a single well plate in the image)
    [radius, index] = max(radii);
    radius = radius - tolerance;
    center = centers(index, :);
end

