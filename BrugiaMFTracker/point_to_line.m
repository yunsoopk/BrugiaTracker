function dist = point_to_line(x, a, b)
%x = [0,0]; %some point
%a = [1,2]; %segment points a,b
%b = [3,5];

slope = (a(1)-b(1))/(a(2)-b(2));
const = a(1)-slope*b(2);

slope2 = -1/slope;
const2 = x(:,1)-slope2*x(:,2);

intersection_x = (const2-const)./(slope-slope2);
intersection_y = slope2.*intersection_x+const2;

dist = sqrt((intersection_x-x(:,2)).^2+(intersection_y-x(:,1)).^2);
