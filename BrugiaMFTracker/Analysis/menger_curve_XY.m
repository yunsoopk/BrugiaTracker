function [ K ] = menger_curve_XY( X, Y )
%MENGER_CURVE_ARRAY calculates the curvature (change in tangent angle along
%the midline of a worm
%X is column vectors of x coordinates for all points
%Y is column vecotrs of y coordinates for all points

%change in x
d1 = diff(X,1,1);
%change in y
d2 = diff(Y,1,1);
%distance between points
s = sqrt(d1.^2+d2.^2);
%ignore first distance
s(1,:) = [];
t = atan(d2./d1);
dt = diff(t,1,1);
dt = dt+double(dt>pi/2).*-pi;
dt = dt+double(dt<-pi/2).*pi;
K = dt./s;
    
end

