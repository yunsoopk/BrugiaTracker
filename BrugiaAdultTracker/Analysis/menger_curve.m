function [ K ] = menger_curve( data )
%Menger Curvature Formula
%C(t1,t2,t3)=1/R=4A/(|t1-t2||t2-t3||t3-t1|)

K=zeros(1,length(data));
for i = 2:length(data)-1

    x1 = data(i-1,1);
    x2 = data(i,1);
    x3 = data(i+1,1);
    y1 = data(i-1,2);
    y2 = data(i,2);
    y3 = data(i+1,2);
    
 K(i) = 2*((x1-x3).*(y2-y1) - (x1-x2).*(y3-y1)) ./ ...
  sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2));

%if any of the points are coincident curvature is defined as zero
if(((x2==x1)&&(y2==y1))||((x3==x1)&&(y3==y1))||((x2==x3)||(y2==y3)))
    K(i) = 0;
end

end

