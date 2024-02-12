function [ K ] = menger_curve_bound( data )
%Menger Curvature Formula
%C(t1,t2,t3)=1/R=4A/(|t1-t2||t2-t3||t3-t1|)

K=zeros(1,length(data));
for i = 1:length(data)
    plus_one = i+1;
    if(plus_one>length(data))
        plus_one = plus_one-length(data);
    end
    plus_two = i+2;
    if(plus_two>length(data))
        plus_two = plus_two-length(data);
    end
    x1 = data(i,1);
    x2 = data(plus_one,1);
    x3 = data(plus_two,1);
    y1 = data(i,2);
    y2 = data(plus_one,2);
    y3 = data(plus_two,2);
    
 K(i) = 2*((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./ ...
  sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2));

%if any of the points are coincident curvature is defined as zero
if(((x2==x1)&&(y2==y1))||((x3==x1)&&(y3==y1))||((x2==x3)||(y2==y3)))
    K(i) = 0;
end

end

