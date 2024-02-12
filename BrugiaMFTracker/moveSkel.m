function [ left, right ] = moveSkel( skel, widths, borders, B, visual )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[left,right] = genEdges(skel,widths/2);
P = [];
for i = 1:length(borders)
    P = [P;borders{i}];
end

if(1)
    [rows,cols] = find(B);
    xmin = min(cols);
    ymin = min(rows);
    xmax = max(cols);
    ymax = max(rows);
end

for i = 3:length(left)-2
    distl = sqrt(sum((P(:,2)-left(i,2)).^2+(P(:,1)-left(i,1)).^2,2));
    %distl(distl>widths(i)) = inf;
    distr = sqrt(sum((P(:,2)-right(i,2)).^2+(P(:,1)-right(i,1)).^2,2));
    %distr(distr>widths(i)) = inf;
    dists = abs(sqrt(sum((P(:,2)-skel(i,2)).^2+(P(:,1)-skel(i,1)).^2,2))-widths(i)/2);
    C = point_to_line(P,left(i,:),right(i,:));
    
    %K = menger_curve(skel);
    %[kval,ind] = max(K);
    if(visual)
        skel_temp = (left+right)/2;
        clf('reset');
        subplot(1,2,1);
        imshow(B(ymin:ymax,xmin:xmax));
        hold on;
        plot(P(:,2)-xmin,P(:,1)-ymin,'c-')
        plot(left(i,2)-xmin,left(i,1)-ymin,'ro');
        plot(right(i,2)-xmin,right(i,1)-ymin,'go');
        plot(skel(i,2)-xmin,skel(i,1)-ymin,'bo');
        plot(skel_temp(:,2)-xmin,skel_temp(:,1)-ymin,'y-');
        subplot(1,2,2);
        plot(distl+abs(distr-widths(i))+dists+C,'b-');
        plot(abs(distl-widths(i))+distr+dists+C,'r-');
        set(gcf, 'Position', get(0,'Screensize'));
        pause(.03)
    end
    %}
    costl = distl+abs(distr-widths(i))+dists+C;
    costr = abs(distl-widths(i))+distr+dists+C;
    [leftv,lidx] = min(costl);
    [rightv,ridx] = min(costr);

    if(leftv<rightv)% && leftv<widths(i))
        temp = round((P(lidx,:)+right(i,:))/2);
        %if(B(temp(1),temp(2)) || leftv<widths(i)/2)
            movement = P(lidx,:)-left(i,:);
            left(i,:) = P(lidx,:);
            right(i,:) = right(i,:)+movement;
            %{
        else
            disp('flip');
            disp(i);
            movement = P(lidx,:)-right(i,:);
            right(i,:) = P(lidx,:);
            left(i,:) = left(i,:)-movement;
        end
            %}
    else%if(rightv < widths(i))
        temp = round((P(ridx,:)+left(i,:))/2);
        %if(B(temp(1),temp(2)) || rightv<widths(i)/2)
        %if(rightv<widths(i)/2)
            movement = P(ridx,:)-right(i,:);
            right(i,:) = P(ridx,:);
            left(i,:) = left(i,:)+movement;  
            %{
        else
            disp('flip');
            disp(i);
            movement = P(ridx,:)-left(i,:);
            left(i,:) = P(ridx,:);
            right(i,:) = right(i,:)-movement;
        end
            %}
    end
end



end


