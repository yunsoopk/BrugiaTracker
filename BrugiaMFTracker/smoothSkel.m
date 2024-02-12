function [ skel ] = smoothSkel( skel, MC )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
K = abs(menger_curve(skel));
iterations = 1;
while(max(K(:))>MC)% && iterations<100)
    x = 2;
    [~,I] = max(K(:));
    skel(I,:) = sum(skel(I-1:I+1,:).*[1,1;x,x;1,1],1)/(x+2);
    K = abs(menger_curve(skel));
    iterations = iterations +1;
end

end

