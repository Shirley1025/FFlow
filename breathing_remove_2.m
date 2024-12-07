function [b_remove,ymid_diff] = breathing_remove_2(filename,breath_threshold,contour_threshold)
obj = VideoReader(filename);
num=obj.NumFrames;
ymid=zeros(num, 1);
ymid_diff=zeros(num, 1);
b_remove=ones(num-1, 1);

for i=1:num  
    Image1=read(obj,i);
    [filtered_contours,x,y,x1,x2,y1,y2] = find_contours(Image1,contour_threshold);   
    ymid(i)=fix((y2 + y1)/2);
end

%% Calculate the position offset
for i=1:num-1 
    ymid_diff(i)=ymid(i+1)-ymid(i);
end

%% Marks frames with position offsets above a threshold as 0
for i=1:num-1 
    if (abs(ymid_diff(i))>=breath_threshold)
        b_remove(i)=0;
    end
end

