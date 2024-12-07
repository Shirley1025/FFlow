close all;
clear all;


%% Specify the required path
folder = 'E:\LXL\LK-flow\video\1129';  % Specify the required path
filename = 'E:\LXL\LK-flow\video\1129\velocity-1129_1.xlsx';  % Specify the required path
frame_rate=6.92;  % Specify the required path

contour_threshold =60;  % Thresholds for outlining intestinal regions in mice
breath_threshold=3;  % Threshold for labeling breathing frames
df_threshold=8;  % Thresholding for smoothing velocity images

%% All video files to be calculated
fileList = dir(fullfile(folder, '*.avi')); %  Specify the suffix of the video file, such as avi or mp4.
numFiles = length(fileList);

%% Cyclically reads the videos in the folder and performs optical flow calculations for each video
for k = 1:numFiles
    file_number = k;
    path = fullfile(folder, fileList(k).name); 
    disp(path);
    
    %% Breathing frame rejection:
    %% Flag as invalid frames that are visibly jittering up and down due to the effects of breathing
    [b_remove,ymid_diff] = breathing_remove_2(path,breath_threshold,contour_threshold);

    %% Determine the contour of the mouse intestinal area 
    %% and the positional information of all pixel points within the contour
    obj=VideoReader(path);
    Image1=read(obj,1);
    Image2=read(obj,2);
    threshold =60; 
    [filtered_contours,x,y,x1,x2,y1,y2] = find_contours(Image1,threshold);
    mask = poly2mask(x, y, size(Image1, 1), size(Image1, 2));  % Creating a Mask
    [inner_x, inner_y] = find(mask);  % Get position information of all points inside the contour

    %% Drawing the outline of mouse intestinal area and rectangular calculation area
    plot_contours;

    %% Use the FF method to cycle through the optical flow every two frames
    num=obj.NumFrames;
    Ux=zeros(abs(y1-y2)+1,abs(x1-x2)+1,num-1);Uy=zeros(abs(y1-y2)+1,abs(x1-x2)+1,num-1);
    FF_Diff=zeros(abs(y1-y2)+1,abs(x1-x2)+1,num-1);

    if (num>100)
        cir_num=100;
    else
        cir_num=num;
    end     

    ymid_diff_new=ymid_diff(1:cir_num-1);

    error_im1_im2=zeros(cir_num-1, 1);
    error_ff=zeros(cir_num-1, 1);
    flow_mean=zeros(cir_num-1, 1);

    for i=1:cir_num-1  
        if(b_remove(i)~=0)
            Im1=read(obj,i);
            Im2=read(obj,i+1);
            [ux,uy,I1_original,I2_original] = every_frame_flow(Im1,Im2,x1,x2,y1,y2);

            u=sqrt(ux.^2+uy.^2);
            x=inner_x-y1+2;y=inner_y-x1+2;
            linear_index = sub2ind(size(u), x, y);
            z=u(linear_index);
            u_mean=mean(z,'all');
            flow_mean(i)=u_mean;

            % test velocity field 
            n=i;
            Im1=uint8(I1_original);Im2=uint8(I2_original);
            Im1_shift_new = optical_flow_testing(Im1,Im2,ux,uy);

            diff_img = abs(double(Im1) - double(Im2));
            x=inner_x-y1+2;y=inner_y-x1+2;
            linear_index = sub2ind(size(diff_img), x, y);
            z=diff_img(linear_index);
            part_error_avg_im1im2=mean(z,'all');       
            error_im1_im2(i)=part_error_avg_im1im2;

            ff_diff_img = abs(double(Im1_shift_new) - double(Im2));
            x=inner_x-y1+2;y=inner_y-x1+2;
            linear_index = sub2ind(size(ff_diff_img), x, y);
            z=ff_diff_img(linear_index);
            part_error_avg_af=mean(z,'all');
            error_ff(i)=part_error_avg_af;
            Ux(:,:,i)=ux;Uy(:,:,i)=uy;FF_Diff(:,:,i)=ff_diff_img;
        end
    end

    [V,Vx,Vy,E] = velocity(Ux,Uy,frame_rate,FF_Diff);
    
    %% Calculation of speed indicators
    linear_index = sub2ind(size(V), x, y);
    z=V(linear_index);
    e=E(linear_index);
    max_velocity = max(z, [], 'all');
    nonzero_elements = z(z ~= 0);
    min_velocity = min(nonzero_elements, [], 'all');
    mean_velocity = round(mean(nonzero_elements, 'all'),2);
    disp(mean_velocity);
    faster = z(z > 0.7*max_velocity);slower = z(z < min_velocity/0.7);
    faster_velocity = round(mean(faster, 'all'),2); slower_velocity = round(mean(slower, 'all'),2);
    disp(['平均速度:',num2str(mean_velocity)]);
    disp(['最大速度:',num2str(max_velocity)]);disp(['最小速度:', num2str(min_velocity)]);
    disp(['快区速度:',num2str(faster_velocity)]);disp(['慢区速度:', num2str(slower_velocity)]);

    %% Create a table and write speed values
    [num, txt, raw] = xlsread(filename);
    number = (file_number - 1) * 7+1;
    raw{1, number} = file_number;
    % 扩展 raw 数组
    numRows = size(x, 1);
    if numRows > size(raw, 1)
        raw(end+1:numRows, :) = {[]}; 
    end
    raw{1, number + 1} = 'x';
    raw(2:numRows+1, number + 1) = num2cell(x);
    raw{1, number + 2} = 'y';
    raw(2:numRows+1, number + 2) = num2cell(y);
    raw{1, number + 3} = 'v';
    raw(2:numRows+1, number + 3) = num2cell(z);
    raw{1, number + 4} = 'point-error';
    raw(2:numRows+1, number + 4) = num2cell(e);
    raw{2, number + 5} = 'avg';
    raw{3, number + 5} = 'max';
    raw{4, number + 5} = 'min';
    raw{5, number + 5} = 'faster';
    raw{6, number + 5} = 'slower';
    raw{2, number + 6} = mean_velocity;
    raw{3, number + 6} = max_velocity;
    raw{4, number + 6} = min_velocity;
    raw{5, number + 6} = faster_velocity;
    raw{6, number + 6} = slower_velocity;
    xlswrite(filename,raw);  

end