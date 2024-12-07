clear all;
close all;

%% Read the two images for which the optical flow needs to be calculated
Image1=imread('E:\LXL\LK-flow\consecutive_frame\6\17.png');
Image2=imread('E:\LXL\LK-flow\consecutive_frame\6\20.png');

lambda_1=20;  % the Horn_schunck estimator for initial field
lambda_2=2000; % the Liu-Shen estimator for refined estimation

no_iteration=1; 
scale_im=0.5; % scaling factor
size_average=20; 
size_filter=4;

%% 确定肠区轮廓及轮廓内所有像素点信息
threshold = 60; % 使用阈值进行二值化处理
[filtered_contours,x,y,x1,x2,y1,y2] = find_contours(Image1,threshold);
% 创建一个掩模
mask = poly2mask(x, y, size(Image1, 1), size(Image1, 2));
% 获取轮廓内部所有点的位置信息
[inner_x, inner_y] = find(mask);
p0 = cat(3, inner_y, inner_x);

%% Determination of the calculation area
I1=double(Image1(y1:y2,x1:x2)); 
I2=double(Image2(y1:y2,x1:x2));

I1_original=I1;
I2_original=I2;

%% correcting the global and local intensity change in images 
[m1,n1]=size(I1);
window_shifting=[1;n1;1;m1]; % [x1,x2,y1,y2] deines a rectangular window for global correction
[I1,I2]=correction_illumination(I1,I2,window_shifting,size_average);


%% pre-processing for reducing random noise,  
%% and downsampling images if displacements are large 
[I1,I2] = pre_processing_a(I1,I2,scale_im,size_filter);
I_region1=I1;
I_region2=I2;

%% initial optical flow calculation for a coarse-grained velocity field   liu-shen方法（hs方法作为初始值）计算粗粒度光流场
%% (ux0,uy0)

[ux0,uy0,vor,ux_horn,uy_horn,error1]=OpticalFlowPhysics_fun(I_region1,I_region2,lambda_1,lambda_2);
%% LK方法做初值
% [ux0,uy0,vor,ux_horn,uy_horn,error1]=lk_OpticalFlowPhysics_fun(I_region1,I_region2,lambda_1,lambda_2,p0,x1,x2,y1,y2);
% ux is the velocity (pixels/unit time) in the image x-coordinate (from the left-up corner to right)
% uy is the velocity (pixels/unit time) in the image y-coordinate (from the left-up corner to bottom)


%% generate the shifted image from Im1 based on the initial coarse-grained velocity field (ux0, uy0),  基于（ux0，uy0）从Im1生成偏移图像
%% and then calculate velocity difference for iterative correction  计算速度差进行迭代矫正
Im1=uint8(I1_original);
Im2=uint8(I2_original);

ux_corr=ux0;
uy_corr=uy0;

%% estimate the displacement vector and make correction in iterations                    估计位移矢量并在迭代中进行校正
k=1;
while k<=no_iteration
    [Im1_shift,uxI,uyI]=shift_image_fun_refine_1(ux_corr,uy_corr,Im1,Im2);
    
    I1=double(Im1_shift);
    I2=double(Im2);
    
    size_filter_1=2;
    [I1,I2] = pre_processing_a(I1,I2,1,size_filter_1);
    
    % calculation of correction of the optical flow 
    [dux,duy,vor,dux_horn,duy_horn,error2]=OpticalFlowPhysics_fun(I1,I2,lambda_1,lambda_2);
    %% LK方法做迭代
%     [dux,duy,vor,dux_horn,duy_horn,error2]=lk_OpticalFlowPhysics_fun(I1,I2,lambda_1,lambda_2,p0,x1,x2,y1,y2);

    % refined optical flow
    ux_corr=uxI+dux;
    uy_corr=uyI+duy;
    
    k=k+1;
end

%% refined velocity field  精细速度场
ux=ux_corr;    %%%%%
uy=uy_corr;    %%%%%

%% test velocity field  验证光流场
Im1=uint8(I1_original);
Im2=uint8(I2_original);
Im1_shift_new = optical_flow_testing(Im1,Im2,ux,uy);

%% show the images and processed results
%% plot the images, velocity vector, and streamlines in the initail and
%% refined estimations
%% Plot  refined  velocity vector field  绘制精细速度矢量场

figure(1);
gx=50; offset=1;
h = vis_flow (ux, uy, gx, offset, 2, 'm');
set(h, 'Color', 'black');
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
set(gca,'YDir','reverse');
% 将图形背景设为白色
set(gcf, 'Color', 'white');
title('Refined Velocity Field');
t = title(gca, 'Refined Velocity Field'); % 获取标题句柄
set(t, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');


%% plot the fields of velocity magnitude, vorticity and the second invariant Q
% plot velocity magnitude field
% calculate the velocity magnitude
u_mag=(ux.^2+uy.^2).^0.5;
u_max=max(max(u_mag));
u_mag=u_mag/u_max;
% calculate vorticity
vor=vorticity(ux, uy);
vor_max=max(max(abs(vor)));
vor=vor/vor_max;
% calculate the 2nd invariant
Q=invariant2_factor(ux, uy, 1, 1);

load('ColorMap_me_purple_1.mat', 'ColorMap_me')
figure(2);
colormap(ColorMap_me)
ulims=[0, 1];
imagesc(u_mag,ulims);
xlabel('x (pixels)');
ylabel('y (pixels)');
set(gca,'LineWidth',1);
% 隐藏刻度线
set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off'); % 关闭小刻度线
set(gca, 'TickLength', [0 0]); % 将刻度线长度设置为零
axis image;
set(gca,'YDir','reverse');
% 将图形背景设为白色
set(gcf, 'Color', 'white');
title('Velocity Magnitude Field');
t = title(gca, 'Velocity Magnitude Field of the FF'); % 获取标题句柄
set(t, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.5,209.5]);ylim([0.5,106.5]); % 将 z 轴的显示范围设置为 0-25
c=colorbar;
set(c,'tickdir','out')  % 朝外
% 设置颜色条的线宽
set(c, 'LineWidth', 1); % 将线宽设置为 2 点
hold on;
% plot streamlines
% figure(2);
[m,n]=size(ux);
[x,y]=meshgrid(1:n,1:m);
dn=10;
dm=10;
[sx,sy]=meshgrid(1:dn:n,1:dm:m);
h=streamslice(x, y, ux, uy, 4);
set(h, 'Color', 'black');
% 将图形背景设为白色
set(gcf, 'Color', 'white');
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
hold off;

% plot Vorticity field
load('ColorMap_me_purple_6.mat', 'ColorMap_me')
% load('ColorMap_me_gray.mat', 'ColorMap_me')
figure(3);
colormap(ColorMap_me)
% vlims=[0, 255];
% imagesc(I1_original,vlims);
imagesc(I1_original);
set(gca,'LineWidth',1);
xlabel('x (pixels)');
ylabel('y (pixels)');
% 隐藏刻度线
set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off'); % 关闭小刻度线
set(gca, 'TickLength', [0 0]); % 将刻度线长度设置为零
axis image;
set(gca,'YDir','reverse');
% 将图形背景设为白色
set(gcf, 'Color', 'white');
title('Vorticity Field');
t = title(gca, 'Velocity Vectors Field of the FF'); % 获取标题句柄
set(t, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.5,209.5]);ylim([0.5,106.5]); % 将 z 轴的显示范围设置为 0-25
c=colorbar;
set(c,'tickdir','out')  % 朝外
% 设置颜色条的线宽
set(c, 'LineWidth', 1); % 将线宽设置为 2 点
hold on;
% Plot  refined  velocity vector field
% figure(3);
gx=50; offset=1;
h = vis_flow (ux, uy, gx, offset, 3, 'm');
set(h, 'Color', 'black');
xlabel('x (pixels)');
ylabel('y (pixels)');
set(gca,'YDir','reverse');
% 将图形背景设为白色
set(gcf, 'Color', 'white');
hold off;
