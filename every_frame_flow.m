function [ux,uy,I1_original,I2_original] = every_frame_flow(Im1,Im2,x1,x2,y1,y2)
    I1=double(Im1(y1:y2,x1:x2)); 
    I2=double(Im2(y1:y2,x1:x2));
      
    I1_original=I1;
    I2_original=I2;
    
    lambda_1=20;  % the Horn_schunck estimator for initial field
    lambda_2=2000; % the Liu-Shen estimator for refined estimation
    size_average=20;
    size_filter=4;
    scale_im=0.5;
    no_iteration=1;
    
    % correcting the global and local intensity change in images 
    [m1,n1]=size(I1);
    window_shifting=[1;n1;1;m1]; % [x1,x2,y1,y2] deines a rectangular window for global correction
    [I1,I2]=correction_illumination(I1,I2,window_shifting,size_average);

    % pre-processing for reducing random noise,  
    % and downsampling images if displacements are large 
    [I1,I2] = pre_processing_a(I1,I2,scale_im,size_filter);

    I_region1=I1;
    I_region2=I2;
    
    % initial optical flow calculation for a coarse-grained velocity field  
    % (ux0,uy0)
    [ux0,uy0,vor,ux_horn,uy_horn,error1]=OpticalFlowPhysics_fun(I_region1,I_region2,lambda_1,lambda_2);

    % generate the shifted image from Im1 based on the initial coarse-grained velocity field (ux0, uy0), 
    % and then calculate velocity difference for iterative correction  
    Im1=uint8(I1_original);
    Im2=uint8(I2_original);

    ux_corr=ux0;
    uy_corr=uy0;
    
    % estimate the displacement vector and make correction in iterations       
    k=1;
    while k<=no_iteration
        [Im1_shift,uxI,uyI]=shift_image_fun_refine_1(ux_corr,uy_corr,Im1,Im2);

        I1=double(Im1_shift);
        I2=double(Im2);

        size_filter_1=2;
        [I1,I2] = pre_processing_a(I1,I2,1,size_filter_1);

        % calculation of correction of the optical flow 
        [dux,duy,vor,dux_horn,duy_horn,error2]=OpticalFlowPhysics_fun(I1,I2,lambda_1,lambda_2);

        % refined optical flow
        ux_corr=uxI+dux;
        uy_corr=uyI+duy;

        k=k+1;
    end

    % refined velocity field  
    ux=ux_corr;    
    uy=uy_corr;    
end

