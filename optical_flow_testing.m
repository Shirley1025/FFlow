function Im1_shift_new = optical_flow_testing(Im1,Im2,ux,uy)

Im1=double(Im1);
Im2=double(Im2);
uxI=double(ux);  
uyI=double(uy);
[m1,n1]=size(Im1);

% generate a shifted image from Im1 based on the velocity field that is rounded

% Rounding to obtain an offset image
Im1_shift0=Im2;
for i=1:m1
    for j=1:n1
        i_shift=i+round(uyI(i,j));
        j_shift=j+round(uxI(i,j));
        if (i_shift<=m1) && (i_shift >=1) && (j_shift<=n1) && (j_shift>=1)
            Im1_shift0(i_shift,j_shift)=Im1(i,j);           
        else
            Im1_shift0(i,j)=Im1(i,j);
        end
    end
end

% Decimal part of velocity field
Im3=Im1_shift0;
Im1_shift1=Im3;
duxI=uxI-round(uxI);
duyI=uyI-round(uyI);

% Applying Gaussian filter H1 for smoothing
mask_size=10;
std=0.6*mask_size;
H1=fspecial('gaussian',mask_size,std);
duxI=imfilter(duxI,H1);
duyI=imfilter(duyI,H1);

dx=1;
dy=1;
dt=1;

for i=1:(m1-1)
    for j=1:(n1-1)
          term1(i,j)=(Im3(i,j+dx)*duxI(i,j+dx)-Im3(i,j)*duxI(i,j))/dx;  
          term2(i,j)=(Im3(i+dy,j)*duyI(i+dy,j)-Im3(i,j)*duyI(i,j))/dy;   
          Im1_shift1(i,j)=Im3(i,j)-(term1(i,j)+term2(i,j))*dt;   
    end
end

Im1_shift_new=uint8(Im1_shift1);



