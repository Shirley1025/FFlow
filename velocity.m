function [V,Vx,Vy,E] = velocity(Ux,Uy,frame_rate,Diff)
sizex=size(Ux, 1);
sizey=size(Ux, 2);
V=zeros(sizex,sizey);
Vx=zeros(sizex,sizey);
Vy=zeros(sizex,sizey);
E=zeros(sizex,sizey);
U=sqrt(Ux.^2+Uy.^2);
for i = 1:sizex  
    for j = 1:sizey
        ux=Ux(i,j,:);uy=Uy(i,j,:);u=U(i,j,:);e=Diff(i,j,:);
        Vx(i,j)=mean(abs(ux(ux~=0)))*frame_rate;
        Vy(i,j)=mean(abs(uy(uy~=0)))*frame_rate;
        V(i,j)=mean(u(u~=0))*frame_rate; 
        E(i,j)=mean(e(e~=0)); 
%        V(i,j)=roundn(mean(sort_data(abs(U(i,j,:)),0.4,1))*frame_rate,-3); 
    end 
end

