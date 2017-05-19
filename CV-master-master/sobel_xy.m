function [Fx,Fy] = sobel_xy(Image)

sobel_y=[1,2,1;0,0,0;-1,-2,-1];
sobel_x=sobel_y';
Fx=conv2(Image,sobel_x,'same');
Fy=conv2(Image,sobel_y,'same');


end






