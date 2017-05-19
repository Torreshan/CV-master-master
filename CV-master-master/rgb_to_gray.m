function [Gray_image] = rgb_to_gray(Image)
Imagesize=size(Image); 
if numel(Imagesize)>2
Image=double(Image);
Gray_image=0.299*Image(:,:,1)+0.587*Image(:,:,2)+0.114*Image(:,:,3);
Gray_image=Gray_image/255;
%disp(size(Gray_image));
else    
    disp('Please give a RGB image!');

% % Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
% das Bild bereits in Graustufen vorliegt, soll es direkt zur?kgegeben werden.

end
