%%Initialize
clear;
close all;
clc;
clf;

%%Read Input Image

Image = imread('Car.jpg');
f = figure(1);
ax = axes(f);
% subplot(2,2,1);
imshow(Image);
title(ax,'Input Image');

ImageGray = rgb2gray(Image);

T = 0.255;
LOW_THRESHOLD  =  max(1,T*(mean2(ImageGray)-std2(ImageGray)))/255;
HIGH_THRESHOLD =  min(254,T*(mean2(ImageGray)+std2(ImageGray)))/255;

ImageCanny = edge(ImageGray,'canny', T);


% subplot(2,2,2);
f = figure(2);
ax = axes(f);
imshow(ImageCanny,[]);
title(ax,'Canny Edge Detection With Built-In Function');
Image = double(ImageGray);

%%Gaussian Filter
 gaussianFilter = fspecial('gaussian');
 
%Convolution
convolutionImage = conv2(Image, gaussianFilter, 'same');

%Gradient calculation
gaussianFilterX = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
gaussianFilterY = [1, 2, 1; 0, 0, 0; -1, -2, -1];


ImageGradientX = conv2(convolutionImage, gaussianFilterX, 'same');
ImageGradientY = conv2(convolutionImage, gaussianFilterY, 'same');

direction = atan2(ImageGradientY, ImageGradientX)*180/pi;

pan = size(convolutionImage,1);
leb = size(convolutionImage,2);

for i=1:pan
    for j=1:leb
        if (direction(i,j)<0) 
            direction(i,j)=360+direction(i,j);
        end
    end
end

direction2=zeros(pan, leb);

for i = 1  : pan
    for j = 1 : leb
        if ((direction(i, j) >= 0 ) && (direction(i, j) < 22.5) || (direction(i, j) >= 157.5) && (direction(i, j) < 202.5) || (direction(i, j) >= 337.5) && (direction(i, j) <= 360))
           direction2(i, j) = 0;
        elseif ((direction(i, j) >= 22.5) && (direction(i, j) < 67.5) || (direction(i, j) >= 202.5) && (direction(i, j) < 247.5))
            direction2(i, j) = 45;
        elseif ((direction(i, j) >= 67.5 && direction(i, j) < 112.5) || (direction(i, j) >= 247.5 && direction(i, j) < 292.5))
            direction2(i, j) = 90;
        elseif ((direction(i, j) >= 112.5 && direction(i, j) < 157.5) || (direction(i, j) >= 292.5 && direction(i, j) < 337.5))
            direction2(i, j) = 135;
        end
    end
end

Gmag = sqrt(ImageGradientX.^2 + ImageGradientY.^2);
BW = zeros (pan, leb);

%Non-Maximum Supression
for i=2:pan-1
    for j=2:leb-1
        if (direction2(i,j)==0)
            BW(i,j) = (Gmag(i,j) == max([Gmag(i,j), Gmag(i,j+1), Gmag(i,j-1)]));
        elseif (direction2(i,j)==45)
            BW(i,j) = (Gmag(i,j) == max([Gmag(i,j), Gmag(i+1,j-1), Gmag(i-1,j+1)]));
        elseif (direction2(i,j)==90)
            BW(i,j) = (Gmag(i,j) == max([Gmag(i,j), Gmag(i+1,j), Gmag(i-1,j)]));
        elseif (direction2(i,j)==135)
            BW(i,j) = (Gmag(i,j) == max([Gmag(i,j), Gmag(i+1,j+1), Gmag(i-1,j-1)]));
        end
    end
end

BW = BW.*Gmag;

%Hysteresis Thresholding
LOW_THRESHOLD = LOW_THRESHOLD * max(max(BW));
HIGH_THRESHOLD = HIGH_THRESHOLD * max(max(BW));
T_res = zeros (pan, leb);
for i = 1  : pan
    for j = 1 : leb
        if (BW(i, j) < LOW_THRESHOLD)
            T_res(i, j) = 0;
        elseif (BW(i, j) > HIGH_THRESHOLD)
            T_res(i, j) = 1;
        %Using 8-connected components
        elseif ( BW(i+1,j)>HIGH_THRESHOLD || BW(i-1,j)>HIGH_THRESHOLD || BW(i,j+1)>HIGH_THRESHOLD || BW(i,j-1)>HIGH_THRESHOLD || BW(i-1, j-1)>HIGH_THRESHOLD || BW(i-1, j+1)>HIGH_THRESHOLD || BW(i+1, j+1)>HIGH_THRESHOLD || BW(i+1, j-1)>HIGH_THRESHOLD)
            T_res(i,j) = 1;
        end
    end
end

edge_final = uint8(T_res.*255);

%Show final edge detection result
% subplot(2,2,3);
f = figure(3);
ax = axes(f);
imshow(edge_final);
title(ax,'Canny Edge Detection Without Built-In Function');