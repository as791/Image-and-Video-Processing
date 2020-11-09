%NAME: Aryaman Sinha
%INST: IIT, Bhubaneswar
%DATE: 09.11.2020
%CATEGORY: BTech
%BRANCH: Electrical Engineering
%Roll Number: 17EE01016

% Image and Video Tutorial 06
clc;clear;close all;
%%
% Q1. Implement affine transformation on the image.
img = imread('sample images/cameraman.tif');
img = double(img)./255;
n = size(img,1);m = size(img,2);
new_img1 = affine_transformation(img,90,1,1,0,-m);
new_img2 = affine_transformation(img,0,0.5,0.5,0,0);
new_img3 = affine_transformation(img,0,1,1,n/2,m/2);
subplot(2,2,1);imshow(img);title('sample img');
subplot(2,2,2);imshow(new_img1);title('90 ACW rotation');
subplot(2,2,3);imshow(new_img2);title('1/2 Scaled version with same origin');
subplot(2,2,4);imshow(new_img3);title('Origin translated to center');

function new_img = affine_transformation(img,theta,sx,sy,tx,ty)
    n = size(img,1);m = size(img,2);
    rot_mat = [[cosd(theta); sind(theta); 0],[-sind(theta); cosd(theta) ;0],[0;0;1]];
    scale_mat = [[sx ;0;0],[0;sy;0],[0;0;1]];
    trans_mat = [[1;0;0],[0;1;0],[tx;ty;1]];
    new_img = zeros(size(img));
    for x=1:n
        for y=1:m
            A = rot_mat*scale_mat*trans_mat;
            cord_mat = inv(A)*[x;y;1];
            pix = bilinear_interpolation(img,cord_mat(1),cord_mat(2));
            new_img(x,y) = pix;
        end
    end
end
function pix = bilinear_interpolation(img,x,y)
    n = size(img,1);m=size(img,2);
    x1 = floor(x); y1 = floor(y);
    x2 = ceil(x); y2 = ceil(y); 
    if(x1<=0)
        x1=1;
    end
    if(y1<=0)
        y1=1;
    end
    if(x1>=n)
        x1=n-1;
    end
    if(y1>=m)
        y1=m-1;
    end
    if(y2<=0)
        y2=1;
    end
    if(x2<=0)
        x2=1;
    end
    if(y2>m)
        y2=m;
    end
    if(x2>n)
        x2=n;
    end
    if(y2==y1&&x2~=x1)
         pix = ((x2-x)/(x2-x1))*img(x1,y1)+((x-x1)/(x2-x1))*img(x2,y1);
    elseif(x2==x1&&y2~=y1)
        pix = ((y2-y)/(y2-y1))*img(x1,y1)+((y-y1)/(y2-y1))*img(x1,y2);
    elseif(x2==x1&&y2==y1)
        pix = img(x1,y1);
    else
        pix_h1 = ((y2-y)/(y2-y1))*img(x1,y1)+((y-y1)/(y2-y1))*img(x1,y2);
        pix_h2 = ((y2-y)/(y2-y1))*img(x2,y1)+((y-y1)/(y2-y1))*img(x2,y2);
        pix = ((x2-x)/(x2-x1))*pix_h1+((x-x1)/(x2-x1))*pix_h2;
    end
end
%%
% *Discussion*: Here we implement the A affine transformation matrix to
% transform image geometrically according to our given parameters here we
% do for 90 anticlock wise accordingly origin translate similarly goes
% normal scaled and simple center translated versions. To appropriately
% copy the image pixel to new image we use inverse mapping method i.e. we
% find the (x,y) in original image for (x',y') in output image so that we
% get the (x,y) in input image accordingly we use simple bilinear
% interpolation to some extant to get those corresponding pixle value for
% output image. Hence, we get some decent results.