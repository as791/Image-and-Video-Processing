%NAME: Aryaman Sinha
%INST: IIT, Bhubaneswar
%DATE: 01.11.2020
%CATEGORY: BTech
%BRANCH: Electrical Engineering
%Roll Number: 17EE01016

% Image and Video Tutorial 04
clc;clear;close all;
%%
% Q1. Assume degradation model derived for motion blurring, apply it to any
% good quality image to create the distorted image. Then apply full inverse
% filtering, radially limited inverse filtering and Wiener filtering to
% restore the original image.
img = imread('./sample images/cameraman.tif');
img = double(img)./255;
subplot(2,3,1);imshow(img);title('original image');
img_blur = motion_blur(img,1,0.05,0.01);
subplot(2,3,2);imshow(img_blur);title('motion blurred image');
img_inv_filt = inverse_filt(img_blur,1,0.05,0.01);
subplot(2,3,3);imshow(img_inv_filt);title('inverse filtered image');
img_rad_inv_filt = radial_inv_filt(img_blur,1,0.05,0.01,70);
subplot(2,3,4);imshow(img_rad_inv_filt);title('radial inverse filterd image');
img_weiner_filt = weiner_filt(img_blur,1,0.05,0.01,0.0067); 
subplot(2,3,5);imshow(img_weiner_filt);title('weiner filtered image');
%%
% Q2. While capturing any face image, ask the person to move, the resultant
% image will be distorted, then try the above degaration model to restore
% the face image.
img = imread('./sample images/face_image.jpg');
img = double(img)./255;
img_blur = img(:,:,1);
subplot(1,2,1);imshow(img_blur);title('motion blurred image');
img_weiner_filt = weiner_filt(img_blur,1,0,0.000000000000001,0.0000000001); 
subplot(1,2,2);imshow(img_weiner_filt);title('weiner filtered image');
%%
function gp = motion_blur(img,T,a,b)
    n = size(img,1); m = size(img,2);
    u = repelem(linspace(0,n-1,n)',1,m);v = repelem(linspace(0,m-1,m)',1,n)';
    F_uv = fft2(img);
    % derived degradation model of motion blurring 
    H_uv = T.*sinc(pi*(a*u+b*v)).*exp(-1i*pi*(a*u+b*v));
    G_uv = H_uv.*F_uv;
    gp = real(ifft2(G_uv));
    gp = (gp-min(gp(:)))./(max(gp(:))-min(gp(:)));
end

function fp = inverse_filt(img,T,a,b)
    n = size(img,1); m = size(img,2);
    u = repelem(linspace(0,n-1,n)',1,m);v = repelem(linspace(0,m-1,m)',1,n)';
    G_uv = fft2(img);
    % derived degradation model of motion blurring 
    H_uv = T.*sinc(pi*(a*u+b*v)).*exp(-1i*pi*(a*u+b*v));
    F_uv = G_uv./H_uv;
    fp = real(ifft2(F_uv));
    fp = (fp-min(fp(:)))./(max(fp(:))-min(fp(:)));
end

function fp = radial_inv_filt(img,T,a,b,D_thresh)
    n = size(img,1); m = size(img,2);
    u = repelem(linspace(0,n-1,n)',1,m);v = repelem(linspace(0,m-1,m)',1,n)';
    G_uv = fft2(img);
    % derived degradation model of motion blurring 
    H_uv = T.*sinc(pi*(a*u+b*v)).*exp(-1i*pi*(a*u+b*v));
    D_uv = sqrt((u-n/2).^2+(v-m/2).^2);
    lpf_H_uv = 1./(1+(D_uv./D_thresh).^20);
    F_uv = G_uv./H_uv;
    F_uv = lpf_H_uv.*F_uv;
    fp = real(ifft2(F_uv));
    fp = (fp-min(fp(:)))./(max(fp(:))-min(fp(:)));
end

function fp = weiner_filt(img,T,a,b,k)
    n = size(img,1); m = size(img,2);
    u = repelem(linspace(0,n-1,n)',1,m);v = repelem(linspace(0,m-1,m)',1,n)';
    G_uv = fft2(img);
    % derived degradation model of motion blurring 
    H_uv = T.*sinc(pi*(a*u+b*v)).*exp(-1i*pi*(a*u+b*v));
    F_uv = G_uv.*(((abs(H_uv).^2)./(abs(H_uv).^2+k))./H_uv);
    fp = real(ifft2(F_uv));
    fp = (fp-min(fp(:)))./(max(fp(:))-min(fp(:)));
end

%%
% *Discussion*: Here we implement the motion blurring using the derived
% degradation model in the class. And the filter that motion blurred image
% using inverse filtering full and radial both we see that they don't get
% the resultant output that we want as saturation happens in inverse
% filtering as H_uv might have low value which may alter the resultant
% output this is appropriately handeled by wiener filtering which give good
% results than other two filtering methods. Here we also try to deblur
% the blur image of face using same degradation model, however its very
% specific to get the appropriate paramater setting to get desired result.