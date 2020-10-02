%NAME: Aryaman Sinha
%INST: IIT, Bhubaneswar
%DATE: 02.10.2020
%CATEGORY: BTech
%BRANCH: Electrical Engineering
%Roll Number: 17EE01016

% Image and Video Tutorial 03
clc;clear;close all;
%%
% Q1. Perform Laplace based sharpening in the frequency domain.
img = imread('./sample images/cameraman.tif');
img = double(img)/255;
subplot(1,2,1);imshow(img);title('original image');
n = size(img,1); m = size(img,2);
p=2*n;q=2*m;
img_pad = zeros(p,q);
img_pad(1:n,1:m) = img;
x = repelem(linspace(0,p-1,p)',1,q);y = repelem(linspace(0,q-1,q)',1,p)';
img_pad = img_pad.*((-1).^(x+y));
ft_img_pad = fft2(img_pad);
F_uv = abs(ft_img_pad);
u=x;v=y;
%laplacian operator in fourier domain 
H_uv = -(4*(pi^2)).*((u-p/2).^2+(v-q/2).^2);
G_uv = H_uv.*F_uv;
gp = real(ifft2(G_uv.*exp(1i*angle(ft_img_pad)))).*((-1).^(x+y));
g_xy = gp(1:n,1:m);
g_xy=g_xy./max(g_xy(:));
sharp_img = img-g_xy;
subplot(1,2,2);imshow(sharp_img);title('sharpened image');figure;
%%
% Q2. Perform Gaussian based Lowpass and Highpass filtering in a
% fingerprint image in frequency domain.
img = imread('./thumb_print.png');
% img = rgb2gray(img);
img=img(:,:,1);
img = double(img)/255;
subplot(2,3,1);imshow(img);title('Smuged thumb print');
n = size(img,1); m = size(img,2);
p=2*n;q=2*m;
img_pad = zeros(p,q);
img_pad(1:n,1:m) = img;
x = repelem(linspace(0,p-1,p)',1,q);y = repelem(linspace(0,q-1,q)',1,p)';
img_pad = img_pad.*((-1).^(x+y));
ft_img_pad = fft2(img_pad);
F_uv = abs(ft_img_pad);
u=x;v=y;
D_uv = sqrt((u-p/2).^2+(v-q/2).^2);
D_thresh = 50;
%Guassian
lpf_H_uv = exp(-(D_uv.^2)./(2*D_thresh*D_thresh));
hpf_H_uv = 1-lpf_H_uv;
%hpf
hpf_G_uv = hpf_H_uv.*F_uv;
gp = real(ifft2(hpf_G_uv.*exp(1i*angle(ft_img_pad)))).*((-1).^(x+y));
g_xy = gp(1:n,1:m);
hpf_img=g_xy./max(g_xy(:));
subplot(2,3,2);imshow(hpf_img);title('Gaussian HPF output');
%lpf
lpf_G_uv = lpf_H_uv.*F_uv;
gp = real(ifft2(lpf_G_uv.*exp(1i*angle(ft_img_pad)))).*((-1).^(x+y));
g_xy = gp(1:n,1:m);
lpf_img=g_xy./max(g_xy(:));
subplot(2,3,3);imshow(lpf_img);title('Gaussian LPF output');
%butterworth
hpf_H_uv = 1./(1+(D_thresh./D_uv).^8);
lpf_H_uv = 1./(1+(D_uv./D_thresh).^8);
%hpf
hpf_G_uv = hpf_H_uv.*F_uv;
gp = real(ifft2(hpf_G_uv.*exp(1i*angle(ft_img_pad)))).*((-1).^(x+y));
g_xy = gp(1:n,1:m);
hpf_img=g_xy./max(g_xy(:));
subplot(2,3,4);imshow(hpf_img);title('ButterWorth HPF output');
%lpf
lpf_G_uv = lpf_H_uv.*F_uv;
gp = real(ifft2(lpf_G_uv.*exp(1i*angle(ft_img_pad)))).*((-1).^(x+y));
g_xy = gp(1:n,1:m);
lpf_img=g_xy./max(g_xy(:));
subplot(2,3,5);imshow(lpf_img);title('ButterWorth LPF output');
out_thresh_img = logical(hpf_img>=0);
subplot(2,3,6);imshow(out_thresh_img);title('Thresholded HPF of image');figure;
%%
% Q3. Take a noisy fingerprint image, then preprocess the image using
% different morphological operators.
img = imread('./fingerprint.png');
% img = imbinarize(rgb2gray(img));
img = imbinarize(img(:,:,1));
subplot(2,3,1);imshow(img);title('noisy fingerprint');
n=3;w=floor(n/2);
element = ones(n,n);
out_img = zeros(size(img));
%erosion
for i=w+1:size(img,1)-w
    for j=w+1:size(img,2)-w
        set_a = img(i-w:i+w,j-w:j+w);
        out_img(i,j)=all(set_a(element>0),'all');
    end
end
subplot(2,3,2);imshow(out_img);title('Erosion output');
img = out_img;
%dilation 
for i=w+1:size(img,1)-w
    for j=w+1:size(img,2)-w
        set_a = img(i-w:i+w,j-w:j+w);
        out_img(i,j)=any(set_a(element>0),'all');
    end
end
subplot(2,3,3);imshow(out_img);title('Opening output');
img = out_img;
%dilation
for i=w+1:size(img,1)-w
    for j=w+1:size(img,2)-w
        set_a = img(i-w:i+w,j-w:j+w);
        out_img(i,j)=any(set_a(element>0),'all');
    end
end
subplot(2,3,4);imshow(out_img);title('Dilation of opening');
img = out_img;
%erosion
for i=w+1:size(img,1)-w
    for j=w+1:size(img,2)-w
        set_a = img(i-w:i+w,j-w:j+w);
        out_img(i,j)=all(set_a(element>0),'all');
    end
end
subplot(2,3,5);imshow(out_img);title('Closing of Opening');