%NAME: Aryaman Sinha
%INST: IIT, Bhubaneswar
%DATE: 26.10.2020
%CATEGORY: BTech
%BRANCH: Electrical Engineering
%Roll Number: 17EE01016

% Image and Video Tutorial 04
clc;clear;close all;
%%
% Q1. Take a fingerprint image, add salt and pepper noise. Perform the
% median filtering and contraharmonic mean filtering. Show the effect of the
% wrong choice of polarity in the order Q.
img = imread('sample images/fingerprint1.jpg');
img = double(img)/255;
img = img(:,:,1);
subplot(1,2,1);imshow(img);title('Sample Image');
Ps = 0.025; Pp = 0.025;%how much prob/perc pixel we want as salt / pepper noise 
noise_img = noise_add(img,Ps,Pp);
subplot(1,2,2);imshow(noise_img);title('Salt and Pepper Noise Image');figure;
%performing median filtering 
filt_img = mean_filt(noise_img);
subplot(1,2,1);imshow(noise_img);title('Noisy Image');
subplot(1,2,2);imshow(filt_img);title('Mean Filtered Image');figure;
%performing contraharmonic mean filtering
filt_img1 = contra_filt(noise_img,-2);
filt_img2 = contra_filt(noise_img,-1);
filt_img3 = contra_filt(noise_img,0);
filt_img4 = contra_filt(noise_img,1);
filt_img5 = contra_filt(noise_img,2);
subplot(2,3,1);imshow(noise_img);title('Noisy Image');
subplot(2,3,2);imshow(filt_img1);title('Contraharmonic Q=-2');
subplot(2,3,3);imshow(filt_img2);title('Contraharmonic Q=-1');
subplot(2,3,4);imshow(filt_img3);title('Contraharmonic Q=0');
subplot(2,3,5);imshow(filt_img4);title('Contraharmonic Q=1');
subplot(2,3,6);imshow(filt_img5);title('Contraharmonic Q=2');

function noisy_img = noise_add(img,Ps,Pp)
    P = Ps+Pp;
    n = size(img,1)*size(img,2);
    m = fix(P*n);
    idx = randperm(n,m); 
    k=fix((Ps/P)*m);
    idx1 = idx(1:k);
    idx2 = idx(k+1:end);
    idx1 = idx1' + n.*(0:size(img,3)-1);idx1 = idx1(:);
    idx2 = idx2' + n.*(0:size(img,3)-1);idx2 = idx2(:);
    img(idx1)=255;%salt noise
    img(idx2)=0;%pepper noise
    noisy_img = img;
end
function filt_img = mean_filt(img)
    filt_img = zeros(size(img));
    n=size(img,1);m=size(img,2);
    for i=2:n-1
        for j=2:m-1
            filt_img(i,j) = median(img(i-1:i+1,j-1:j+1),'all');
        end
    end
end
function filt_img = contra_filt(img,Q)
    filt_img = zeros(size(img));
    n=size(img,1);m=size(img,2);
    for i=2:n-1
        for j=2:m-1
            num = sum(img(i-1:i+1,j-1:j+1).^(Q+1),'all');
            den = sum(img(i-1:i+1,j-1:j+1).^Q,'all');
            filt_img(i,j)=num./den;
        end
    end
end
%%
% *Discussion*: Here we have simply added the salt and pepper noise
% according to the what much percentage pixels we want as noise, here
% divided the salt and pepper in 50:50 ratio. Then after adding noise
% simple median filtering is done which is simply checking N8 neighborhood
% and taking median in that space and replace center pixel by that, spatial
% filtering method. Then, we also saw the change of factor Q in
% contraharmoic mean filtering in which it was seen that the for Q=0 it is
% same as mean filtering and for Q=-1 its harmonic filteer. The Q>0 pepper
% noise gets eliminated whereas for Q<0 salt noise gets eliminated.
