%NAME: Aryaman Sinha
%INST: IIT, Bhubaneswar
%DATE: 19.09.2020
%CATEGORY: BTech
%BRANCH: Electrical Engineering
%Roll Number: 17EE01016

% Image and Video Tutorial 01
clc;clear;close all;
%%                               
% Q1.Given RGB Image separate intensity & color information using different 
% color models & display them separately & then combine them back into 
% original RGB and compare with original.
%Read the RGB Image (Mandril color as sample image)
img = imread('sample images/mandril_color.tif');
%Convert to [0,1] range
img = double(img)/255;
imshow(img);title('Sample RGB Image');figure;
%Fetching RGB planes separately for getting each planes intensity
%values/gray levels which can be used to find intensity info of image by
%avgering them
img_r = img(:,:,1);img_g=img(:,:,2);img_b=img(:,:,3);
I_img = (img_r+img_g+img_b)/3;
imshow(I_img);title('Intensity Info of Image');figure;
%Fetching Hue and Saturation info using HIS to RGB model universal conversion
H_img = acosd((((img_r-img_g)+(img_r-img_b))/2)./sqrt((img_r-img_g).^2+(img_r-img_b).*(img_g-img_b)));
H_img(img_b>img_g)=360-H_img(img_b>img_g);
H_img = H_img/360; %converting [0,1] angle convension 
S_img = 1-(1./I_img).*min(img,[],3);
%Adding Hue and Saturation infos to get chromatic info/color info of the image
chromac_img = cat(3,H_img,S_img,zeros(size(img,1)));
imshow(chromac_img);title('Chromatic Info of Image');figure;
hsi_img = cat(3,H_img,S_img,I_img);
imshow(hsi_img);title('HIS Model Representation of Image');figure;
H_img = 2*pi*H_img; %converting to [0,2pi] range for HIS to RGB conversion
rec_r = zeros(size(H_img));rec_g = zeros(size(H_img));rec_b = zeros(size(H_img));
%RG sector 
idx = find((0<=H_img)&(H_img<2*pi/3));
rec_b(idx) = I_img(idx).*(1-S_img(idx));
rec_r(idx) = I_img(idx).*(1+S_img(idx).*cos(H_img(idx))./cos(pi/3-H_img(idx)));
rec_g(idx) = 3*I_img(idx)-(rec_r(idx)+rec_b(idx));
%BG sector
idx = find((2*pi/3<=H_img)&(H_img<4*pi/3));
rec_r(idx) = I_img(idx).*(1-S_img(idx));
rec_g(idx) = I_img(idx).*(1+S_img(idx).*cos(H_img(idx)-2*pi/3)./cos(pi-H_img(idx)));
rec_b(idx) = 3*I_img(idx)-(rec_r(idx)+rec_g(idx));
%BR sector
idx = find((4*pi/3<=H_img)&(H_img<=2*pi));
rec_g(idx) = I_img(idx).*(1-S_img(idx));
rec_b(idx) = I_img(idx).*(1+S_img(idx).*cos(H_img(idx)-4*pi/3)./cos(5*pi/3-H_img(idx)));
rec_r(idx) = 3*I_img(idx)-(rec_g(idx)+rec_b(idx));
%Concating recovered RGB planes to get Recovered Image
rec_img = cat(3,rec_r,rec_g,rec_b);
imshow(rec_img);title('Recovered Image');figure;
%%                              
% Q2.Given a gray scale image find the image negative
img = imread('sample images/lena_gray_256.tif');%Lena Gray as sample image
imshow(img);title('Sample Gray Image');figure;
img_neg = 255-img;
imshow(img_neg);title('Negative of Image');figure;
%%                              
% Q3.Given a gray scale image find & display magnitude of fourier specturm and
% apply log transformation on it. 
img = imread('sample images/cameraman.tif');%Cameraman gray as sample image
n = size(img,1);m = size(img,2);
img = double(img)/255; %to [0,1] range
imshow(img);title('Sample Gray Image');figure;
%Using 2D-DFT i.e. first applying M-DFT to each one by one and then N-DFT 
%to result to get final DFT => DFT[x (i.e. NXM matrix)]  = N-DFT[M-DFT[x]] 
%or X(n,m)=sum{k,0...n-1}[sum{l,0...m-1} x(k,l)*exp(-2*i*pi*k*n/N)*exp(-2*i*pi*l*m/M) 
weights_row = exp(-1i*2*pi*linspace(0,n-1,n)'*linspace(0,n-1,n)/n);
weights_col = exp(-1i*2*pi*linspace(0,m-1,m)'*linspace(0,m-1,m)/m);
dft_img = (weights_row*(img*weights_col));
mag_dft_img = abs(dft_img);%magnitude of DFT(x)
imshow(mag_dft_img);title('DFT Magnitude of Image(Manual)');figure;
%Using Inbuilt func to verify results
imshow(abs(fft2(img)));title('DFT Magnitude of Image(using Inbuilt)');figure;
log_tf_dft = log10(1+mag_dft_img); %log transformation on mag of DFT
imshow(log_tf_dft);title('Log Transformation of DFT Mag.');figure;
%%                              
% Q4.Same as Problem-3 but apply Power law transformation
%Continuation from Problem-3
gamma_corr_dft = mag_dft_img.^0.4; %Gamma Correction to Image.
imshow(gamma_corr_dft);title('Power Transformation of DFT Mag.');figure;
%%                             
% Q5.Given  a gray scale image do a histogram equalization enhancement on it.
img = imread('sample images/peppers_gray.tif');%peppers_gray image as sample
img = img(:,:,1); %removing unncessary alpha channel plane as we don't need it here 
img_copy = img;
imshow(img);title('Sample Gray Image');figure;
%Histogram of image pixel values
histogram(img);title('Histogram of Image'); xlabel('Unique Pixel Values');
ylabel('Frequency of Pixel Values');figure;
vec_img = img(:); %converting image to 1-D vector 
%getting transformation results, get p(r), calculating num of times each
%pixel came / total num of pixels in image
arr = zeros(256,1);
for i=1:size(vec_img,1)
    arr(vec_img(i)+1)= arr(vec_img(i)+1)+1;
end
arr = arr/size(vec_img,1);
%p(s) = for i=1...N sum{p(r)for j<=i} (N=total number of unique pixel values i.e. 256 here) 
new_arr = zeros(256,1);
new_arr(1)=arr(1);
for i=2:256
    new_arr(i)=new_arr(i-1)+arr(i);
end
%converting to s pixel values using p(s) dist.
for i=1:size(img,1)
    for j=1:size(img,2)
        img(i,j)=round(new_arr(img(i,j)+1)*255);
    end
end
%histogram of transformed image
histogram(img);title('Histogram of Transformed Image'); xlabel('Unique Pixel Values');
ylabel('Frequency of Pixel Values');figure;
%transformed image
imshow(img);title('Transformed Image(Manual)');figure;
%using inbuilt histeq func to verify the results
imshow(histeq(img_copy));title('Transformed Image(using Inbuilt)');
%%
% *Results & Discussion*: As we can see in Problem-1 results the intensity
% information is nothing but the average of gray level of the resepctive
% channels which depicts the light intensity levels varying from white to
% gray to black color with light being the maximum intensity value. And as
% we know Hue and Saturation combination can tell us about the color
% information in the resultant output as we can see the level of redness,
% green color are dominant rather than blue color. RGB to HIS model
% conversion & RGB to HIS model formaulas were taken directly from the
% net and perfect recovery was seen. 
% For Problem-2 negative was easily fetched by L-1-r transformation as
% 8-bit so 255-r. 
% For Problem-3 and 4 using the 2D DFT formula we got the log and power
% transformation easly and seen that high value range were perfectly
% clipped back by log conversion do have valid gray level range. i.e. high
% values like ~10^6 -> ~6 and as for gamma converion with gamma = 2.5 i.e.
% power=0.4 we gamma correction and image FT was more light intensity near
% to white or near to 255 value. 
% For Problem-5 appyling the direct results we got in class i.e. taking
% running cummalative sum of p(r) to get p(s) to have almost uniform
% distribution. Results verifies with inbuilt function.