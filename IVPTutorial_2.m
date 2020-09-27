%NAME: Aryaman Sinha
%INST: IIT, Bhubaneswar
%DATE: 26.09.2020
%CATEGORY: BTech
%BRANCH: Electrical Engineering
%Roll Number: 17EE01016

% Image and Video Tutorial 02
clc;clear;close all;
%%
% Q1. Implementation of Marr-Hildreth Edge detector and Canny Edge
% detector on given input image.
img = imread('sample images/lena_gray_256.tif');
img = double(img)/255;
subplot(2,3,1);imshow(img);title('Sample Image');
% Marr-Hildreth Edge Detector
%Designing LOG & Gaussian Filter using sigma=2
sigma=2;
n=ceil(3*sigma)*2+1;
w = floor(n/2);
[x,y] = meshgrid(-w:w,-w:w);
g_filter = exp(-((x.^2+y.^2)/(2*sigma*sigma)));
g_filter = g_filter./sum(g_filter(:));
log_filter = ((x.^2+y.^2-2*sigma^2)/(sigma^4)).*g_filter;
subplot(2,3,2);surf(g_filter);title('Gaussain Filter(sigma=2)');
subplot(2,3,3);surf(-1*log_filter);title('LOG Filter(sigma=2)');
%using g_filter approach
img_pad = zeros(size(img,1)+2*w,size(img,2)+2*w);
img_pad(w+1:w+size(img,1),w+1:w+size(img,2))=img;
filter_out = zeros(size(img,1),size(img,2));
%zero padded image convolved to smooth out using low pass G_filter 
for i=w+1:size(img,1)-w
    for j=w+1:size(img,2)-w
        filter_out(i,j) = sum(g_filter.*img_pad(i-w:i+w,j-w:j+w),'all');
    end
end
%Laplacian for second order derivative computation 
laplacian_filter = [[-1,-1,-1];[-1,8,-1];[-1,-1,-1]];
lap_filt_out = zeros(size(img,1),size(img,2));
%finding second order derivative using  Laplacian filter
for i=2:size(img,1)-1
    for j=2:size(img,2)-1
        lap_filt_out(i,j) = sum(laplacian_filter.*filter_out(i-1:i+1,j-1:j+1),'all');
    end
end
lap_filt_out = lap_filt_out./max(lap_filt_out(:));
zero_crossing = zeros(size(img,1),size(img,2));
%changing 0->1 whenever zero crossing occurs iff absolute is >=threshold 
thresh = 0.0053;
for i=2:size(img,1)-1
    for j=2:size(img,2)-1
        if lap_filt_out(i,j)>0
            if (lap_filt_out(i,j-1)*lap_filt_out(i,j+1)<0)&&(abs(lap_filt_out(i,j-1))>=thresh)&&(abs(lap_filt_out(i,j+1))>=thresh)
                zero_crossing(i,j)=1;
            elseif (lap_filt_out(i-1,j)*lap_filt_out(i+1,j)<0)&&(abs(lap_filt_out(i-1,j))>=thresh)&&(abs(lap_filt_out(i+1,j))>=thresh)
                zero_crossing(i,j)=1;
            elseif (lap_filt_out(i-1,j-1)*lap_filt_out(i+1,j+1)<0)&&(abs(lap_filt_out(i-1,j-1))>=thresh)&&(abs(lap_filt_out(i+1,j+1))>=thresh)
                zero_crossing(i,j)=1;
            elseif (lap_filt_out(i+1,j-1)*lap_filt_out(i-1,j+1)<0)&&(abs(lap_filt_out(i+1,j-1))>=thresh)&&(abs(lap_filt_out(i-1,j+1))>=thresh)
                zero_crossing(i,j)=1;
            end
        end
    end
end
subplot(2,3,4);imshow(zero_crossing);title('Marr-Hildreth(Manual)');
%using direct log_filter approach
img_pad = zeros(size(img,1)+2*w,size(img,2)+2*w);
img_pad(w+1:w+size(img,1),w+1:w+size(img,2))=img;
filter_out = zeros(size(img,1),size(img,2));
%using LOG filter directly to convolove and smooth out the image
for i=w+1:size(img,1)-w
    for j=w+1:size(img,2)-w
        filter_out(i,j) = sum(log_filter.*img_pad(i-w:i+w,j-w:j+w),'all');
    end
end
%changing 0->1 whenever zero crossing occurs iff absolute is >=threshold 
zero_crossing = zeros(size(img,1),size(img,2));
thresh = 0.0053;
for i=2:size(img,1)-1
    for j=2:size(img,2)-1
        if filter_out(i,j)>0
            if (filter_out(i,j-1)*filter_out(i,j+1)<0)&&(abs(filter_out(i,j-1))>=thresh)&&(abs(filter_out(i,j+1))>=thresh)
                zero_crossing(i,j)=1;
            elseif (filter_out(i-1,j)*filter_out(i+1,j)<0)&&(abs(filter_out(i-1,j))>=thresh)&&(abs(filter_out(i+1,j))>=thresh)
                zero_crossing(i,j)=1;
            elseif (filter_out(i-1,j-1)*filter_out(i+1,j+1)<0)&&(abs(filter_out(i-1,j-1))>=thresh)&&(abs(filter_out(i+1,j+1))>=thresh)
                zero_crossing(i,j)=1;
            elseif (filter_out(i+1,j-1)*filter_out(i-1,j+1)<0)&&(abs(filter_out(i+1,j-1))>=thresh)&&(abs(filter_out(i-1,j+1))>=thresh)
                zero_crossing(i,j)=1;
            end
        end
    end
end
subplot(2,3,5);imshow(zero_crossing);title('Using LOG directly');
%using inbuilt marr-hildreth detector 
subplot(2,3,6);imshow(edge(img,'log'));title('Marr-Hildreth (Inbuilt)');figure;
% Canny Edge Detector
%Designing Gaussian filter with sigma=sqrt(2)
sigma=sqrt(2);
n=ceil(3*sigma)*2+1;
w = floor(n/2);
[x,y] = meshgrid(-w:w,-w:w);
g_filter = exp(-((x.^2+y.^2)/(2*sigma*sigma)));
g_filter = g_filter./sum(g_filter(:));
subplot(2,2,1);imshow(img);title('Sample Image');
subplot(2,2,2);surf(g_filter);title('Gaussain Filter(sigma=sqrt(2))');
img_pad = zeros(size(img,1)+2*w,size(img,2)+2*w);
img_pad(w+1:w+size(img,1),w+1:w+size(img,2))=img;
filter_out = zeros(size(img,1),size(img,2));
%getting smoothed out filter output 
for i=w+1:size(img,1)-w
    for j=w+1:size(img,2)-w
        filter_out(i,j) = sum(g_filter.*img_pad(i-w:i+w,j-w:j+w),'all');
    end
end
%using sobel filter to compute gradient values for x and y direction 
sobel_y = [[1,2,1];[0,0,0];[-1,-2,-1]];
sobel_x = [[-1,0,1];[-2,0,2];[-1,0,1]];
grad_x = zeros(size(img,1),size(img,2));
grad_y = zeros(size(img,1),size(img,2));
for i=2:size(img,1)-1
    for j=2:size(img,2)-1
        grad_x(i,j) = sum(sobel_x.*filter_out(i-1:i+1,j-1:j+1),'all');
        grad_y(i,j) = sum(sobel_y.*filter_out(i-1:i+1,j-1:j+1),'all');
    end
end
%mag and phase values for gradient 
mag_grad = sqrt(grad_x.^2+grad_y.^2);
mag_grad = mag_grad./max(mag_grad(:));
alpha_grad = atan2(grad_y,grad_x).*(180/pi);
idx = find(alpha_grad<0);
alpha_grad(idx)=alpha_grad(idx)+180;
%using non_maxima suppression here to check four directions respectively
%and get the Gnh and Gnl i.e. weak and strong edges of image Gnl will
%include Gnh. Threshold for each weak and strong edges accordingly for
%hysterisis thresholding operation 
grad_nh=zeros(size(img,1),size(img,2));
grad_nl=zeros(size(img,1),size(img,2));
thresh_h = 0.055;
thresh_l = 0.4*thresh_h;
for i=2:size(img,1)-1
    for j=2:size(img,2)-1
        if (67.5<=alpha_grad(i,j)&&alpha_grad(i,j)<112.5) %horizontal sector
            if mag_grad(i-1,j)<=mag_grad(i,j)&&mag_grad(i+1,j)<=mag_grad(i,j)
                if mag_grad(i,j)>=thresh_l
                    grad_nl(i,j)=1;
                    if mag_grad(i,j)>=thresh_h
                        grad_nh(i,j)=1;
                    end
                end
            end
        elseif (112.5<=alpha_grad(i,j)&&alpha_grad(i,j)<157.5)%+45 deg sector
            if mag_grad(i-1,j-1)<=mag_grad(i,j)&&mag_grad(i+1,j+1)<=mag_grad(i,j)
               if mag_grad(i,j)>=thresh_l
                    grad_nl(i,j)=1;
                    if mag_grad(i,j)>=thresh_h
                        grad_nh(i,j)=1;
                    end
                end
            end
        elseif (0<=alpha_grad(i,j)&&alpha_grad(i,j)<22.5)||(157.5<=alpha_grad(i,j)&&alpha_grad(i,j)<=180)%
            if mag_grad(i,j-1)<=mag_grad(i,j)&&mag_grad(i,j+1)<=mag_grad(i,j)
                if mag_grad(i,j)>=thresh_l
                    grad_nl(i,j)=1;
                    if mag_grad(i,j)>=thresh_h
                        grad_nh(i,j)=1;
                    end
                end
            end            
        elseif (22.5<=alpha_grad(i,j)&&alpha_grad(i,j)<67.5) %+45 deg sector
            if mag_grad(i-1,j+1)<=mag_grad(i,j)&&mag_grad(i+1,j-1)<=mag_grad(i,j)
               if mag_grad(i,j)>=thresh_l
                    grad_nl(i,j)=1;
                    if mag_grad(i,j)>=thresh_h
                        grad_nh(i,j)=1;
                    end
                end
            end
        end
    end
end
%as Gnl include Gnh we will 
grad_nl = grad_nl - grad_nh;
[r,c] = find(grad_nl>0);
%now applying 8-connectivity rule i.e. for remainaing weak edges is any 
%of the 8 neighbourhoods is strong then we will include that weak edge in
%strong edge.
for t=1:size(c)
    if (grad_nh(r(t)+1,c(t))>0||grad_nh(r(t)-1,c(t))>0||grad_nh(r(t),c(t)-1)>0||grad_nh(r(t),c(t)+1)>0||grad_nh(r(t)+1,c(t)+1)>0||grad_nh(r(t)-1,c(t)-1)>0||grad_nh(r(t)-1,c(t)+1)>0||grad_nh(r(t)+1,c(t)-1)>0)
        grad_nh(r(t),c(t))=1;
    end
end
subplot(2,2,3);imshow(grad_nh);title('Canny Edge(Manual)');
%using inbuilt canny edge detector 
subplot(2,2,4);imshow(edge(img,'canny'));title('Canny Edge(Inbuilt)');figure;
%%
% Q2. Peform phase only reconstruction using two images.
img1 = imread('sample images/cameraman.tif');
img2 = imread('sample images/mandril_gray.tif');
img1 = double(img1)/255;img2 = double(img2)/255;
n=size(img1,1);m=size(img1,2);
subplot(1,2,1);imshow(img1);title('Image-1');
subplot(1,2,2);imshow(img2);title('Image-2');figure;
%Taking DFT for each image, respective phase and mag spectrum
weights_row = exp(-1i*2*pi*(linspace(0,n-1,n))'*(linspace(0,n-1,n))/n);
weights_col = exp(-1i*2*pi*(linspace(0,m-1,m))'*(linspace(0,m-1,m))/m);
dft_img1 = (weights_row*(img1*weights_col));
phase_dft_img1 = angle(dft_img1);
mag_dft_img1 = abs(dft_img1);
dft_img2 = (weights_row*(img2*weights_col));
phase_dft_img2 = angle(dft_img2);
mag_dft_img2 = abs(dft_img2);
subplot(2,2,1);imshow(mag_dft_img1);title('Image-1 Mag Spectrum');
subplot(2,2,2);imshow(mag_dft_img2);title('Image-2 Mag Spectrum');
subplot(2,2,3);imshow(phase_dft_img1);title('Image-1 Phase Spectrum');
subplot(2,2,4);imshow(phase_dft_img2);title('Image-2 Phase Spectrum');figure;
phase_only_img1 = exp(1i*phase_dft_img1);%phase only for img-1
phase_only_img2 = exp(1i*phase_dft_img2);%phase only for img-2
dft_img21 = mag_dft_img2.*exp(1i*phase_dft_img1);%mixture mag(img2)<phase(img1)
dft_img12 = mag_dft_img1.*exp(1i*phase_dft_img2);%mixture mag(img1)<phase(img2)
%Doing inverse DFT i.e. x = (W^-1)X/N or (W*)X/N for 1-D 
%similarly for 2D x = (W_col*)((W_row*)X)/(N*M); 
idft_phase_only_img1 = real(conj(weights_col)*(phase_only_img1*conj(weights_row)))./(n*m);
idft_phase_only_img2 = real(conj(weights_col)*(phase_only_img2*conj(weights_row)))./(n*m);
idft_img21 = real(conj(weights_col)*(dft_img21*conj(weights_row)))./(n*m);
idft_img12 = real(conj(weights_col)*(dft_img12*conj(weights_row)))./(n*m);
subplot(2,2,1);imshow(idft_phase_only_img1);title('Img-1 Phase Only');
subplot(2,2,2);imshow(idft_phase_only_img2);title('Img-2 Phase Only');
subplot(2,2,3);imshow(idft_img21);title('Img-1 phase & Img-2 Mag');
subplot(2,2,4);imshow(idft_img12);title('Img-2 phase & Img-1 Mag');figure;
%%
% Q3. Compute 2D fourier spectrum of an image and center the magnitude
% spectrum and apply log transformation.
img = imread('sample images/cameraman.tif');
img = double(img)/255;
n=size(img,1);m=size(img,2);
subplot(1,3,1);imshow(img);title('Sample Image');
s1=floor(n/2);s2=floor(m/2);
%shifted DFT with N/2(in fourier domain) as we want to center the
%magnitude spectrum
weights_row = exp(-1i*2*pi*(linspace(0,n-1,n)-s1)'*(linspace(0,n-1,n)-s1)/n);
weights_col = exp(-1i*2*pi*(linspace(0,m-1,m)-s2)'*(linspace(0,m-1,m)-s2)/m);
dft_img = (weights_row*(img*weights_col));
mag_dft = abs(dft_img);
subplot(1,3,2);imshow(mag_dft);title('Center Shifted Mag.Spectrum');
log_tf_dft = log10(1+mag_dft);%log transformation
subplot(1,3,3);imshow(log_tf_dft);title('Log Transformed');
%%
% *Result & Discussion*:  In problem-1 we used the basic spatial filtering
% operation technique to convolove in time domain and marr-hildreth's
% algorithm to design manually all the steps each and every step is
% from scratch and also give us the almost same results as matlab's inbuilt
% result, here as we used zero padding to convolve , we can also use
% replicate padding for more better output, the threshold values are
% selected manually to have almost same result as inbuilt result with same
% defualt sigma values for filter as used by matlab inbuilt funtion.Same
% strategy was followed for canny edge detection using basic algorithm as
% discussed. Canny Edge detector gave more better results than
% Marr-Hildreth detector, more smooth and less false edges, as we see in
% the output. 
% In problem-2 we see the effect of reconstruction of image using fourier
% altered spectrum in one case phase only reconsturction was done where
% as in other case magnitude of one image was swapped with other. The
% observsation is usefull to note that even if we use other magnitude
% specturm the phase of the original image retains it features in time
% domain thus phase is very useful rather only using phase
% reconstruction won't give us much information back to time domain.
% In problem-3 we see the effect of shifting of fourier spectum in fourier
% domain and how log transformation can be useful to highlight it more
% significantly.