% measure texture similarity in range of [0, 1] based on Gaussian derivative hostogram
% kambiz.rahbar@gmail.com, 2019.
%
% based on paper
% J.R.R.Uijlings, K.E.A.van de Sande, T.Gevers, A.W.M.Smeulders,"Selective
% Search for Object Recognition, " International Journal of Computer
% Vision, September 2013, Volume 104, Issue 2, pp 154–171
%

% clc
% clear
% close all
% 
% img1 = imread('texture1.jpg');
% img2 = imread('texture2.jpg');
% img3 = imread('texture3.jpg');
% img4 = imread('texture4.jpg');
% 
% s_texture1 = measure_texture_similarity(img1,img1);
% s_texture2 = measure_texture_similarity(img1,img2);
% s_texture3 = measure_texture_similarity(img1,img3);
% s_texture4 = measure_texture_similarity(img1,img4);
% 
% 
% % show the results
% figure(1);
% subplot(2,2,1);
% imshow(img1);
% title('image #1');
% 
% subplot(2,2,2); 
% imshow(img2);
% title('image #2');
% 
% subplot(2,2,3); 
% imshow(img3);
% title('image #3');
% 
% subplot(2,2,4); 
% imshow(img4);
% title('image #4');
% 
% 
% fprintf('texture similarity between image #1 and itself: %f\n',s_texture1);
% fprintf('texture similarity between image #1 and image #2: %f\n',s_texture2);
% fprintf('texture similarity between image #1 and image #3: %f\n',s_texture3);
% fprintf('texture similarity between image #1 and image #4: %f\n',s_texture4);

%%
function [s_texture] = measure_texture_similarity(img1, img2)
% measure texture similarity in range of [0, 1];

% calculate texture similarity vector for image #1 and #2
t1 = calculate_texture_similarity_vector(img1);
t2 = calculate_texture_similarity_vector(img2);

s_texture = sum(min(t1,t2))/24;
end


function [t] = calculate_texture_similarity_vector(img)

img = double(img);

% get RGB channels data
red_channel = img(:,:,1);
green_channel = img(:,:,2);
blue_channel = img(:,:,3);

% calculate texture similarity vector for each channel
t_red = calculate_channel_texture_vector(red_channel);
t_green = calculate_channel_texture_vector(green_channel);
t_blue = calculate_channel_texture_vector(blue_channel);

% aggrigate similarity vectors
t = [t_red; t_green; t_blue];

end


function [t] = calculate_channel_texture_vector(channel)
% define Gaussian derivatives in eight orientations
% using σ = 1 for each colour channel.
sigma = 1;
mask_size = 5;

h = fspecial('gauss',[round(mask_size*sigma), round(mask_size*sigma)], sigma);
[hx,hy] = gradient(h);

% apply Gaussian filter on the image
Hx = imfilter(channel,hx);
Hy = imfilter(channel,hy);

H_magnitude = sqrt(Hx.^2+Hy.^2);
H_orientation = atand(Hy./Hx)+90;
H_orientation(H_orientation==0) = eps;

% calculate bins histogram for each of eight orientations
channel_hist = zeros(8,10);
Theta = 0:180/8:180;
for k = 1:length(Theta)-1
    selected_magnitudes = H_magnitude(H_orientation > Theta(k) & H_orientation <= Theta(k+1));
    channel_hist(k,:) = histcounts(selected_magnitudes,10);
end

% normalize each histograms
normalized_channel_hist = channel_hist./sum(channel_hist,2);

% aggrigate the channel similarity vector
normalized_channel_hist = normalized_channel_hist';
t = normalized_channel_hist(:);
end

