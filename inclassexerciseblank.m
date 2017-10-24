%% step 1: write a few lines of code or use FIJI to separately save the
% nuclear channel of the image Colony1.tif for segmentation in Ilastik


%% step 2: train a classifier on the nuclei
% try to get the get nuclei completely but separe them where you can
% save as both simple segmentation and probabilities

%% step 3: use h5read to read your Ilastik simple segmentation
% and display the binary masks produced by Ilastik 

% (datasetname = '/exported_data')
% Ilastik has the image transposed relative to matlab
% values are integers corresponding to segmentation classes you defined,
% figure out which value corresponds to nuclei

data=h5read('Segmentation.h5', '/exported_data');
data=squeeze(data);
figure(1);
imshow(data, []);

%% step 3.1: show segmentation as overlay on raw data

img=imread('48hColony1_DAPI.tif');
img_overlay=cat(3, imadjust(img), data, zeros(size(img)));
figure(2);
imshow(img_overlay, []);

%% step 4: visualize the connected components using label2rgb
% probably a lot of nuclei will be connected into large objects

BW=im2bw(data, graythresh(data));
CC=bwconncomp(BW);
L=labelmatrix(CC);
rgb=label2rgb(L,'jet','k');
figure(3);
imshow(rgb);

%% step 5: use h5read to read your Ilastik probabilities and visualize

% it will have a channel for each segmentation class you defined

prob=h5read('Prediction.h5','/exported_data');
prob=squeeze(prob);
figure(4);
imshow(prob,[]);

%% step 6: threshold probabilities to separate nuclei better

prob_mask=prob>quantile(quantile(prob, 0.75),0.75);
imshow(prob_mask);
BW2=im2bw(prob_mask);
CC2=bwconncomp(BW2);
L2=labelmatrix(CC2);
rgb2=label2rgb(L2, 'jet', 'k');
figure(5);
imshow(rgb2);
merge=imfuse(prob, BW2);
figure(6);
imshow(merge);

%% step 7: watershed to fill in the original segmentation (~hysteresis threshold)

img=imread('48hColony1_DAPI.tif');
mask = img > quantile(quantile(img,0.75),0.75);
figure(7);
imshow(mask, []);
CC=bwconncomp(mask);
stats=regionprops(CC, 'Area');
area=[stats.Area];
fused=area>mean(area)+std(area);
sublist=CC.PixelIdxList(fused);
sublist=cat(1, sublist{:});
fusedMask=false(size(mask));
fusedMask(sublist)=1;
figure(8);
imshow(fusedMask,'InitialMagnification','fit');

s=round(1.2*sqrt(mean(area))/pi);
nucmin=imerode(fusedMask, strel('disk',s));
figure(9);
imshow(nucmin,'InitialMagnification','fit');

out=~imdilate(fusedMask,strel('disk',1));
figure(10);
imshow(out,'InitialMagnification','fit');

basin=imcomplement(bwdist(out));
basin=imimposemin(basin,nucmin|out);
pcolor(basin);
shading flat;

L=watershed(basin);
figure(11);
imshow(L);
colormap('jet');
caxis([0 20]);

newNuclearMask=L>1|(mask-fusedMask);
figure(12);
imshow(newNuclearMask,'InitialMagnification','fit');

%% step 8: perform hysteresis thresholding in Ilastik and compare the results
% explain the differences

%% Ilastik is better at separating nuclei and better at edge detection.

%% step 9: clean up the results more if you have time 
% using bwmorph, imopen, imclose etc

newNuclearMask_open=imopen(newNuclearMask,strel('disk', 5));
figure(13);
imshow(newNuclearMask_open,[]);

newNuclearMask_morph=bwmorph(newNuclearMask,'remove');
figure(14);
imshow(newNuclearMask_morph);

newNuclearMask_close = imclose(newNuclearMask,strel('disk', 5));
figure(15);
imshow(newNuclearMask_open,[]);