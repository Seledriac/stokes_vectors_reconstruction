clear all
close all
clc
I=imread('radiance_image.bmp');

I=double(I);
I=I/max(max(I));

I45=I(1:2:end,1:2:end);
I135=I(2:2:end,2:2:end);
I0=I(1:2:end,2:2:end);
I90=I(2:2:end,1:2:end);

save('I0.mat','I0');
save('I90.mat','I90');
save('I45.mat','I45');
save('I135.mat','I135');

figure;imshow(I0,[]);colorbar
print -dpng I0.png;
figure;imshow(I90,[]);colorbar
print -dpng I90.png;
figure;imshow(I45,[]);colorbar
print -dpng I45.png;
figure;imshow(I135,[]);colorbar
print -dpng I135.png;