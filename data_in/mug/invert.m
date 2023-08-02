clear all
close all
clc

I45=imread("I45.png");
I135=imread("I135.png");
I0=imread("I0.png");
I90=imread("I90.png");

I45=double(I45);
I45=I45/max(max(I45));
I135=double(I135);
I135=I135/max(max(I135));
I0=double(I0);
I0=I0/max(max(I0));
I90=double(I90);
I90=I90/max(max(I90));

M = size(I45,1);
N = size(I45,2);
I=zeros(M*2,N*2);
I(1:2:end,1:2:end)=I45;
I(2:2:end,2:2:end)=I135;
I(1:2:end,2:2:end)=I0;
I(2:2:end,1:2:end)=I90;

imwrite(I,'radiance_image.tiff','tiff');
