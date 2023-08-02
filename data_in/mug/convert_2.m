clear all
close all
clc

I45=imread("I45.png");
I135=imread("I135.png");
I0=imread("I0.png");
I90=imread("I90.png");

I45=double(I45);
I135=double(I135);
I0=double(I0);
I90=double(I90);

save('I0.mat','I0');
save('I90.mat','I90');
save('I45.mat','I45');
save('I135.mat','I135');