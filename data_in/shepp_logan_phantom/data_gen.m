clear all
close all
clc

% Génération du fantôme de shepp-logan

tab_sig=[0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
sig=tab_sig(5); % Intensité du bruit


% I1, I2, I3, I4 non bruité
I1 = phantom('Modified Shepp-Logan',255);
I2 = phantom('Modified Shepp-Logan',255);
I3 = phantom('Modified Shepp-Logan',255);
I4 = phantom('Modified Shepp-Logan',255);
IF(:,:,1) = I1;
IF(:,:,2) = I2;
IF(:,:,3) = I3;
IF(:,:,4) = I4;

% IF = non bruité, IN = bruité
for i = 1:4
    IN(:,:,i) = IF(:,:,i) + sqrt(sig)*randn(size(IF(:,:,i)));
end

% I0, I90, I45, I135 bruité
I0=IN(:,:,1);
I90=IN(:,:,2);
I45=IN(:,:,3);
I135=IN(:,:,4);

save('I0.mat','I0');
save('I90.mat','I90');
save('I45.mat','I45');
save('I135.mat','I135');

save('I0_clean.mat','I1');
save('I90_clean.mat','I2');
save('I45_clean.mat','I3');
save('I135_clean.mat','I4');

% figure 1, canaux I0, I90, I45, I135 non bruités vs bruités
figure;subplot(241), imshow(I1,[]);colorbar
subplot(242),imshow(I2,[]);colorbar
subplot(243),imshow(I3,[]);colorbar
subplot(244),imshow(I4,[]);colorbar

subplot(245),imshow(I0,[]);colorbar
subplot(246),imshow(I90,[]);colorbar
subplot(247),imshow(I45,[]);colorbar
subplot(248),imshow(I135,[]);colorbar

print -dpng shepp_logan_phantom.png;
