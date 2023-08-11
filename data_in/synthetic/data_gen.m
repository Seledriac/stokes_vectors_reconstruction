clear all
close all
clc

% Génération de l'image polarimétrique synthétique

tab_sig=[0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
sig=tab_sig(5); % Intensité du bruit
sig = 0;

rand('state',sum(100*clock));
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];

mask = ones(256);
Ic = A*[1. 1/sqrt(2) 1/sqrt(2)]';
I1 = mask*Ic(1);%I0
I2 = mask*Ic(2);%I90
I3 = mask*Ic(3);%I45
I4 = mask*Ic(4);%I135

theta = (0:1/255:1)*pi/180;
epsilon = (0:45/255:45)*pi/180;
[T,E] = meshgrid(theta,epsilon);
Sv0 = ones(256);
Sv1 = cos(2*T).*cos(E);
Sv2 = sin(2*T).*cos(E);

% Sur le fond de l'image, le vecteur de Stokes en chaque pixel
% est [ 1, 1/sqrt(2), 1/sqrt(2) ] (polarisation complète)
real_S0=ones(256);
real_S1=(1/sqrt(2))*ones(256);
real_S2=(1/sqrt(2))*ones(256);


%Admissibility criterion

%Sv0.^2>=Sv1.^2+Sv2.^2

% Évolution verticale graduelle de la polarisation au sein du disque centré de rayon 100 pixels 
% haut = polarisation élevée, S1 haut
% bas = polarisation basse, S1 bas
% S2 : reste bas car sin(2T) est très faible car T est faible
for i = -128:128
    for j = -128:128
    	if sqrt(i*i+j*j) < 100
        S(1) = 1.;
        S(2) = Sv1(i+128,j+128);
        S(3) = Sv2(i+128,j+128);
        
        real_S0(128+i,128+j)=S(1);
        real_S1(128+i,128+j)=S(2);
        real_S2(128+i,128+j)=S(3);
        
        % Conversion des vecteurs de Stokes en image polarimétrique
        Ic = A*S';
        I1(128+i,128+j) = Ic(1);
        I2(128+i,128+j) = Ic(2);
        I3(128+i,128+j) = Ic(3);
        I4(128+i,128+j) = Ic(4);
	end      
    end
end


% figure 1, image polarimétrique de synthèse, canaux I0, I90, I45, I135
figure;subplot(221),imshow(I1,[]);colorbar 
subplot(222),imshow(I2,[]);colorbar
subplot(223),imshow(I3,[]);colorbar
subplot(224),imshow(I4,[]);colorbar
% 
% 
% figure;subplot(221),imagesc(I1);colorbar
% subplot(222),imagesc(I2);colorbar
% subplot(223),imagesc(I3);colorbar
% subplot(224),imagesc(I4);colorbar
% 
% figure;subplot(221),imshow(real_S0,[]);colorbar
% subplot(222),imshow(real_S1,[]);colorbar
% subplot(223),imshow(real_S2,[]);colorbar
% 

% I1, I2, I3, I4 non bruité
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

% figure 2, canaux I0, I90, I45, I135 non bruités vs bruités
figure;subplot(241), imshow(I1,[]);colorbar
subplot(242),imshow(I2,[]);colorbar
subplot(243),imshow(I3,[]);colorbar
subplot(244),imshow(I4,[]);colorbar

subplot(245),imshow(I0,[]);colorbar
subplot(246),imshow(I90,[]);colorbar
subplot(247),imshow(I45,[]);colorbar
subplot(248),imshow(I135,[]);colorbar

save('I0.mat','I0');
save('I90.mat','I90');
save('I45.mat','I45');
save('I135.mat','I135');

save('I0_clean.mat','I1');
save('I90_clean.mat','I2');
save('I45_clean.mat','I3');
save('I135_clean.mat','I4');

save('S0_synthetic_data.mat','real_S0');
save('S1_synthetic_data.mat','real_S1');
save('S2_synthetic_data.mat','real_S2');

print -dpng synthetic_data_radiance_images.png;