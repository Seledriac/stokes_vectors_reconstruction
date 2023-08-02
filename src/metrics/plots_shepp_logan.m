clear all;
close all;
clc;


% RÃ©cupÃ©ration du shepp logan phantom
fI1=matfile('../../data_in/shepp_logan_phantom/I0_clean.mat');
I1=fI1.I1;
fI2=matfile('../../data_in/shepp_logan_phantom/I90_clean.mat');
I2=fI2.I2;
fI3=matfile('../../data_in/shepp_logan_phantom/I45_clean.mat');
I3=fI3.I3;
fI4=matfile('../../data_in/shepp_logan_phantom/I135_clean.mat');
I4=fI4.I4;
M=size(I1,1);
N=size(I1,2);

% Paramètres
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];
mu_bi_r=1:0.1:2;
mu_re_r=1:0.1:2;

% metrics
PSNR_arr=zeros(5);
PSNR_refitted_arr=zeros(5);
SSIM_arr=zeros(5);
SSIM_refitted_arr=zeros(5);

% flags
refitting = true;

% Modèle considéré
model = 'cm';


% ########## CALCUL #########

if refitting == true

    for i_mu_bi=6:6
        for i_mu_re=1:5

            mu_bi = mu_bi_r(i_mu_bi);
            mu_re = mu_re_r(i_mu_re);
            
            metrics_shepp_logan;

            PSNR_arr(i_mu_bi,i_mu_re)=PSNR;
            PSNR_refitted_arr(i_mu_bi,i_mu_re)=PSNR_refitted;

            SSIM_arr(i_mu_bi,i_mu_re)=mean(SSIM); % SSIM moyen sur les 4 channels
            SSIM_refitted_arr(i_mu_bi,i_mu_re)=mean(SSIM_refitted); % moyen sur les 4

        end
    end
    
    % % surfaces
    % figure;
    % surf(mu_bi_r,mu_re_r,PSNR_arr);
    % figure;
    % surf(mu_bi_r,mu_re_r,PSNR_refitted_arr);
    % figure;
    % surf(mu_bi_r,mu_re_r,SSIM_arr);
    % figure;
    % surf(mu_bi_r,mu_re_r,SSIM_refitted_arr);

    % plots
    figure;
    title('PSNR');
    plot(mu_bi_r,diag(PSNR_arr));
    hold on;
    plot(mu_bi_r,diag(PSNR_refitted_arr));
    figure;
    title('SSIM');
    plot(mu_bi_r,diag(SSIM_arr));
    hold on;
    plot(mu_bi_r,diag(SSIM_refitted_arr));
    
else
    
    for i_mu_bi=1:length(mu_bi_r)

        mu_bi = mu_bi_r(i_mu_bi);

        metrics_shepp_logan;

        PSNR_arr(i_mu_bi,i_mu_bi)=PSNR;
        SSIM_arr(i_mu_bi,i_mu_bi)=mean(SSIM); % SSIM moyen sur les 4 channels
        
    end
    
    % plots
    figure;
    title('PSNR');
    plot(mu_bi_r,diag(PSNR_arr));
    figure;
    title('SSIM');
    plot(mu_bi_r,diag(SSIM_arr));
        
end
