clear all;
close all;
clc;


% RÃ©cupÃ©ration du synthetic data
freal_S0=matfile('../../data_in/synthetic/S0_synthetic_data.mat');
real_S0=freal_S0.real_S0;
freal_S1=matfile('../../data_in/synthetic/S1_synthetic_data.mat');
real_S1=freal_S1.real_S1;
freal_S2=matfile('../../data_in/synthetic/S2_synthetic_data.mat');
real_S2=freal_S2.real_S2;
fI1=matfile('../../data_in/synthetic/I0_clean.mat');
I1=fI1.I1;
fI2=matfile('../../data_in/synthetic/I90_clean.mat');
I2=fI2.I2;
fI3=matfile('../../data_in/synthetic/I45_clean.mat');
I3=fI3.I3;
fI4=matfile('../../data_in/synthetic/I135_clean.mat');
I4=fI4.I4;
M=size(I0,1);
N=size(I0,2);

% Paramètres
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];
mu_bi_r=1:0.1:2;
mu_re_r=1:0.1:2;

% metrics
error_SV_arr=zeros(5);
error_SV_refitted_arr=zeros(5);
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
            
            metrics_synthetic_data;
            
            error_SV_arr(i_mu_bi,i_mu_re)=error_SV;
            error_SV_refitted_arr(i_mu_bi,i_mu_re)=error_SV_refitted;

            PSNR_arr(i_mu_bi,i_mu_re)=PSNR;
            PSNR_refitted_arr(i_mu_bi,i_mu_re)=PSNR_refitted;

            SSIM_arr(i_mu_bi,i_mu_re)=mean(SSIM); % SSIM moyen sur les 4 channels
            SSIM_refitted_arr(i_mu_bi,i_mu_re)=mean(SSIM_refitted); % moyen sur les 4

        end
    end
    
    % % surfaces
    % figure;
    % surf(mu_bi_r,mu_re_r,error_SV_arr);
    % figure;
    % surf(mu_bi_r,mu_re_r,error_SV_refitted_arr);
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
    title('error_SV');
    plot(mu_bi_r,diag(error_SV_arr));
    hold on;
    plot(mu_bi_r,diag(error_SV_refitted_arr));
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

        metrics_synthetic_data;
        
        error_SV_arr(i_mu_bi,i_mu_re)=error_SV;
        PSNR_arr(i_mu_bi,i_mu_bi)=PSNR;
        SSIM_arr(i_mu_bi,i_mu_bi)=mean(SSIM); % SSIM moyen sur les 4 channels
        
    end
    
    % plots
    figure;
    title('error_SV');
    plot(mu_bi_r,diag(error_SV_arr));
    figure;
    title('PSNR');
    plot(mu_bi_r,diag(PSNR_arr));
    figure;
    title('SSIM');
    plot(mu_bi_r,diag(SSIM_arr));
        
end
