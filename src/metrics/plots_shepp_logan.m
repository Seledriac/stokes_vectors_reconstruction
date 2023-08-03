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
PSNR_arr=zeros(length(mu_bi_r),length(mu_re_r));
PSNR_refitted_arr=zeros(length(mu_bi_r),length(mu_re_r));
SSIM_arr=zeros(length(mu_bi_r),length(mu_re_r));
SSIM_refitted_arr=zeros(length(mu_bi_r),length(mu_re_r));

% flags
refitting = true;

% Modèle considéré
model = 'cm';


% ########## CALCUL #########

if refitting == true

    for i_mu_bi=1:length(mu_bi_r)
        for i_mu_re=1:length(mu_re_r)

            mu_bi = mu_bi_r(i_mu_bi);
            mu_re = mu_re_r(i_mu_re);
            
            metrics_shepp_logan;

            PSNR_arr(i_mu_bi,i_mu_re)=PSNR;
            PSNR_refitted_arr(i_mu_bi,i_mu_re)=PSNR_refitted;

            SSIM_arr(i_mu_bi,i_mu_re)=mean(SSIM); % SSIM moyen sur les 4 channels
            SSIM_refitted_arr(i_mu_bi,i_mu_re)=mean(SSIM_refitted); % moyen sur les 4

        end
    end
    
    % surfaces
    figure;
    surf(mu_re_r,mu_bi_r,PSNR_arr);
    title('PSNR');
    xlabel('mu re');
    ylabel('mu bi');
    hold on;
    maxval = max(PSNR_arr(:));
    [row,col] = find(PSNR_arr==maxval);
    row = row(end);
    col = col(end);
    h = scatter3(mu_re_r(col),mu_bi_r(row),PSNR_arr(row,col),'filled');
    h.CData = [1 0 0];
    h.SizeData = 300;
    disp(strcat('max PSNR : mu_bi = ',sprintf('%.1f',mu_bi_r(row)),'; mu_re = ',sprintf('%.1f',mu_re_r(col))));

    figure;
    surf(mu_re_r,mu_bi_r,PSNR_refitted_arr);
    title('PSNR refitted');
    xlabel('mu re');
    ylabel('mu bi');
    hold on;
    maxval = max(PSNR_refitted_arr(:));
    [row,col] = find(PSNR_refitted_arr==maxval);
    row = row(end);
    col = col(end);
    h = scatter3(mu_re_r(col),mu_bi_r(row),PSNR_refitted_arr(row,col),'filled');
    h.CData = [1 0 0];
    h.SizeData = 300;
    disp(strcat('max PSNR refitted : mu_bi = ',sprintf('%.1f',mu_bi_r(row)),'; mu_re = ',sprintf('%.1f',mu_re_r(col))));
    
    figure;
    surf(mu_re_r,mu_bi_r,SSIM_arr);
    title('SSIM');
    xlabel('mu re');
    ylabel('mu bi');
    hold on;
    maxval = max(SSIM_arr(:));
    [row,col] = find(SSIM_arr==maxval);
    row = row(end);
    col = col(end);
    h = scatter3(mu_re_r(col),mu_bi_r(row),SSIM_arr(row,col),'filled');
    h.CData = [1 0 0];
    h.SizeData = 300;
    disp(strcat('max SSIM : mu_bi = ',sprintf('%.1f',mu_bi_r(row)),'; mu_re = ',sprintf('%.1f',mu_re_r(col))));
    
    figure;
    surf(mu_re_r,mu_bi_r,SSIM_refitted_arr);
    title('SSIM refitted');
    xlabel('mu re');
    ylabel('mu bi');
    hold on;
    maxval = max(SSIM_refitted_arr(:));
    [row,col] = find(SSIM_refitted_arr==maxval);
    row = row(end);
    col = col(end);
    h = scatter3(mu_re_r(col),mu_bi_r(row),SSIM_refitted_arr(row,col),'filled');
    h.CData = [1 0 0];
    h.SizeData = 300;
    disp(strcat('max SSIM refitted : mu_bi = ',sprintf('%.1f',mu_bi_r(row)),'; mu_re = ',sprintf('%.1f',mu_re_r(col))));

    % plots
    figure;
    plot(mu_bi_r,diag(PSNR_arr));
    [mxy,mxi]=max(diag(PSNR_arr));
    text(mu_bi_r(mxi),mxy,'max PSNR');
    hold on;
    plot(mu_bi_r(mxi),mxy,'*');
    hold on;
    plot(mu_bi_r,diag(PSNR_refitted_arr));
    [mxy,mxi]=max(diag(PSNR_refitted_arr));
    text(mu_bi_r(mxi),mxy,'max PSNR refitted');
    hold on;
    plot(mu_bi_r(mxi),mxy,'*');
    title('PSNR');
    xlabel('mu');
    legend('PSNR','PSNR refitted');
    
    figure;
    plot(mu_bi_r,diag(SSIM_arr));
    [mxy,mxi]=max(diag(SSIM_arr));
    text(mu_bi_r(mxi),mxy,'max SSIM');
    hold on;
    plot(mu_bi_r(mxi),mxy,'*');
    hold on;
    plot(mu_bi_r,diag(SSIM_refitted_arr));
    [mxy,mxi]=max(diag(SSIM_refitted_arr));
    text(mu_bi_r(mxi),mxy,'max SSIM refitted');
    hold on;
    plot(mu_bi_r(mxi),mxy,'*');
    title('SSIM');
    xlabel('mu');
    legend('SSIM','SSIM refitted');
    
else
    
    for i_mu_bi=1:length(mu_bi_r)

        mu_bi = mu_bi_r(i_mu_bi);

        metrics_shepp_logan;

        PSNR_arr(i_mu_bi,i_mu_bi)=PSNR;
        SSIM_arr(i_mu_bi,i_mu_bi)=mean(SSIM); % SSIM moyen sur les 4 channels
        
    end
    
    figure;
    plot(mu_bi_r,diag(PSNR_arr));
    [mxy,mxi]=max(diag(PSNR_arr));
    text(mu_bi_r(mxi),mxy,'max PSNR');
    hold on;
    plot(mu_bi_r(mxi),mxy,'*');    
    title('PSNR');
    xlabel('mu');
    legend('PSNR');
    
    figure;
    plot(mu_bi_r,diag(SSIM_arr));
    [mxy,mxi]=max(diag(SSIM_arr));
    text(mu_bi_r(mxi),mxy,'max SSIM');
    hold on;
    plot(mu_bi_r(mxi),mxy,'*');    
    title('SSIM');
    xlabel('mu');
    legend('SSIM');
        
end
