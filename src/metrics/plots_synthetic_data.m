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
M=size(I1,1);
N=size(I1,2);

% Paramètres
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];
% à faire varier
mu_bi_r=0.0:0.5:10.0; % poids du terme d'attache aux données dans la fonctionnelle du problème à solution biaisée
mu_re_r=0.0:0.5:10.0; % poids du terme d'attache aux données dans la fonctionnelle du problème de refitting

% metrics
error_SV_arr=zeros(length(mu_bi_r),length(mu_re_r));
error_SV_refitted_arr=zeros(length(mu_bi_r),length(mu_re_r));
PSNR_arr=zeros(length(mu_bi_r),length(mu_re_r));
PSNR_refitted_arr=zeros(length(mu_bi_r),length(mu_re_r));
SSIM_arr=zeros(length(mu_bi_r),length(mu_re_r));
SSIM_refitted_arr=zeros(length(mu_bi_r),length(mu_re_r));

% flags
refitting = true;

% Modèle considéré
model = 'ls';


% ########## CALCUL #########

if refitting == true

    if (strcmp(model,'ls'))
        for i_mu_re=1:length(mu_re_r)

            mu_re = mu_re_r(i_mu_re);
            
            metrics_synthetic_data;
            
            error_SV_arr(i_mu_re,i_mu_re)=error_SV;
            error_SV_refitted_arr(i_mu_re,i_mu_re)=error_SV_refitted;

            PSNR_arr(i_mu_re,i_mu_re)=PSNR;
            PSNR_refitted_arr(i_mu_re,i_mu_re)=PSNR_refitted;

            SSIM_arr(i_mu_re,i_mu_re)=mean(abs(SSIM)); % SSIM moyen sur les 4 channels
            SSIM_refitted_arr(i_mu_re,i_mu_re)=mean(abs(SSIM_refitted)); % moyen sur les 4

        end
    else
        for i_mu_bi=1:length(mu_bi_r)
            for i_mu_re=1:length(mu_re_r)
    
                mu_bi = mu_bi_r(i_mu_bi);
                mu_re = mu_re_r(i_mu_re);
                
                metrics_synthetic_data;
                
                error_SV_arr(i_mu_bi,i_mu_re)=error_SV;
                error_SV_refitted_arr(i_mu_bi,i_mu_re)=error_SV_refitted;
    
                PSNR_arr(i_mu_bi,i_mu_re)=PSNR;
                PSNR_refitted_arr(i_mu_bi,i_mu_re)=PSNR_refitted;
    
                SSIM_arr(i_mu_bi,i_mu_re)=mean(abs(SSIM)); % SSIM moyen sur les 4 channels
                SSIM_refitted_arr(i_mu_bi,i_mu_re)=mean(abs(SSIM_refitted)); % moyen sur les 4
    
            end
        end
    end
    
    if(~strcmp(model,'ls'))
        
        % surfaces
        f=figure();
        surf(mu_re_r,mu_bi_r,error_SV_arr);
        title('error SV','FontSize',24);
        xlabel('mu re','FontSize',24);
        ylabel('mu bi','FontSize',24);
        hold on;
        minval = min(error_SV_arr(:));
        [row,col] = find(error_SV_arr==minval);
        row = row(end);
        col = col(end);
        h = scatter3(mu_re_r(col),mu_bi_r(row),error_SV_arr(row,col),'filled');
        h.CData = [1 0 0];
        h.SizeData = 300;
        disp(strcat('min error SV : mu_bi = ',sprintf('%.1f',mu_bi_r(col)),'; mu_re = ',sprintf('%.1f',mu_re_r(row))));
        set(gca,"FontSize",24);
        set(gcf, 'Position', get(0, 'Screensize'));
        filename = fullfile('../../figures/',strcat(model,'_errorSV_surf_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
        exportgraphics(f,filename);
    
        f=figure();
        surf(mu_re_r,mu_bi_r,error_SV_refitted_arr);
        title('error SV refitted','FontSize',24);
        xlabel('mu re','FontSize',24);
        ylabel('mu bi','FontSize',24);
        hold on;
        minval = min(error_SV_refitted_arr(:));
        [row,col] = find(error_SV_refitted_arr==minval);
        row = row(end);
        col = col(end);
        h = scatter3(mu_re_r(col),mu_bi_r(row),error_SV_refitted_arr(row,col),'filled');
        h.CData = [1 0 0];
        h.SizeData = 300;
        disp(strcat('min error SV refitted : mu_bi = ',sprintf('%.1f',mu_bi_r(row)),'; mu_re = ',sprintf('%.1f',mu_re_r(col))));
        set(gca,"FontSize",24);
        set(gcf, 'Position', get(0, 'Screensize'));
        filename = fullfile('../../figures/',strcat(model,'_errorSV_r_surf_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
        exportgraphics(f,filename);
    
        f=figure();
        surf(mu_re_r,mu_bi_r,PSNR_arr);
        title('PSNR','FontSize',24);
        xlabel('mu re','FontSize',24);
        ylabel('mu bi','FontSize',24);
        hold on;
        maxval = max(PSNR_arr(:));
        [row,col] = find(PSNR_arr==maxval);
        row = row(end);
        col = col(end);
        h = scatter3(mu_re_r(col),mu_bi_r(row),PSNR_arr(row,col),'filled');
        h.CData = [1 0 0];
        h.SizeData = 300;
        disp(strcat('max PSNR : mu_bi = ',sprintf('%.1f',mu_bi_r(row)),'; mu_re = ',sprintf('%.1f',mu_re_r(col))));
        set(gca,"FontSize",24);
        set(gcf, 'Position', get(0, 'Screensize'));
        filename = fullfile('../../figures/',strcat(model,'_PSNR_surf_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
        exportgraphics(f,filename);
    
        f=figure();
        surf(mu_re_r,mu_bi_r,PSNR_refitted_arr);
        title('PSNR refitted','FontSize',24);
        xlabel('mu re','FontSize',24);
        ylabel('mu bi','FontSize',24);
        hold on;
        maxval = max(PSNR_refitted_arr(:));
        [row,col] = find(PSNR_refitted_arr==maxval);
        row = row(end);
        col = col(end);
        h = scatter3(mu_re_r(col),mu_bi_r(row),PSNR_refitted_arr(row,col),'filled');
        h.CData = [1 0 0];
        h.SizeData = 300;
        disp(strcat('max PSNR refitted : mu_bi = ',sprintf('%.1f',mu_bi_r(row)),'; mu_re = ',sprintf('%.1f',mu_re_r(col))));
        set(gca,"FontSize",24);
        set(gcf, 'Position', get(0, 'Screensize'));
        filename = fullfile('../../figures/',strcat(model,'_PSNR_r_surf_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
        exportgraphics(f,filename);
        
        f=figure();
        surf(mu_re_r,mu_bi_r,SSIM_arr);
        title('SSIM','FontSize',24);
        xlabel('mu re','FontSize',24);
        ylabel('mu bi','FontSize',24);
        hold on;
        maxval = max(SSIM_arr(:));
        [row,col] = find(SSIM_arr==maxval);
        row = row(end);
        col = col(end);
        h = scatter3(mu_re_r(col),mu_bi_r(row),SSIM_arr(row,col),'filled');
        h.CData = [1 0 0];
        h.SizeData = 300;
        disp(strcat('max SSIM : mu_bi = ',sprintf('%.1f',mu_bi_r(row)),'; mu_re = ',sprintf('%.1f',mu_re_r(col))));
        set(gca,"FontSize",24);
        set(gcf, 'Position', get(0, 'Screensize'));
        filename = fullfile('../../figures/',strcat(model,'_SSIM_surf_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
        exportgraphics(f,filename);
        
        f=figure();
        surf(mu_re_r,mu_bi_r,SSIM_refitted_arr);
        title('SSIM refitted','FontSize',24);
        xlabel('mu re','FontSize',24);
        ylabel('mu bi','FontSize',24);
        hold on;
        maxval = max(SSIM_refitted_arr(:));
        [row,col] = find(SSIM_refitted_arr==maxval);
        row = row(end);
        col = col(end);
        h = scatter3(mu_re_r(col),mu_bi_r(row),SSIM_refitted_arr(row,col),'filled');
        h.CData = [1 0 0];
        h.SizeData = 300;
        disp(strcat('max SSIM refitted : mu_bi = ',sprintf('%.1f',mu_bi_r(row)),'; mu_re = ',sprintf('%.1f',mu_re_r(col))));
        set(gca,"FontSize",24);
        set(gcf, 'Position', get(0, 'Screensize'));
        filename = fullfile('../../figures/',strcat(model,'_SSIM_r_surf_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
        exportgraphics(f,filename);
    end

    % plots
    f=figure();
    plot(mu_bi_r,diag(error_SV_arr),'LineWidth',3);
    [mny,mni]=min(diag(error_SV_arr));
    text(mu_bi_r(mni),mny,'min error SV','FontSize',24);
    disp(strcat('min error SV : mu = ',sprintf('%.1f',mu_bi_r(mni))));
    hold on;
    plot(mu_bi_r(mni),mny,'*','LineWidth',10);
    hold on;
    plot(mu_bi_r,diag(error_SV_refitted_arr),'LineWidth',3);
    [mny,mni]=min(diag(error_SV_refitted_arr));
    text(mu_bi_r(mni),mny,'min error SV refitted','FontSize',24);
    disp(strcat('min error SV refitted : mu = ',sprintf('%.1f',mu_bi_r(mni))));
    hold on;
    plot(mu_bi_r(mni),mny,'*','LineWidth',10);
    title('error SV','FontSize',24);
    xlabel('mu','FontSize',24);
    legend('error SV','error SV refitted','FontSize',24);
    set(gca,"FontSize",24);
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = fullfile('../../figures/',strcat(model,'_errorSV_plot_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
    exportgraphics(f,filename);
    

    f=figure();
    plot(mu_bi_r,diag(PSNR_arr),'LineWidth',3);
    [mxy,mxi]=max(diag(PSNR_arr));
    text(mu_bi_r(mxi),mxy,'max PSNR','FontSize',24);
    disp(strcat('max PSNR : mu = ',sprintf('%.1f',mu_bi_r(mxi))));
    hold on;
    plot(mu_bi_r(mxi),mxy,'*','LineWidth',10);
    xline(mu_bi_r(mxi),':','LineWidth',3);
    hold on;
    plot(mu_bi_r,diag(PSNR_refitted_arr),'LineWidth',3);
    [mxy,mxi]=max(diag(PSNR_refitted_arr));
    text(mu_bi_r(mxi),mxy,'max PSNR refitted','FontSize',24);
    disp(strcat('max PSNR refitted : mu = ',sprintf('%.1f',mu_bi_r(mxi))));
    hold on;
    plot(mu_bi_r(mxi),mxy,'*','LineWidth',10);
    xline(mu_bi_r(mxi),':','LineWidth',3);
    title('PSNR','FontSize',24);
    xlabel('mu','FontSize',24);
    legend('PSNR','PSNR refitted','FontSize',24);
    set(gca,"FontSize",24);
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = fullfile('../../figures/',strcat(model,'_PSNR_plot_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
    exportgraphics(f,filename);
    
    f=figure();
    plot(mu_bi_r,diag(SSIM_arr),'LineWidth',3);
    [mxy,mxi]=max(diag(SSIM_arr));
    text(mu_bi_r(mxi),mxy,'max SSIM','FontSize',24);
    disp(strcat('max SSIM : mu = ',sprintf('%.1f',mu_bi_r(mxi))));
    hold on;
    plot(mu_bi_r(mxi),mxy,'*','LineWidth',10);
    xline(mu_bi_r(mxi),':','LineWidth',3);
    hold on;
    plot(mu_bi_r,diag(SSIM_refitted_arr),'LineWidth',3);
    [mxy,mxi]=max(diag(SSIM_refitted_arr));
    text(mu_bi_r(mxi),mxy,'max SSIM refitted','FontSize',24);
    disp(strcat('max SSIM refitted : mu = ',sprintf('%.1f',mu_bi_r(mxi))));
    hold on;
    plot(mu_bi_r(mxi),mxy,'*','LineWidth',10);
    xline(mu_bi_r(mxi),':','LineWidth',3);
    title('SSIM','FontSize',24);
    xlabel('mu','FontSize',24);
    legend('SSIM','SSIM refitted','FontSize',24);
    set(gca,"FontSize",24);
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = fullfile('../../figures/',strcat(model,'_SSIM_plot_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
    exportgraphics(f,filename);
    
else
    if (strcmp(model,'ls'))
        metrics_synthetic_data;
        disp(strcat('PSNR = ', sprintf('%.3f',PSNR), '; SSIM = ', sprintf('%.3f',mean(abs(SSIM)))));
    else    
        for i_mu_bi=1:length(mu_bi_r)
    
            mu_bi = mu_bi_r(i_mu_bi);
    
            metrics_synthetic_data;
            
            error_SV_arr(i_mu_bi,i_mu_bi)=error_SV;
            PSNR_arr(i_mu_bi,i_mu_bi)=PSNR;
            SSIM_arr(i_mu_bi,i_mu_bi)=mean(abs(SSIM)); % SSIM moyen sur les 4 channels
            
        end
        
       % plots
        f=figure();
        plot(mu_bi_r,diag(error_SV_arr),'LineWidth',3);
        [mny,mni]=min(diag(error_SV_arr));
        text(mu_bi_r(mni),mny,'min error SV','FontSize',24);
        disp(strcat('min error SV : mu = ',sprintf('%.1f',mu_bi_r(mni))));
        hold on;
        plot(mu_bi_r(mni),mny,'*','LineWidth',10);
        title('error SV','FontSize',24);
        xlabel('mu','FontSize',24);
        legend('error SV','FontSize',24);
        set(gca,"FontSize",24);
        set(gcf, 'Position', get(0, 'Screensize'));
        filename = fullfile('../../figures/',strcat(model,'_errorSV_plot_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'.png'));
        exportgraphics(f,filename);
    
        f=figure();
        plot(mu_bi_r,diag(PSNR_arr),'LineWidth',3);
        [mxy,mxi]=max(diag(PSNR_arr));
        text(mu_bi_r(mxi),mxy,'Max PSNR','FontSize',24);
        disp(strcat('max PSNR : mu = ',sprintf('%.1f',mu_bi_r(mxi))));
        hold on;
        plot(mu_bi_r(mxi),mxy,'*','LineWidth',10);
        xline(mu_bi_r(mxi),':','LineWidth',3);
        title('PSNR','FontSize',24);
        xlabel('mu','FontSize',24);
        legend('PSNR','FontSize',24);
        set(gca,"FontSize",24);
        set(gcf, 'Position', get(0, 'Screensize'));
        filename = fullfile('../../figures/',strcat(model,'_PSNR_plot_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'.png'));
        exportgraphics(f,filename);
        
        f=figure();
        plot(mu_bi_r,diag(SSIM_arr),'LineWidth',3);
        [mxy,mxi]=max(diag(SSIM_arr));
        text(mu_bi_r(mxi),mxy,'Max SSIM','FontSize',24);
        disp(strcat('max SSIM : mu = ',sprintf('%.1f',mu_bi_r(mxi))));
        hold on;
        plot(mu_bi_r(mxi),mxy,'*','LineWidth',10);
        xline(mu_bi_r(mxi),':','LineWidth',3);
        title('SSIM','FontSize',24);
        xlabel('mu','FontSize',24);
        legend('SSIM','FontSize',24);
        set(gca,"FontSize",24);
        set(gcf, 'Position', get(0, 'Screensize'));
        filename = fullfile('../../figures/',strcat(model,'_SSIM_plot_synthetic_mu_bi_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'.png'));
        exportgraphics(f,filename);

    end
        
end
