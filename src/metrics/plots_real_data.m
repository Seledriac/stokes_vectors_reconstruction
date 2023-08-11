% clear all;
% close all;
% clc;


% Paramètres
A=(1/2)*[1 1 0;1 -1 0;1 0 1;1 0 -1];
mu_bi_r=99.0:1.0:100.0;
mu_re_r=1:0.1:2;

% metrics
PSNR_arr=zeros(length(mu_bi_r),length(mu_re_r));
PSNR_refitted_arr=zeros(length(mu_bi_r),length(mu_re_r));
SSIM_arr=zeros(length(mu_bi_r),length(mu_re_r));
SSIM_refitted_arr=zeros(length(mu_bi_r),length(mu_re_r));

% flags
refitting = false;
% image polarimétrique, en choisir une à la fois
hands = false;
hazy_road = false;
key_ball_rubber = false;
mug = false;
mug_cafe = false;
pot = false;
road_experiment = false;
various_objects = true;

model='ls';

if hands == true
    data_source_in = '../../data_in/hands/';
    data_source_out = strcat('../../data_out/',model,'/hands/');
    img_name='hands';
end

if hazy_road == true
    data_source_in = '../../data_in/hazy_road/';
    data_source_out = strcat('../../data_out/',model,'/hazy_road/');
    img_name='hazy_road';
end

if key_ball_rubber == true
    data_source_in = '../../data_in/key_ball_rubber/';
    data_source_out = strcat('../../data_out/',model,'/key_ball_rubber/');                
    img_name='key_ball_rubber';
end

if mug == true
    data_source_in = '../../data_in/mug/';
    data_source_out = strcat('../../data_out/',model,'/mug/');       
    img_name='mug';
end

if mug_cafe == true
    data_source_in = '../../data_in/mug_cafe/';
    data_source_out = strcat('../../data_out/',model,'/mug_cafe/');                
    img_name='mug_cafe';
end

if pot == true
    data_source_in = '../../data_in/pot/';
    data_source_out = strcat('../../data_out/',model,'/pot/');                
    img_name='pot';
end

if road_experiment == true
    data_source_in = '../../data_in/road_experiment/';
    data_source_out = strcat('../../data_out/',model,'/road_experiment/');                
    img_name='road_experiment';
end

if various_objects == true
    data_source_in = '../../data_in/various_objects/';
    data_source_out = strcat('../../data_out/',model,'/various_objects/');                
    img_name='various_objects';
end

% model='cm'; % modèle considéré


% ########## CALCUL #########

if refitting == true

    if (strcmp(model,'ls'))
        for i_mu_re=1:length(mu_re_r)

            mu_re = mu_re_r(i_mu_re);
            
            metrics_real_data;

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
                
                metrics_real_data;
    
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
        filename = fullfile('../../figures/',strcat(model,'_PSNR_surf_',img_name,'_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
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
        filename = fullfile('../../figures/',strcat(model,'_PSNR_r_surf_',img_name,'_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
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
        filename = fullfile('../../figures/',strcat(model,'_SSIM_surf_',img_name,'_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
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
        filename = fullfile('../../figures/',strcat(model,'_SSIM_r_surf_',img_name,'_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
        exportgraphics(f,filename);
    
    end

    % plots    

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
    filename = fullfile('../../figures/',strcat(model,'_PSNR_plot_',img_name,'_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
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
    filename = fullfile('../../figures/',strcat(model,'_SSIM_plot_',img_name,'_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'_mu_re_',sprintf('%.1f',mu_re_r(1)),'_',sprintf('%.1f',mu_re_r(end)),'.png'));
    exportgraphics(f,filename);
    
else
    if (strcmp(model,'ls'))
        metrics_real_data;
        disp(strcat('PSNR = ', sprintf('%.3f',PSNR), '; SSIM = ', sprintf('%.3f',mean(abs(SSIM)))));
    else    
        for i_mu_bi=1:length(mu_bi_r)
    
            mu_bi = mu_bi_r(i_mu_bi);
    
            metrics_real_data;
            
            PSNR_arr(i_mu_bi,i_mu_bi)=PSNR;
            SSIM_arr(i_mu_bi,i_mu_bi)=mean(abs(SSIM)); % SSIM moyen sur les 4 channels
            
        end
        
       % plots
    
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
        filename = fullfile('../../figures/',strcat(model,'_PSNR_plot_',img_name,'_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'.png'));
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
        filename = fullfile('../../figures/',strcat(model,'_SSIM_plot_',img_name,'_',sprintf('%.1f',mu_bi_r(1)),'_',sprintf('%.1f',mu_bi_r(end)),'.png'));
        exportgraphics(f,filename);

    end
        
end