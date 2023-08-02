clear all;
close all;
clc;


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
hands = true;
hazy_road = true;
key_ball_rubber = true;
mug = true;
mug_cafe = true;
pot = true;
road_experiment = true;
various_objects = true;

model='cm'; % modèle considéré


% ########## CALCUL #########

if refitting == true

    for i_mu_bi=6:6
        for i_mu_re=1:5

            mu_bi = mu_bi_r(i_mu_bi);
            mu_re = mu_re_r(i_mu_re);
            
            if hands == true
                data_source_in = '../../data_in/hands/';
                data_source_out = strcat('../../data_out/',model,'/hands/');
            end

            if hazy_road == true
                data_source_in = '../../data_in/hazy_road/';
                data_source_out = strcat('../../data_out/',model,'/hazy_road/');                
            end

            if key_ball_rubber == true
                data_source_in = '../../data_in/key_ball_rubber/';
                data_source_out = strcat('../../data_out/',model,'/key_ball_rubber/');                
            end

            if mug == true
                data_source_in = '../../data_in/mug/';
                data_source_out = strcat('../../data_out/',model,'/mug/');                
            end

            if mug_cafe == true
                data_source_in = '../../data_in/mug_cafe/';
                data_source_out = strcat('../../data_out/',model,'/mug_cafe/');                
            end

            if pot == true
                data_source_in = '../../data_in/pot/';
                data_source_out = strcat('../../data_out/',model,'/pot/');                
            end

            if road_experiment == true
                data_source_in = '../../data_in/road_experiment/';
                data_source_out = strcat('../../data_out/',model,'/road_experiment/');                
            end

            if various_objects == true
                data_source_in = '../../data_in/various_objects/';
                data_source_out = strcat('../../data_out/',model,'/various_objects/');                
            end
            
            metrics_real_data;

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

        if hands == true
            data_source_in = '../../data_in/hands/';
            data_source_out = strcat('../../data_out/',model,'/hands/');                
        end

        if hazy_road == true
            data_source_in = '../../data_in/hazy_road/';
            data_source_out = strcat('../../data_out/',model,'/hazy_road/');                
        end

        if key_ball_rubber == true
            data_source_in = '../../data_in/key_ball_rubber/';
            data_source_out = strcat('../../data_out/',model,'/key_ball_rubber/');                
        end

        if mug == true
            data_source_in = '../../data_in/mug/';
            data_source_out = strcat('../../data_out/',model,'/mug/');                
        end

        if mug_cafe == true
            data_source_in = '../../data_in/mug_cafe/';
            data_source_out = strcat('../../data_out/',model,'/mug_cafe/');                
        end

        if pot == true
            data_source_in = '../../data_in/pot/';
            data_source_out = strcat('../../data_out/',model,'/pot/');               
        end

        if road_experiment == true
            data_source_in = '../../data_in/road_experiment/';
            data_source_out = strcat('../../data_out/',model,'/road_experiment/');                
        end

        if various_objects == true
            data_source_in = '../../data_in/various_objects/';
            data_source_out = strcat('../../data_out/',model,'/various_objects/');                
        end

        metrics_real_data;

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
