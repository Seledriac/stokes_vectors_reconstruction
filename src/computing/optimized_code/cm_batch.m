
clear all; % à éviter si on veut garder le pool
clc;


% ###### PARAMETRES #####

% fixés
Nb_iter_max=1000;
gamma_bi=1.0; % poids des TV dans la fonctionnelle du problème à solution biaisée, reste à 1
gamma_re=1.0; % poids de phi_SD dans la fonctionnelle du problème de refitting, reste à 1
sigma_bi=0.3; % paramètre CP itération biaisée, reste à 0.3
sigma_re=0.3; % paramètre CP itération refitting, reste à 0.3
tau_bi=1.0/(8.0*2*sigma_bi*gamma_bi^2); % tau_bi * sigma_bi * gamma_bi^2 = 1/16 < 1/8
tau_re=1.0/(8.0*2*sigma_re*gamma_re^2); % tau_re * sigma_re * gamma_re^2 = 1/16 < 1/8
theta_bi=1.0; % paramètre CP itération biaisée, reste à 1.0
theta_re=1.0; % paramètre CP itération refitting, reste à 1.0

% à faire varier
mu_bi_r=100.0:1.0:100.0; % poids du terme d'attache aux données dans la fonctionnelle du problème à solution biaisée
mu_re_r=100.0:1.0:100.0; % poids du terme d'attache aux données dans la fonctionnelle du problème de refitting

% flags
refitting = true;
hands = true;
hazy_road = false;
key_ball_rubber = false;
mug = false;
mug_cafe = false;
pot = true;
road_experiment = false;
shepp_logan_phantom = false;
synthetic = false;
various_objects = false;


% ########## CALCUL #########

if refitting == true

    for i_mu_bi=1:length(mu_bi_r)   
        for i_mu_re=1:length(mu_re_r)

            mu_bi = mu_bi_r(i_mu_bi);
            mu_re = mu_re_r(i_mu_re);
            
            if hands == true
                data_source = '../../../data_in/hands/';
                data_dest = '../../../data_out/cm/hands/';
                cm_refitted;
            end

            if hazy_road == true
                data_source = '../../../data_in/hazy_road/';
                data_dest = '../../../data_out/cm/hazy_road/';
                cm_refitted;
            end

            if key_ball_rubber == true
                data_source = '../../../data_in/key_ball_rubber/';
                data_dest = '../../../data_out/cm/key_ball_rubber/';
                cm_refitted;
            end

            if mug == true
                data_source = '../../../data_in/mug/';
                data_dest = '../../../data_out/cm/mug/';
                cm_refitted;
            end

            if mug_cafe == true
                data_source = '../../../data_in/mug_cafe/';
                data_dest = '../../../data_out/cm/mug_cafe/';
                cm_refitted;
            end

            if pot == true
                data_source = '../../../data_in/pot/';
                data_dest = '../../../data_out/cm/pot/';
                cm_refitted;
            end

            if road_experiment == true
                data_source = '../../../data_in/road_experiment/';
                data_dest = '../../../data_out/cm/road_experiment/';
                cm_refitted;
            end

            if shepp_logan_phantom == true
                data_source = '../../../data_in/shepp_logan_phantom/';
                data_dest = '../../../data_out/cm/shepp_logan_phantom/';
                cm_refitted;
            end

            if synthetic == true
                data_source = '../../../data_in/synthetic/';
                data_dest = '../../../data_out/cm/synthetic/';
                cm_refitted;
            end

            if various_objects == true
                data_source = '../../../data_in/various_objects/';
                data_dest = '../../../data_out/cm/various_objects/';
                cm_refitted;
            end

        end
    end
    
else
    
    for i_mu_bi=1:length(mu_bi_r)

        mu_bi = mu_bi_r(i_mu_bi);

        if hands == true
            data_source = '../../../data_in/hands/';
            data_dest = '../../../data_out/cm/hands/';
            cm_standard;
        end

        if hazy_road == true
            data_source = '../../../data_in/hazy_road/';
            data_dest = '../../../data_out/cm/hazy_road/';
            cm_standard;
        end

        if key_ball_rubber == true
            data_source = '../../../data_in/key_ball_rubber/';
            data_dest = '../../../data_out/cm/key_ball_rubber/';
            cm_standard;
        end

        if mug == true
            data_source = '../../../data_in/mug/';
            data_dest = '../../../data_out/cm/mug/';
            cm_standard;
        end

        if mug_cafe == true
            data_source = '../../../data_in/mug_cafe/';
            data_dest = '../../../data_out/cm/mug_cafe/';
            cm_standard;
        end

        if pot == true
            data_source = '../../../data_in/pot/';
            data_dest = '../../../data_out/cm/pot/';
            cm_standard;
        end

        if road_experiment == true
            data_source = '../../../data_in/road_experiment/';
            data_dest = '../../../data_out/cm/road_experiment/';
            cm_standard;
        end

        if shepp_logan_phantom == true
            data_source = '../../../data_in/shepp_logan_phantom/';
            data_dest = '../../../data_out/cm/shepp_logan_phantom/';
            cm_standard;
        end

        if synthetic == true
            data_source = '../../../data_in/synthetic/';
            data_dest = '../../../data_out/cm/synthetic/';
            cm_standard;
        end

        if various_objects == true
            data_source = '../../../data_in/various_objects/';
            data_dest = '../../../data_out/cm/various_objects/';
            cm_standard;
        end
        
    end
        
end
