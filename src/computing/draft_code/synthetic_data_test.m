%% 
clear all
close all
clc


% ######### Génération de l'image polarimétrique synthétique ########

tab_sig=[0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
sig=tab_sig(5); % Intensité du bruit

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

print -dpng synthetic_data_radiance_images.png;
%%


% ############################### CALCUL ###################################


% Mesure du temps de calcul
tic 

% Algorithme CP avec refitting pour reconstruction des vecteurs de Stokes
% à partir de l'image polarimétrique de synthèse bruitée

gamma_bi=1.0; % poids de TV dans la fonctionnelle du problème à solution biaisée
mu_bi=1.5; % poids du terme d'attache aux données dans la fonctionnelle du problème à solution biaisée

% sigma_bi * tau_bi * gamma_bi^2 = 1/16 < 1/8, convergence de CP pour problème à solution biaisée ok
sigma_bi=0.3;
tau_bi=1.0/(8.0*2*sigma_bi); 
theta_bi=1.0; % convergence ok

gamma_re=1.0; % poids de phi_SD dans la fonctionnelle du problème de refitting
mu_re=1.5; % poids du terme d'attache aux données dans la fonctionnelle du problème de refitting

% sigma_re * tau_re * gamma_re^2 = 1/16 < 1/8, convergence de CP pour problème de refitting ok
sigma_re=0.3;
tau_re=1.0/(8.0*2*sigma_re); 
theta_re=1.0; % convergence ok

M=size(I0,1);
N=size(I0,2);

% variables primales, problème à solution biaisée
S0_hat=zeros(M,N);
S1_hat=zeros(M,N);
S2_hat=zeros(M,N);
S0_new_hat=zeros(M,N);
S1_new_hat=zeros(M,N);
S2_new_hat=zeros(M,N);
S0_bar_hat=S0_hat;
S1_bar_hat=S1_hat;
S2_bar_hat=S2_hat;
S0_bar_new_hat=S0_hat;
S1_bar_new_hat=S1_hat;
S2_bar_new_hat=S2_hat;

% variables primales, problème de refitting
S0_til=zeros(M,N);
S1_til=zeros(M,N);
S2_til=zeros(M,N);
S0_new_til=zeros(M,N);
S1_new_til=zeros(M,N);
S2_new_til=zeros(M,N);
S0_bar_til=S0_til;
S1_bar_til=S1_til;
S2_bar_til=S2_til;

% variables duales, problème à solution biaisée
p1_S0_hat=zeros(M,N);
p2_S0_hat=zeros(M,N);
p1_S1_hat=zeros(M,N);
p2_S1_hat=zeros(M,N);
p1_S2_hat=zeros(M,N);
p2_S2_hat=zeros(M,N);

% variables duales, problème de refitting
p1_S0_til=zeros(M,N);
p2_S0_til=zeros(M,N);
p1_S1_til=zeros(M,N);
p2_S1_til=zeros(M,N);
p1_S2_til=zeros(M,N);
p2_S2_til=zeros(M,N);


Nb_iter_max=1000;
for k=1:Nb_iter_max


    % ################## ITÉRATION PROBLÈME A SOLUTION BIAISÉE ##################

    
    % ######### PREMIER PROX ##########

    % variables intermédiaires z_hat, argument du premier prox
    z1_S0_hat=p1_S0_hat+sigma_bi*gamma_bi*gradient_x(S0_bar_hat);
    z2_S0_hat=p2_S0_hat+sigma_bi*gamma_bi*gradient_y(S0_bar_hat);
    z1_S1_hat=p1_S1_hat+sigma_bi*gamma_bi*gradient_x(S1_bar_hat);
    z2_S1_hat=p2_S1_hat+sigma_bi*gamma_bi*gradient_y(S1_bar_hat);
    z1_S2_hat=p1_S2_hat+sigma_bi*gamma_bi*gradient_x(S2_bar_hat);
    z2_S2_hat=p2_S2_hat+sigma_bi*gamma_bi*gradient_y(S2_bar_hat);        

    % projection sur la boule de rayon 1 dans R^2 pour z_0_hat
    p1_S0_hat=z1_S0_hat./max(1,sqrt(z1_S0_hat.^2+z2_S0_hat.^2));
    p2_S0_hat=z2_S0_hat./max(1,sqrt(z1_S0_hat.^2+z2_S0_hat.^2)); 
    % projection sur la boule de rayon 1 dans R^4 pour z_tilde_hat
    p1_S1_hat=z1_S1_hat./max(1,sqrt(z1_S1_hat.^2+z2_S1_hat.^2+z1_S2_hat.^2+z2_S2_hat.^2));
    p2_S1_hat=z2_S1_hat./max(1,sqrt(z1_S1_hat.^2+z2_S1_hat.^2+z1_S2_hat.^2+z2_S2_hat.^2));   
    p1_S2_hat=z1_S2_hat./max(1,sqrt(z1_S1_hat.^2+z2_S1_hat.^2+z1_S2_hat.^2+z2_S2_hat.^2));
    p2_S2_hat=z2_S2_hat./max(1,sqrt(z1_S1_hat.^2+z2_S1_hat.^2+z1_S2_hat.^2+z2_S2_hat.^2));
   

    % ######### DEUXIÈME PROX ##########

    % argument du 2ème prox
    s0_arg_hat=S0_hat+tau_bi*gamma_bi*(divergence_x(p1_S0_hat)+divergence_y(p2_S0_hat));
    s1_arg_hat=S1_hat+tau_bi*gamma_bi*(divergence_x(p1_S1_hat)+divergence_y(p2_S1_hat));
    s2_arg_hat=S2_hat+tau_bi*gamma_bi*(divergence_x(p1_S2_hat)+divergence_y(p2_S2_hat));
    
    % calcul du prox séparé : pixel par pixel
    for i=1:M
        for j=1:N

            % variable intermédiaire b pour le pixel (i,j)
            b=A'*[I0(i,j);I90(i,j);I45(i,j);I135(i,j)];
            % argument s du prox pour le pixel (i,j)
            s=[s0_arg_hat(i,j);s1_arg_hat(i,j);s2_arg_hat(i,j)];
            % Calcul du deuxième prox en (i,j) :
            % fonction attache aux données + indicatrice de la contrainte
            % d'admissibilité physique
            sol_ex=constrained_optimization_part(b,s,tau_bi,mu_bi);
            
            S0_new_hat(i,j)=sol_ex(1);
            S1_new_hat(i,j)=sol_ex(2);
            S2_new_hat(i,j)=sol_ex(3);

        end
    end
    
    % Accélération de la convergence pour le problème à solution biaisée en continuant dans la direction
    % donnée par la mise à jour de S_hat
    S0_bar_new_hat=S0_new_hat+theta_bi*(S0_new_hat-S0_hat);  
    S1_bar_new_hat=S1_new_hat+theta_bi*(S1_new_hat-S1_hat);  
    S2_bar_new_hat=S2_new_hat+theta_bi*(S2_new_hat-S2_hat);


    % ################## ITÉRATION PROBLÈME DE REFITTING ##################

    
    % ######### CALCUL DU CO-SUPPORT ##########
   
    % En chaque pixel (i,j), le canal 1 est pour la variable duale z_0
    % et le canal 2 pour les variables duales (z_1,Z_2)
    I_hat=zeros(M,N,2);
    for i=1:M
        for j=1:N
            % on remplace gamma_bi par 1 dans inégalité car
            % la projection se fait sur la boule unité dans
            % l'itération du problème à solution biaisée pour la variable duale
            if(sqrt(z1_S0_hat(i,j)^2+z2_S0_hat(i,j)^2) > 1)
                I_hat(i,j,1) = 1;
            end
            if(sqrt(z1_S1_hat(i,j)^2+z2_S1_hat(i,j)^2+z1_S2_hat(i,j)^2+z2_S2_hat(i,j)^2) > 1)
                I_hat(i,j,2) = 1;
            end
        end
    end


    % ######### PREMIER PROX ##########

    % variables intermédiaires z_til, argument du premier prox
    z1_S0_til=p1_S0_til+sigma_re*gradient_x(S0_bar_til);
    z2_S0_til=p2_S0_til+sigma_re*gradient_y(S0_bar_til);
    z1_S1_til=p1_S1_til+sigma_re*gradient_x(S1_bar_til);
    z2_S1_til=p2_S1_til+sigma_re*gradient_y(S1_bar_til);
    z1_S2_til=p1_S2_til+sigma_re*gradient_x(S2_bar_til);
    z2_S2_til=p2_S2_til+sigma_re*gradient_y(S2_bar_til);

    % calcul du prox séparé : pixel par pixel
    for i=1:M
        for j=1:N

            % Canal z_0
            z_til_1 = [
                z1_S0_til(i,j)
                z2_S0_til(i,j)
            ];
            if(I_hat(i,j,1) == 1) % pixel du co-support de grad(S_0)
                % variable intermédiaire estimant grad(S_0)
                nu_1 = [
                    z1_S0_hat(i,j)
                    z2_S0_hat(i,j)
                ];
                norm_nu_1 = sqrt(z1_S0_hat(i,j)^2+z2_S0_hat(i,j)^2);
                % on remplace encore gamma_bi par 1 car
                % la projection se fait sur la boule unité dans
                % l'itération du problème à solution biaisée pour la variable duale
                psi_1 = ( ( norm_nu_1 - 1 ) / ( sigma_bi * norm_nu_1 ) ) * nu_1;
                % calcul de prox_{phi_SD^*}(z_til_1,psi_1,I_hat)
                sol_ex_1 = refitting_part(z_til_1,psi_1,gamma_re);
            else % pixel hors du co-support                
                sol_ex_1 = z_til_1;
            end
            p1_S0_til(i,j)=sol_ex_1(1);
            p2_S0_til(i,j)=sol_ex_1(2);

            % Canaux (z_1,z_2)
            z_til_2 = [
                z1_S1_til(i,j)
                z2_S1_til(i,j)
                z1_S2_til(i,j)
                z2_S2_til(i,j)
            ];            
            if(I_hat(i,j,2) == 1) % pixel du co-support de (grad(S_1),grad(S_2))
                % variable intermédiaire estimant (grad(S_1),grad(S_2))
                nu_2 = [
                    z1_S1_hat(i,j)
                    z2_S1_hat(i,j)
                    z1_S2_hat(i,j)
                    z2_S2_hat(i,j)
                ];
                norm_nu_2 = sqrt(z1_S1_hat(i,j)^2+z2_S1_hat(i,j)^2+z1_S2_hat(i,j)^2+z2_S2_hat(i,j)^2);
                % on remplace encore gamma_bi par 1 car
                % la projection se fait sur la boule unité dans
                % l'itération du problème à solution biaisée pour la variable duale
                psi_2 = ( ( norm_nu_2 - 1 ) / ( sigma_bi * norm_nu_2 ) ) * nu_2;
                % calcul de prox_{phi_SD^*}(z_til_2,psi_2,I_hat)
                sol_ex_2 = refitting_part(z_til_2,psi_2,gamma_re);                    
            else % pixel hors du co-support
                sol_ex_2 = z_til_2;
            end
            p1_S1_til(i,j)=sol_ex_2(1);
            p2_S1_til(i,j)=sol_ex_2(2);
            p1_S2_til(i,j)=sol_ex_2(3);
            p2_S2_til(i,j)=sol_ex_2(4);
        end
    end


    % ######### DEUXIÈME PROX ##########

    % argument du 2ème prox
    s0_arg_til=S0_til+tau_re*(divergence_x(p1_S0_til)+divergence_y(p2_S0_til));
    s1_arg_til=S1_til+tau_re*(divergence_x(p1_S1_til)+divergence_y(p2_S1_til));
    s2_arg_til=S2_til+tau_re*(divergence_x(p1_S2_til)+divergence_y(p2_S2_til));

    % calcul du prox séparé : pixel par pixel
    for i=1:M
        for j=1:N            
            % variable intermédiaire b pour le pixel (i,j)
            b=A'*[I0(i,j);I90(i,j);I45(i,j);I135(i,j)];
            % argument du prox s pour le pixel (i,j)
            s=[s0_arg_til(i,j);s1_arg_til(i,j);s2_arg_til(i,j)];
            % Calcul du deuxième prox en (i,j) :
            % fonction attache aux données + indicatrice de la contrainte
            % d'admissibilité physique
            sol_ex=constrained_optimization_part(b,s,tau_re,mu_re);
            
            S0_new_til(i,j)=sol_ex(1);
            S1_new_til(i,j)=sol_ex(2);
            S2_new_til(i,j)=sol_ex(3);
        end
    end
    
    % Accélération de la convergence pour le problème de refitting en continuant dans la direction
    % donnée par la mise à jour de S_til
    S0_bar_til=S0_new_til+theta_re*(S0_new_til-S0_til);  
    S1_bar_til=S1_new_til+theta_re*(S1_new_til-S1_til);  
    S2_bar_til=S2_new_til+theta_re*(S2_new_til-S2_til);
   
    S0_hat=S0_new_hat;
    S1_hat=S1_new_hat;
    S2_hat=S2_new_hat;
    S0_bar_hat=S0_bar_new_hat;
    S1_bar_hat=S1_bar_new_hat;
    S2_bar_hat=S2_bar_new_hat;
    S0_til=S0_new_til;
    S1_til=S1_new_til;
    S2_til=S2_new_til;
    
end

% Solution biaisée de l'algorithme CP appliqué à
% l'image polarimétrique d'entrée I
S0_KKT_coupled_regularization=S0_hat;
S1_KKT_coupled_regularization=S1_hat;
S2_KKT_coupled_regularization=S2_hat;

% Solution débiaisée de l'algorithme CP appliqué à
% l'image polarimétrique d'entrée I
S0_KKT_coupled_regularization_refitted=S0_til;
S1_KKT_coupled_regularization_refitted=S1_til;
S2_KKT_coupled_regularization_refitted=S2_til;

% Mesure du temps de calcul
toc 

%%

% ############################### PLOTTTING ###################################


% figure 3, S reconstruit biaisé VS S reconstruit débiaisé VS S vérité
% terrain, S0 S1 S2, 3x3
c_min=min([min(min(S0_KKT_coupled_regularization)), min(min(S1_KKT_coupled_regularization)), min(min(S2_KKT_coupled_regularization)), min(min(S0_KKT_coupled_regularization_refitted)), min(min(S1_KKT_coupled_regularization_refitted)), min(min(S2_KKT_coupled_regularization_refitted)), min(min(real_S0)), min(min(real_S1)), min(min(real_S2))]);
c_max=max([max(max(S0_KKT_coupled_regularization)), max(max(S1_KKT_coupled_regularization)), max(max(S2_KKT_coupled_regularization)), max(max(S0_KKT_coupled_regularization_refitted)), max(max(S1_KKT_coupled_regularization_refitted)), max(max(S2_KKT_coupled_regularization_refitted)), max(max(real_S0)), max(max(real_S1)), max(max(real_S2))]);
figure;
subplot(311),imshow([S0_KKT_coupled_regularization S1_KKT_coupled_regularization S2_KKT_coupled_regularization],[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(312),imshow([S0_KKT_coupled_regularization_refitted S1_KKT_coupled_regularization_refitted S2_KKT_coupled_regularization_refitted],[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(313),imshow([real_S0 real_S1 real_S2],[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];

% % figure 4, S1 reconstruit biaisé VS S1 reconstruit débiaisé VS S1 vérité terrain
% c_min=min([min(min(S1_KKT_coupled_regularization)), min(min(S1_KKT_coupled_regularization_refitted)), min(min(real_S1))]);
% c_max=max([max(max(S1_KKT_coupled_regularization)), max(max(S1_KKT_coupled_regularization_refitted)), max(max(real_S1))]);
% figure;
% subplot(131),imshow(S1_KKT_coupled_regularization,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(132),imshow(S1_KKT_coupled_regularization_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(133),imshow(real_S1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];

print -dpng synthetic_data_S1_reconstruction_comparison.png;

% Calcul de la métrique de précision : distance entre
% S vérité terrain VS S reconstruit biaisé
error_SV=0.0;
for i=1:M
    for j=1:N
        Sgt=[real_S0(i,j);real_S1(i,j);real_S2(i,j)];
        S_star=[S0_KKT_coupled_regularization(i,j);S1_KKT_coupled_regularization(i,j);S2_KKT_coupled_regularization(i,j)];
        error_SV=error_SV+(norm(Sgt-S_star,2)*norm(Sgt-S_star,2))/(norm(Sgt,2)*norm(Sgt,2));
        
    end
end
error_SV=100*sqrt(error_SV*(1/(M*N)))

% Calcul de la métrique de précision : distance entre
% S vérité terrain VS S reconstruit débiaisé
error_SV_refitted=0.0;
for i=1:M
    for j=1:N
        Sgt=[real_S0(i,j);real_S1(i,j);real_S2(i,j)];
        S_star_refitted=[S0_KKT_coupled_regularization_refitted(i,j);S1_KKT_coupled_regularization_refitted(i,j);S2_KKT_coupled_regularization_refitted(i,j)];
        error_SV_refitted=error_SV_refitted+(norm(Sgt-S_star_refitted,2)*norm(Sgt-S_star_refitted,2))/(norm(Sgt,2)*norm(Sgt,2));
        
    end
end
error_SV_refitted=100*sqrt(error_SV_refitted*(1/(M*N)))

% % figure 5, S vérité terrain VS S reconstruit biaisé, canaux séparés
% c_min=min([min(min(real_S0)), min(min(real_S1)), min(min(real_S2)), min(min(S0_KKT_coupled_regularization)), min(min(S1_KKT_coupled_regularization)), min(min(S2_KKT_coupled_regularization))]);
% c_max=max([max(max(real_S0)), max(max(real_S1)), max(max(real_S2)), max(max(S0_KKT_coupled_regularization)), max(max(S1_KKT_coupled_regularization)), max(max(S2_KKT_coupled_regularization))]);
% figure;
% subplot(231),imshow(real_S0,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(232),imshow(real_S1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(233),imshow(real_S2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(234),imshow(S0_KKT_coupled_regularization,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(235),imshow(S1_KKT_coupled_regularization,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(236),imshow(S2_KKT_coupled_regularization,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% 
% print -dpng synthetic_data_stokes_reconstruction.png;

% % figure 6, S vérité terrain VS S reconstruit débiaisé, canaux séparés
% c_min=min([min(min(real_S0)), min(min(real_S1)), min(min(real_S2)), min(min(S0_KKT_coupled_regularization_refitted)), min(min(S1_KKT_coupled_regularization_refitted)), min(min(S2_KKT_coupled_regularization_refitted))]);
% c_max=max([max(max(real_S0)), max(max(real_S1)), max(max(real_S2)), max(max(S0_KKT_coupled_regularization_refitted)), max(max(S1_KKT_coupled_regularization_refitted)), max(max(S2_KKT_coupled_regularization_refitted))]);
% figure;
% subplot(231),imshow(real_S0,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(232),imshow(real_S1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(233),imshow(real_S2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(234),imshow(S0_KKT_coupled_regularization_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(235),imshow(S1_KKT_coupled_regularization_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(236),imshow(S2_KKT_coupled_regularization_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% 
% print -dpng synthetic_data_stokes_reconstruction_refitted.png;

% % figure 7, S reconstruit biaisé VS S reconstruit débiaisé, canaux séparés
% c_min=min([min(min(S0_KKT_coupled_regularization)), min(min(S1_KKT_coupled_regularization)), min(min(S2_KKT_coupled_regularization)), min(min(S0_KKT_coupled_regularization_refitted)), min(min(S1_KKT_coupled_regularization_refitted)), min(min(S2_KKT_coupled_regularization_refitted))]);
% c_max=max([max(max(S0_KKT_coupled_regularization)), max(max(S1_KKT_coupled_regularization)), max(max(S2_KKT_coupled_regularization)), max(max(S0_KKT_coupled_regularization_refitted)), max(max(S1_KKT_coupled_regularization_refitted)), max(max(S2_KKT_coupled_regularization_refitted))]);
% figure;
% subplot(231),imshow(S0_KKT_coupled_regularization,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(232),imshow(S1_KKT_coupled_regularization,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(233),imshow(S2_KKT_coupled_regularization,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(234),imshow(S0_KKT_coupled_regularization_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(235),imshow(S1_KKT_coupled_regularization_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(236),imshow(S2_KKT_coupled_regularization_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% 
% print -dpng comparison_synthetic_data_stokes_reconstruction_refitting.png;


m=1000;
cm_inferno=inferno(m);
colormap(inferno);

% % figure 8, DOLP S vérité terrain VS DOLP S reconstruit biaisé VS DOLP S reconstruit débiaisé, imshow
% c_min=min([min(min(sqrt(real_S1.^2+real_S2.^2)./real_S0)), min(min(sqrt(S1_KKT_coupled_regularization.^2+S2_KKT_coupled_regularization.^2)./S0_KKT_coupled_regularization)), min(min(sqrt(S1_KKT_coupled_regularization_refitted.^2+S2_KKT_coupled_regularization_refitted.^2)./S0_KKT_coupled_regularization_refitted))]);
% c_max=max([max(max(sqrt(real_S1.^2+real_S2.^2)./real_S0)), max(max(sqrt(S1_KKT_coupled_regularization.^2+S2_KKT_coupled_regularization.^2)./S0_KKT_coupled_regularization)), max(max(sqrt(S1_KKT_coupled_regularization_refitted.^2+S2_KKT_coupled_regularization_refitted.^2)./S0_KKT_coupled_regularization_refitted))]);
% figure;
% subplot(131),imshow(sqrt(real_S1.^2+real_S2.^2)./real_S0,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(132),imshow(sqrt(S1_KKT_coupled_regularization.^2+S2_KKT_coupled_regularization.^2)./S0_KKT_coupled_regularization,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(132),imshow(sqrt(S1_KKT_coupled_regularization_refitted.^2+S2_KKT_coupled_regularization_refitted.^2)./S0_KKT_coupled_regularization_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];

% figure 9, DOLP S vérité terrain VS DOLP S reconstruit biaisé VS DOLP S
% reconstruit débiaisé, imagesc
figure;
subplot(131),imagesc(sqrt(real_S1.^2+real_S2.^2)./real_S0,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(132),imagesc(sqrt(S1_KKT_coupled_regularization.^2+S2_KKT_coupled_regularization.^2)./S0_KKT_coupled_regularization,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(133),imagesc(sqrt(S1_KKT_coupled_regularization_refitted.^2+S2_KKT_coupled_regularization_refitted.^2)./S0_KKT_coupled_regularization_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];

% figure 10, AOLP S vérité terrain VS AOLP S reconstruit biaisé VS AOLP S
% reconstruit débiaisé, imagesc
c_min=min([min(min(0.5*atan2(real_S2,real_S1))), min(min(0.5*atan2(S2_KKT_coupled_regularization,S1_KKT_coupled_regularization))), min(min(0.5*atan2(S2_KKT_coupled_regularization_refitted,S1_KKT_coupled_regularization_refitted)))]);
c_max=max([max(max(0.5*atan2(real_S2,real_S1))), max(max(0.5*atan2(S2_KKT_coupled_regularization,S1_KKT_coupled_regularization))), max(max(0.5*atan2(S2_KKT_coupled_regularization_refitted,S1_KKT_coupled_regularization_refitted)))]);
figure;
subplot(131),imagesc(0.5*atan2(real_S2,real_S1),[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(132),imagesc(0.5*atan2(S2_KKT_coupled_regularization,S1_KKT_coupled_regularization),[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(133),imagesc(0.5*atan2(S2_KKT_coupled_regularization_refitted,S1_KKT_coupled_regularization_refitted),[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];

% Reconstruction de l'image polarimétrique de synthèse à partir du
% S reconstruit biaisé
for i=1:M
    for j=1:N
        I_recovered=A*[S0_KKT_coupled_regularization(i,j);S1_KKT_coupled_regularization(i,j);S2_KKT_coupled_regularization(i,j)];
        I0_recovered(i,j)=I_recovered(1);
        I90_recovered(i,j)=I_recovered(2);
        I45_recovered(i,j)=I_recovered(3);
        I135_recovered(i,j)=I_recovered(4);
    end
end

% % figure 11, image polarimétrique : non bruitée VS bruitée VS reconstruite
% % biaisée
% c_min=min([min(min(I1)), min(min(I2)), min(min(I3)), min(min(I4)), min(min(I0)), min(min(I90)), min(min(I45)), min(min(I135)), min(min(I0_recovered)), min(min(I90_recovered)), min(min(I45_recovered)), min(min(I135_recovered))]);
% c_max=max([max(max(I1)), max(max(I2)), max(max(I3)), max(max(I4)), max(max(I0)), max(max(I90)), max(max(I45)), max(max(I135)), max(max(I0_recovered)), max(max(I90_recovered)), max(max(I45_recovered)), max(max(I135_recovered))]);
% figure;
% subplot(341), imshow(I1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(342),imshow(I2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(343),imshow(I3,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(344),imshow(I4,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(345),imshow(I0,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(346),imshow(I90,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(347),imshow(I45,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(348),imshow(I135,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(349),imshow(I0_recovered,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(3,4,10),imshow(I90_recovered,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(3,4,11),imshow(I45_recovered,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(3,4,12),imshow(I135_recovered,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];

% PSNR entre image polarimétrique non bruitée et reconstruite biaisée
PSNR=10*log10(1/((1/(4*M*N))*(sum(sum((I1-I0_recovered).^2))+sum(sum((I2-I90_recovered).^2))+...
    sum(sum((I3-I45_recovered).^2))+sum(sum((I4-I135_recovered).^2)))))

% Reconstruction de l'image polarimétrique de synthèse à partir du
% S reconstruit débiaisé
for i=1:M
    for j=1:N
        I_recovered_refitted=A*[S0_KKT_coupled_regularization_refitted(i,j);S1_KKT_coupled_regularization_refitted(i,j);S2_KKT_coupled_regularization_refitted(i,j)];
        I0_recovered_refitted(i,j)=I_recovered_refitted(1);
        I90_recovered_refitted(i,j)=I_recovered_refitted(2);
        I45_recovered_refitted(i,j)=I_recovered_refitted(3);
        I135_recovered_refitted(i,j)=I_recovered_refitted(4);
    end
end

% % figure 12, image polarimétrique : non bruitée VS bruitée VS reconstruite
% % débiaisée
% c_min=min([min(min(I1)), min(min(I2)), min(min(I3)), min(min(I4)), min(min(I0)), min(min(I90)), min(min(I45)), min(min(I135)), min(min(I0_recovered_refitted)), min(min(I90_recovered_refitted)), min(min(I45_recovered_refitted)), min(min(I135_recovered_refitted))]);
% c_max=max([max(max(I1)), max(max(I2)), max(max(I3)), max(max(I4)), max(max(I0)), max(max(I90)), max(max(I45)), max(max(I135)), max(max(I0_recovered_refitted)), max(max(I90_recovered_refitted)), max(max(I45_recovered_refitted)), max(max(I135_recovered_refitted))]);
% figure;
% subplot(341), imshow(I1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(342),imshow(I2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(343),imshow(I3,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(344),imshow(I4,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(345),imshow(I0,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(346),imshow(I90,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(347),imshow(I45,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(348),imshow(I135,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(349),imshow(I0_recovered_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(3,4,10),imshow(I90_recovered_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(3,4,11),imshow(I45_recovered_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
% subplot(3,4,12),imshow(I135_recovered_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];

% PSNR entre image polarimétrique non bruitée et reconstruite débiaisée
PSNR_refitted=10*log10(1/((1/(4*M*N))*(sum(sum((I1-I0_recovered_refitted).^2))+sum(sum((I2-I90_recovered_refitted).^2))+...
    sum(sum((I3-I45_recovered_refitted).^2))+sum(sum((I4-I135_recovered_refitted).^2)))))


% figure 13, image polarimétrique : non bruitée VS bruitée VS reconstruite
% biaisée VS reconstruite débiaisée, 4x4
c_min=min([min(min(I1)), min(min(I2)), min(min(I3)), min(min(I4)), min(min(I0)), min(min(I90)), min(min(I45)), min(min(I135)), min(min(I0_recovered)), min(min(I90_recovered)), min(min(I45_recovered)), min(min(I135_recovered)), min(min(I0_recovered_refitted)), min(min(I90_recovered_refitted)), min(min(I45_recovered_refitted)), min(min(I135_recovered_refitted))]);
c_max=max([max(max(I1)), max(max(I2)), max(max(I3)), max(max(I4)), max(max(I0)), max(max(I90)), max(max(I45)), max(max(I135)), max(max(I0_recovered)), max(max(I90_recovered)), max(max(I45_recovered)), max(max(I135_recovered)), max(max(I0_recovered_refitted)), max(max(I90_recovered_refitted)), max(max(I45_recovered_refitted)), max(max(I135_recovered_refitted))]);
figure;
subplot(441), imshow(I1,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(442),imshow(I2,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(443),imshow(I3,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(444),imshow(I4,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(445),imshow(I0,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(446),imshow(I90,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(447),imshow(I45,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(448),imshow(I135,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(449),imshow(I0_recovered,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(4,4,10),imshow(I90_recovered,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(4,4,11),imshow(I45_recovered,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(4,4,12),imshow(I135_recovered,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(4,4,13),imshow(I0_recovered_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(4,4,14),imshow(I90_recovered_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(4,4,15),imshow(I45_recovered_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];
subplot(4,4,16),imshow(I135_recovered_refitted,[c_min,c_max]);colorbar;axis off;caxis=[c_min,c_max];


% for i=1:256
%   
%     for j=1:256
%         
%         b=A\[I1(i,j);I2(i,j);I3(i,j);I4(i,j)];
%         
%         check_S0(i,j)=b(1);
%         check_S1(i,j)=b(2);
%         check_S2(i,j)=b(3);
%     end
% end
% 
% figure;subplot(221),imshow(check_S0,[]);colorbar
% subplot(222),imshow(check_S1,[]);colorbar
% subplot(223),imshow(check_S2,[]);colorbar

